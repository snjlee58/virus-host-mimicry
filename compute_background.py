#!/usr/bin/env python3
"""Count proteins in the Foldseek target DB per species of interest.

Produces the background (N, K_species) for the hypergeometric taxonomic
enrichment test (Lasso et al. 2021, Figure 1B).

Reads the Foldseek DB's `_mapping` file (one line per DB entry: internal_id\\ttaxid)
and walks the NCBI taxonomy tree from each entry's taxid up to root, marking
the entry as belonging to a species if that species' taxid appears in the
ancestor chain.

With --divisions, instead counts entries per host taxonomic division (Fig 1C)
using the same bucket logic as classify_hits.py (bacteria/plantfungi/invertebrate/
vertebrate/human, plus virus/phage/unknown). This is the background for the 5x5
Fig 1C enrichment matrix.

This is a one-time computation per target DB. The output TSV is then consumed by
the enrichment step. IMPORTANT: --db must match the DB searched in run_foldseek.sh,
otherwise N/K do not correspond to the hits being tested.

Example:
    # Fig 1B (species: ecoli/yeast/human)
    python compute_background.py
        -> tests/ecoli/background_counts.tsv

    # Fig 1C (the five host divisions)
    python compute_background.py --divisions
        -> tests/host_divisions/background_counts.tsv

    python compute_background.py --divisions --db /path/to/foldseek_db -o bg.tsv
"""
import argparse
from collections import defaultdict
from pathlib import Path
import sys

from classify_hits import load_nodes, classify

# Species taxids we test in Figure 1B.
# Counting at species level means a strain-level taxid (e.g. E. coli K-12 = 83333)
# still counts toward "E. coli" because 562 appears in its ancestor chain.
SPECIES = {
    "ecoli": 562,    # Escherichia coli
    "yeast": 4932,   # Saccharomyces cerevisiae
    "human": 9606,   # Homo sapiens
}

# Output order for --divisions mode (Fig 1C). The first five are the matrix
# rows; virus/phage/unknown are part of the DB total N but are not Fig 1C rows.
BUCKET_ORDER = [
    "bacteria", "plantfungi", "invertebrate", "vertebrate", "human",
    "virus", "phage", "unknown",
]

# Foldseek target DB that searches run against (DEDUPLICATED, LIKE PDB100)
# TARGET_DB, or N/K won't correspond to the hits being tested).
DEFAULT_DB = "/fast2/yewon1/AFCDB_analysis_data/foldseek_search_PDBe/foldseek_pdb_db/gpu_pdb"


def _split_dmp(line):
    return line.rstrip("\n").rstrip("|").rstrip("\t").split("\t|\t")


def load_parents(path):
    """nodes.dmp -> {taxid: parent_taxid}"""
    parents = {}
    with open(path) as f:
        for line in f:
            fields = _split_dmp(line)
            parents[int(fields[0])] = int(fields[1])
    return parents


def load_merged(path):
    """merged.dmp -> {old_taxid: new_taxid}"""
    merged = {}
    with open(path) as f:
        for line in f:
            fields = _split_dmp(line)
            merged[int(fields[0])] = int(fields[1])
    return merged


def species_of(taxid, parents, merged, species_roots, cache):
    """Return the species key (e.g. 'ecoli') for taxid, or None."""
    if taxid in cache:
        return cache[taxid]
    resolved = merged.get(taxid, taxid)
    current = resolved
    seen = set()
    while current and current not in seen:
        for sp_name, sp_root in species_roots.items():
            if current == sp_root:
                cache[taxid] = sp_name
                return sp_name
        seen.add(current)
        parent = parents.get(current)
        if parent is None or parent == current:
            break
        current = parent
    cache[taxid] = None
    return None


# ---------------------------------------------------------------------------
# TODO (revisit for a faithful Fig 1C): two known divergences from Lasso et al.
#
# 1. GRANULARITY — chain-level, not UniProt-level.
#    This counts raw Foldseek DB entries (PDB chains, e.g. "3syy-assembly1_A").
#    The paper maps PDB -> UniProt via SIFTS and counts UNIQUE PROTEINS, "to
#    minimize experimental bias of multiple PDB entries for the same protein"
#    (their human background is 5,841 UniProt proteins). Heavily-crystallized
#    organisms (human, E. coli) are over-counted here non-uniformly, which biases
#    the hypergeometric (makes them look LESS enriched). To make faithful: map
#    gpu_pdb targets -> UniProt (SIFTS) and dedup, here AND in the neighbor space
#    (keep best E-value per UniProt). Until then this is a chain-level approximation.
#
# 2. VERTEBRATE NESTING — handled downstream, not here (intentional).
#    The paper's "vertebrate" category INCLUDES human (vertebrate = human +
#    nonHumanMammal + nonMammalVertebrate). We emit DISJOINT buckets as building
#    blocks; the enrichment step must form the Fig 1C vertebrate row/background as
#    (vertebrate + human). Do NOT "fix" that by merging here.
# ---------------------------------------------------------------------------
def count_divisions(mapping_path, nodes, merged):
    """Count Foldseek DB entries per host taxonomic division bucket.

    Returns (N, counts, no_taxid): N = total DB entries, counts = {bucket: n},
    no_taxid = entries with taxid 0/missing. Buckets follow classify_hits.classify()
    (bacteria, plantfungi, invertebrate, vertebrate, human, virus, phage, unknown),
    so the division definitions match exactly what classify_hits.py assigns to the
    hits — same DIV_TO_BUCKET mapping and the same human-9606 carve-out.

    NOTE: counts are at PDB-chain level, not UniProt level — see TODO above.
    """
    counts = defaultdict(int)
    N = 0
    no_taxid = 0
    with open(mapping_path) as f:
        for line in f:
            # _mapping is whitespace-separated (space, not tab) in current foldseek
            parts = line.split()
            if len(parts) < 2:
                continue
            N += 1
            taxid_str = parts[1]
            if taxid_str in ("", "0"):
                no_taxid += 1
            counts[classify(taxid_str, nodes, merged)] += 1
    return N, counts, no_taxid


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--divisions", action="store_true",
        help="Count DB entries per host taxonomic division (Fig 1C), using the same "
             "bucket logic as classify_hits.py, instead of per species (Fig 1B).",
    )
    ap.add_argument(
        "--db", default=DEFAULT_DB,
        help="Foldseek target DB path prefix; MUST match the DB searched in "
             "run_foldseek.sh (default: %(default)s)",
    )
    ap.add_argument(
        "--taxdump", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/taxdump"),
        help="Directory containing nodes.dmp and merged.dmp",
    )
    ap.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output TSV. Default: tests/host_divisions/background_counts.tsv "
             "(--divisions) or tests/ecoli/background_counts.tsv (species mode).",
    )
    args = ap.parse_args()

    if args.output is None:
        base = Path("/fast/sunny/virus-host-mimicry/tests")
        args.output = (base / "host_divisions/background_counts.tsv" if args.divisions
                       else base / "ecoli/background_counts.tsv")

    mapping_path = Path(args.db + "_mapping")
    nodes_path = args.taxdump / "nodes.dmp"
    merged_path = args.taxdump / "merged.dmp"

    for p in (mapping_path, nodes_path, merged_path):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")

    print(f"Loading {merged_path}...", file=sys.stderr)
    merged = load_merged(merged_path)
    print(f"  {len(merged):,} merged taxid redirects", file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.divisions:
        print(f"Loading {nodes_path} (division map)...", file=sys.stderr)
        nodes = load_nodes(nodes_path)
        print(f"  {len(nodes):,} taxid entries", file=sys.stderr)

        print(f"Iterating {mapping_path}...", file=sys.stderr)
        N, counts, no_taxid = count_divisions(mapping_path, nodes, merged)

        print(f"\nN (total DB entries):   {N:,}", file=sys.stderr)
        print(f"  no/invalid taxid:     {no_taxid:,}", file=sys.stderr)
        for b in BUCKET_ORDER:
            tag = "  <- Fig 1C row" if b in BUCKET_ORDER[:5] else ""
            print(f"  {b:>12}: {counts.get(b, 0):>9,}{tag}", file=sys.stderr)

        with open(args.output, "w") as f:
            f.write("label\tcount\n")
            f.write(f"_total\t{N}\n")
            for b in BUCKET_ORDER:
                f.write(f"{b}\t{counts.get(b, 0)}\n")
        print(f"\nWrote {args.output}", file=sys.stderr)
        return

    # --- Species mode (Fig 1B; e.g. the E. coli reproduction) ---
    print(f"Loading {nodes_path} (parent map)...", file=sys.stderr)
    parents = load_parents(nodes_path)
    print(f"  {len(parents):,} taxid entries", file=sys.stderr)

    print(f"Iterating {mapping_path}...", file=sys.stderr)
    cache = {}
    counts = defaultdict(int)
    N = 0
    no_taxid = 0

    with open(mapping_path) as f:
        for line in f:
            # _mapping is whitespace-separated (space, not tab) in current foldseek
            parts = line.split()
            if len(parts) < 2:
                continue
            N += 1
            try:
                taxid = int(parts[1])
            except ValueError:
                no_taxid += 1
                continue
            if taxid == 0:
                no_taxid += 1
                continue
            sp = species_of(taxid, parents, merged, SPECIES, cache)
            if sp is not None:
                counts[sp] += 1

    print(f"\nN (total DB entries):   {N:,}", file=sys.stderr)
    print(f"  no/invalid taxid:     {no_taxid:,}", file=sys.stderr)
    for sp in SPECIES:
        print(f"  {sp:>6}: {counts[sp]:>8,}", file=sys.stderr)

    with open(args.output, "w") as f:
        f.write("species\tcount\n")
        f.write(f"_total\t{N}\n")
        for sp in SPECIES:
            f.write(f"{sp}\t{counts[sp]}\n")
    print(f"\nWrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
