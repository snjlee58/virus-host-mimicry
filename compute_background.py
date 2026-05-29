#!/usr/bin/env python3
"""Count proteins in the Foldseek target DB per species of interest.

Produces the background (N, K_species) for the hypergeometric taxonomic
enrichment test (Lasso et al. 2021, Figure 1B).

Reads the Foldseek DB's `_mapping` file (one line per DB entry: internal_id\\ttaxid)
and walks the NCBI taxonomy tree from each entry's taxid up to root, marking
the entry as belonging to a species if that species' taxid appears in the
ancestor chain.

This is a one-time computation. The output TSV is then consumed by
enrichment_fig1b.py.

Example:
    python3 compute_background.py
        -> writes background_counts.tsv next to the script

    python3 compute_background.py --db /path/to/foldseek_db -o bg.tsv
"""
import argparse
from collections import defaultdict
from pathlib import Path
import sys

# Species taxids we test in Figure 1B.
# Counting at species level means a strain-level taxid (e.g. E. coli K-12 = 83333)
# still counts toward "E. coli" because 562 appears in its ancestor chain.
SPECIES = {
    "ecoli": 562,    # Escherichia coli
    "yeast": 4932,   # Saccharomyces cerevisiae
    "human": 9606,   # Homo sapiens
}


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


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--db", default="/fast/databases/foldseek/pdb/pdb2025/pdb",
        help="Foldseek DB path prefix (default: %(default)s)",
    )
    ap.add_argument(
        "--taxdump", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/taxdump"),
        help="Directory containing nodes.dmp and merged.dmp",
    )
    ap.add_argument(
        "-o", "--output", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/taxonomic_enrichment/background_counts.tsv"),
    )
    args = ap.parse_args()

    mapping_path = Path(args.db + "_mapping")
    nodes_path = args.taxdump / "nodes.dmp"
    merged_path = args.taxdump / "merged.dmp"

    for p in (mapping_path, nodes_path, merged_path):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")

    print(f"Loading {nodes_path}...", file=sys.stderr)
    parents = load_parents(nodes_path)
    print(f"  {len(parents):,} taxid entries", file=sys.stderr)

    print(f"Loading {merged_path}...", file=sys.stderr)
    merged = load_merged(merged_path)
    print(f"  {len(merged):,} merged taxid redirects", file=sys.stderr)

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
