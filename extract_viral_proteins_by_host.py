#!/usr/bin/env python3
"""Extract all viral proteins encoded by viruses infecting a given host.

Workflow:
  1. Read virushostdb.tsv -> find viruses whose host_tax_id matches the input.
  2. Collect all source genome accessions for those viruses (column 4 in
     virushostdb.tsv, comma-separated for multipartite viruses).
  3. Stream virushostdb.formatted.cds.faa.gz -> emit FASTA records whose
     source genome accession (header field 6) is in the set.

Output: a FASTA file with the matching viral proteins, ready for foldseek.

With --all-divisions, instead classifies every host in virushostdb.tsv into one
of the five Fig 1C taxonomic divisions and emits one FASTA per division
(bacteria/plantfungi/invertebrate/vertebrate/human) in a single pass over the
big FASTA. A virus with hosts in several divisions contributes to each.

Example:
    # single host (Fig 1B style)
    python3 extract_viral_proteins_by_host.py 562 \\
        -o /fast/sunny/virus-host-mimicry/queries/ecoli_phages.faa

    # all five divisions (Fig 1C)
    python3 extract_viral_proteins_by_host.py --all-divisions \\
        -o /fast/sunny/virus-host-mimicry/tests/host_divisions
"""
import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

from classify_hits import load_nodes, load_merged, classify

# virushostdb.tsv column indices (0-based)
VHDB_VIRUS_TAXID = 0
VHDB_REFSEQ_IDS = 3       # comma-separated for multipartite viruses
VHDB_HOST_TAXID = 7

# FASTA header pipe-field index for the source genome accession (0-based)
FASTA_GENOME_FIELD = 5

# The five host taxonomic divisions reproduced as columns/rows of Fig 1C.
# (classify() may also return virus/phage/unknown; those hosts are skipped.)
PAPER_DIVISIONS = ["bacteria", "plantfungi", "invertebrate", "vertebrate", "human"]


def collect_genome_accessions(vhdb_path, host_taxid):
    """Return (virus_taxids: set[int], genome_accs: set[str]) for the given host."""
    virus_taxids = set()
    genome_accs = set()
    with open(vhdb_path) as f:
        header = f.readline()
        if "host tax id" not in header.lower():
            print(f"WARN: virushostdb.tsv header looks unexpected: {header.rstrip()[:80]}",
                  file=sys.stderr)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= VHDB_HOST_TAXID:
                continue
            try:
                if int(fields[VHDB_HOST_TAXID]) != host_taxid:
                    continue
            except ValueError:
                continue
            try:
                virus_taxids.add(int(fields[VHDB_VIRUS_TAXID]))
            except ValueError:
                pass
            for acc in fields[VHDB_REFSEQ_IDS].split(", "):
                acc = acc.strip()
                if acc:
                    genome_accs.add(acc)
    return virus_taxids, genome_accs


def extract_proteins(fasta_path, genome_accs, out_path):
    """Stream FASTA, emit records whose source genome accession is in genome_accs.

    Returns the number of records written.
    """
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    n_records = 0
    keep = False
    with opener(fasta_path, "rt") as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fields = line[1:].split("|")
                genome = (fields[FASTA_GENOME_FIELD].strip()
                          if len(fields) > FASTA_GENOME_FIELD else "")
                keep = genome in genome_accs
                if keep:
                    n_records += 1
            if keep:
                fout.write(line)
    return n_records


def classify_hosts_to_accessions(vhdb_path, nodes, merged, divisions):
    """Single pass over virushostdb.tsv -> {division: set(genome_accessions)}.

    Each virus-host row is classified by its host's NCBI division (via the same
    classify() used elsewhere). The virus's genome accessions are added to the
    bucket for that host's division. Hosts whose bucket is not one of `divisions`
    (e.g. virus/phage/unknown) are skipped. A virus with hosts in several
    divisions (e.g. an arbovirus: human + mosquito) contributes its accessions to
    each of those divisions.
    """
    accs = {d: set() for d in divisions}
    skipped = defaultdict(int)
    with open(vhdb_path) as f:
        header = f.readline()
        if "host tax id" not in header.lower():
            print(f"WARN: virushostdb.tsv header looks unexpected: {header.rstrip()[:80]}",
                  file=sys.stderr)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= VHDB_HOST_TAXID:
                continue
            bucket = classify(fields[VHDB_HOST_TAXID], nodes, merged)
            if bucket not in accs:
                skipped[bucket] += 1
                continue
            for acc in fields[VHDB_REFSEQ_IDS].split(", "):
                acc = acc.strip()
                if acc:
                    accs[bucket].add(acc)
    return accs, skipped


def extract_proteins_multi(fasta_path, accs_by_division, out_dir):
    """Stream the FASTA once, routing each record to every division it belongs to.

    Writes <out_dir>/<division>.faa for each division. Returns {division: n_records}.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    handles = {d: open(out_dir / f"{d}.faa", "w") for d in accs_by_division}
    counts = {d: 0 for d in accs_by_division}
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    targets = []  # divisions the current record is being written to
    try:
        with opener(fasta_path, "rt") as fin:
            for line in fin:
                if line.startswith(">"):
                    fields = line[1:].split("|")
                    genome = (fields[FASTA_GENOME_FIELD].strip()
                              if len(fields) > FASTA_GENOME_FIELD else "")
                    targets = [d for d, accset in accs_by_division.items()
                               if genome in accset]
                    for d in targets:
                        counts[d] += 1
                for d in targets:
                    handles[d].write(line)
    finally:
        for h in handles.values():
            h.close()
    return counts


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("host_taxid", type=int, nargs="?",
                    help="NCBI taxid of the host (e.g., 562 for E. coli). "
                         "Omit when using --all-divisions.")
    ap.add_argument(
        "--all-divisions", action="store_true",
        help="Emit one FASTA per host taxonomic division "
             f"({', '.join(PAPER_DIVISIONS)}) in a single pass over the FASTA. "
             "-o is then treated as an output directory.",
    )
    ap.add_argument(
        "--virushostdb", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/virushostdb/virushostdb.tsv"),
        help="Path to virushostdb.tsv",
    )
    ap.add_argument(
        "--fasta", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/virushostdb/virushostdb.formatted.cds.faa.gz"),
        help="Path to virushostdb.formatted.cds.faa.gz",
    )
    ap.add_argument(
        "--taxdump", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/taxdump"),
        help="Directory with nodes.dmp and merged.dmp (only used by --all-divisions)",
    )
    ap.add_argument("-o", "--output", type=Path, required=True,
                    help="Output FASTA path (single-host mode) or "
                         "output directory (--all-divisions)")
    args = ap.parse_args()

    if args.all_divisions == (args.host_taxid is not None):
        sys.exit("ERROR: provide exactly one of: a host_taxid (single-host mode) "
                 "OR --all-divisions")

    for p in (args.virushostdb, args.fasta):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")

    if args.all_divisions:
        nodes_path = args.taxdump / "nodes.dmp"
        merged_path = args.taxdump / "merged.dmp"
        for p in (nodes_path, merged_path):
            if not p.exists():
                sys.exit(f"ERROR: file not found: {p}")

        print(f"Loading {nodes_path}...", file=sys.stderr)
        nodes = load_nodes(nodes_path)
        print(f"  {len(nodes):,} taxid entries", file=sys.stderr)
        print(f"Loading {merged_path}...", file=sys.stderr)
        merged = load_merged(merged_path)
        print(f"  {len(merged):,} merged taxid redirects", file=sys.stderr)

        print(f"Classifying hosts in {args.virushostdb}...", file=sys.stderr)
        accs, skipped = classify_hosts_to_accessions(
            args.virushostdb, nodes, merged, PAPER_DIVISIONS)
        for d in PAPER_DIVISIONS:
            print(f"  {d:>12}: {len(accs[d]):,} genome accessions", file=sys.stderr)
        if skipped:
            skip_str = ", ".join(f"{b}={n:,}" for b, n in sorted(skipped.items()))
            print(f"  (skipped host rows: {skip_str})", file=sys.stderr)

        print(f"Extracting proteins from {args.fasta} -> {args.output}/",
              file=sys.stderr)
        counts = extract_proteins_multi(args.fasta, accs, args.output)
        print("Done:", file=sys.stderr)
        for d in PAPER_DIVISIONS:
            print(f"  {d}.faa: {counts[d]:,} proteins", file=sys.stderr)
        return

    print(f"Scanning {args.virushostdb} for host_taxid={args.host_taxid}...",
          file=sys.stderr)
    virus_taxids, genome_accs = collect_genome_accessions(args.virushostdb, args.host_taxid)
    print(f"  found {len(virus_taxids):,} viruses with {len(genome_accs):,} genome accessions",
          file=sys.stderr)
    if not genome_accs:
        sys.exit(f"ERROR: no viruses found with host_taxid={args.host_taxid}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Extracting proteins from {args.fasta} -> {args.output}", file=sys.stderr)
    n = extract_proteins(args.fasta, genome_accs, args.output)
    print(f"Done: {n:,} protein sequences written", file=sys.stderr)


if __name__ == "__main__":
    main()
