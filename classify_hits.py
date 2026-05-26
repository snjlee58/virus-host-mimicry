#!/usr/bin/env python3
"""Classify foldseek hits.m8 rows by host taxonomic bucket.

Buckets follow the Lasso et al. (Cell Systems 2021) Figure 1C scheme:
bacteria, plantfungi, invertebrate, vertebrate, human, virus, phage, unknown.

Classification uses NCBI taxonomic divisions from nodes.dmp (the paper's
method), with merged.dmp applied to rescue stale strain taxids.

Example:
    python3 classify_hits.py virprot1_output/foldseek/hits.m8
        -> writes virprot1_output/foldseek/hits.classified.tsv

    python3 classify_hits.py hits.m8 --taxdump taxdump/ -o out.tsv
"""
import argparse
import sys
from pathlib import Path

# NCBI division_id -> bucket (per division.dmp)
DIV_TO_BUCKET = {
    0:  "bacteria",      # BCT
    1:  "invertebrate",  # INV
    2:  "vertebrate",    # MAM (non-primate, non-rodent)
    3:  "phage",         # PHG
    4:  "plantfungi",    # PLN
    5:  "vertebrate",    # PRI (human handled separately by taxid)
    6:  "vertebrate",    # ROD
    7:  "unknown",       # SYN
    8:  "unknown",       # UNA
    9:  "virus",         # VRL
    10: "vertebrate",    # VRT (non-mammal vertebrates)
    11: "unknown",       # ENV
}

HUMAN_TAXID = 9606
TAXID_COL = 4


def _split_dmp(line):
    # NCBI .dmp lines end with "\t|\n"; strip that then split inner "\t|\t".
    return line.rstrip("\n").rstrip("|").rstrip("\t").split("\t|\t")


def load_nodes(path):
    """Parse nodes.dmp -> {taxid: division_id}."""
    d = {}
    with open(path) as f:
        for line in f:
            fields = _split_dmp(line)
            d[int(fields[0])] = int(fields[4])
    return d


def load_merged(path):
    """Parse merged.dmp -> {old_taxid: new_taxid}."""
    d = {}
    with open(path) as f:
        for line in f:
            fields = _split_dmp(line)
            d[int(fields[0])] = int(fields[1])
    return d


def classify(taxid_str, nodes, merged):
    if not taxid_str or taxid_str in ("0", "-"):
        return "unknown"
    try:
        taxid = int(taxid_str)
    except ValueError:
        return "unknown"
    taxid = merged.get(taxid, taxid)
    if taxid == HUMAN_TAXID:
        return "human"
    div = nodes.get(taxid)
    if div is None:
        return "unknown"
    return DIV_TO_BUCKET.get(div, "unknown")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("hits", type=Path, help="Path to foldseek hits.m8")
    ap.add_argument(
        "--taxdump", type=Path, default=Path("taxdump"),
        help="Directory containing nodes.dmp and merged.dmp (default: taxdump/)",
    )
    ap.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output TSV path (default: <hits_stem>.classified.tsv next to input)",
    )
    args = ap.parse_args()

    if args.output is None:
        args.output = args.hits.parent / (args.hits.stem + ".classified.tsv")

    nodes_path = args.taxdump / "nodes.dmp"
    merged_path = args.taxdump / "merged.dmp"
    for p in (args.hits, nodes_path, merged_path):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")

    print(f"Loading {nodes_path}...", file=sys.stderr)
    nodes = load_nodes(nodes_path)
    print(f"  {len(nodes):,} taxid entries", file=sys.stderr)

    print(f"Loading {merged_path}...", file=sys.stderr)
    merged = load_merged(merged_path)
    print(f"  {len(merged):,} merged taxid redirects", file=sys.stderr)

    print(f"Classifying {args.hits} -> {args.output}", file=sys.stderr)
    counts = {}
    n_rows = 0
    with open(args.hits) as fin, open(args.output, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            taxid_str = fields[TAXID_COL] if len(fields) > TAXID_COL else ""
            bucket = classify(taxid_str, nodes, merged)
            counts[bucket] = counts.get(bucket, 0) + 1
            n_rows += 1
            fout.write(line + "\t" + bucket + "\n")

    print(f"Done: {n_rows:,} rows", file=sys.stderr)
    print("Bucket distribution:", file=sys.stderr)
    for bucket, count in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {bucket}: {count:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
