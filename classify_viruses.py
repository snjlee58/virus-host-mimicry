#!/usr/bin/env python3
"""Classify viruses by the taxonomic division of their host(s).

For each virus taxid, looks up host(s) in Virus-Host DB (virushostdb.tsv),
then classifies each host using NCBI division codes (same logic as
classify_hits.py). Output is long-format TSV, one row per virus-host pair:

    virus_taxid  virus_name  host_taxid  host_name  host_bucket

Multi-host viruses produce multiple rows.

Example:
    python3 classify_viruses.py 12345 67890
    python3 classify_viruses.py --from-file query_viruses.txt -o virus_buckets.tsv
"""
import argparse
import sys
from pathlib import Path

from classify_hits import load_nodes, load_merged, classify

# virushostdb.tsv column indices (0-based), per README:
# virus tax id | virus name | virus lineage | refseq id | KEGG GENOME |
# KEGG DISEASE | DISEASE | host tax id | host name | host lineage | ...
VHDB_VIRUS_TAXID = 0
VHDB_VIRUS_NAME = 1
VHDB_HOST_TAXID = 7
VHDB_HOST_NAME = 8


def load_virushostdb(path):
    """Parse virushostdb.tsv -> {virus_taxid: [(host_taxid, virus_name, host_name), ...]}."""
    d = {}
    with open(path) as f:
        header = f.readline()  # skip header
        if not header.lower().startswith("# virus tax id") and "virus tax id" not in header.lower():
            print(f"WARN: virushostdb.tsv header looks unexpected: {header.rstrip()[:80]}",
                  file=sys.stderr)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= VHDB_HOST_NAME:
                continue
            try:
                virus_taxid = int(fields[VHDB_VIRUS_TAXID])
                host_taxid = int(fields[VHDB_HOST_TAXID])
            except ValueError:
                continue  # missing/non-numeric taxid (e.g. metagenome rows)
            virus_name = fields[VHDB_VIRUS_NAME]
            host_name = fields[VHDB_HOST_NAME]
            d.setdefault(virus_taxid, []).append((host_taxid, virus_name, host_name))
    return d


def classify_virus(virus_taxid, vhdb, nodes, merged):
    """Return list of (virus_name, host_taxid, host_name, host_bucket).

    Empty list if virus_taxid not in virushostdb.
    """
    hosts = vhdb.get(virus_taxid, [])
    return [
        (vname, htid, hname, classify(str(htid), nodes, merged))
        for (htid, vname, hname) in hosts
    ]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("virus_taxids", nargs="*", type=int,
                    help="One or more virus taxids to classify")
    ap.add_argument("--from-file", type=Path, default=None,
                    help="File with one virus taxid per line")
    ap.add_argument(
        "--virushostdb", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/virushostdb/virushostdb.tsv"),
        help="Path to virushostdb.tsv",
    )
    ap.add_argument(
        "--taxdump", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/taxdump"),
        help="Directory containing nodes.dmp and merged.dmp",
    )
    ap.add_argument("-o", "--output", type=Path, default=None,
                    help="Output TSV (default: stdout)")
    args = ap.parse_args()

    # Collect input taxids
    taxids = list(args.virus_taxids)
    if args.from_file:
        with open(args.from_file) as f:
            taxids.extend(int(line.strip()) for line in f if line.strip())
    if not taxids:
        sys.exit("ERROR: provide virus taxids as args or via --from-file")

    # Sanity-check inputs exist
    for p in (args.virushostdb, args.taxdump / "nodes.dmp", args.taxdump / "merged.dmp"):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")

    print(f"Loading {args.taxdump}/nodes.dmp...", file=sys.stderr)
    nodes = load_nodes(args.taxdump / "nodes.dmp")
    print(f"  {len(nodes):,} taxid entries", file=sys.stderr)

    print(f"Loading {args.taxdump}/merged.dmp...", file=sys.stderr)
    merged = load_merged(args.taxdump / "merged.dmp")
    print(f"  {len(merged):,} merged taxid redirects", file=sys.stderr)

    print(f"Loading {args.virushostdb}...", file=sys.stderr)
    vhdb = load_virushostdb(args.virushostdb)
    print(f"  {len(vhdb):,} viruses with at least one host", file=sys.stderr)

    out = open(args.output, "w") if args.output else sys.stdout
    out.write("virus_taxid\tvirus_name\thost_taxid\thost_name\thost_bucket\n")

    n_found = n_missing = 0
    for vid in taxids:
        rows = classify_virus(vid, vhdb, nodes, merged)
        if not rows:
            n_missing += 1
            out.write(f"{vid}\t-\t-\t-\tunknown\n")
            continue
        n_found += 1
        for vname, htid, hname, bucket in rows:
            out.write(f"{vid}\t{vname}\t{htid}\t{hname}\t{bucket}\n")

    if args.output:
        out.close()

    print(f"Done: {n_found:,} viruses classified, {n_missing:,} not found in virushostdb",
          file=sys.stderr)


if __name__ == "__main__":
    main()
