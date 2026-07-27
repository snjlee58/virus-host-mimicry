#!/usr/bin/env python3
"""Map each BFVD entry to its virus's host taxonomic division(s).

This is the bridge between the BFVD query set (Fig 1C reproduction) and the
Lasso et al. host-division grouping: only viral proteins whose virus has a
known host in Virus-Host DB can be placed in a host-infecting group.

Chain, per BFVD entry:
    entry --(col2 of bfvd_taxid_rank_scientificname_lineage.tsv)--> virus taxid
          --(virushostdb.tsv)-->                                    host taxid(s)
          --(nodes.dmp division, classify_hits.classify)-->         host bucket

Buckets match classify_hits.py (bacteria, invertebrate, vertebrate, human,
plantfungi, phage, virus, unknown) plus two markers specific to this step:
  - "no_host"     : virus absent from Virus-Host DB (no record at all)
  - "unspecified" : host recorded as taxid 1 ("root"), VHDB's placeholder for
                    an unspecified host -> not usable for grouping

Only entries with >=1 host in a REAL division (bacteria/invertebrate/vertebrate/
human/plantfungi) are "groupable" and enter the Fig 1C enrichment; those are the
ones written to --query-list. Multi-host viruses produce one row per host and can
land in more than one group (as in the paper's non-exclusive grouping).

Output (TSV, long format):
    bfvd_entry  accession  virus_taxid  virus_name  host_taxid  host_name  host_bucket

Example (local paths):
    python3 map_bfvd_to_host.py \
        --bfvd-taxids bfvd_2023_02_v2/bfvd_taxid_rank_scientificname_lineage.tsv \
        --virushostdb virushostdb/virushostdb.tsv \
        --taxdump taxdump \
        -o bfvd_2023_02_v2/bfvd_host_mapping.tsv \
        --query-list bfvd_2023_02_v2/bfvd_host_mapped_entries.txt
"""
import argparse
import sys
from pathlib import Path

from classify_hits import load_nodes, load_merged, classify
from classify_viruses import load_virushostdb

# Real host divisions usable in the Fig 1C grouping (human is also nested under
# vertebrate at the enrichment step; kept separate here, per classify_hits).
REAL_DIVISIONS = {"bacteria", "invertebrate", "vertebrate", "human", "plantfungi"}
ROOT_TAXID = 1  # VHDB placeholder for an unspecified host


def parse_bfvd_taxids(path):
    """Yield (entry, accession, virus_taxid, virus_name). File has no header.

    Columns: entry_filename, taxid, rank, scientific_name, lineage.
    UniProt accession is the entry-name prefix before the first '_'.
    """
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            entry = fields[0]
            try:
                vtax = int(fields[1])
            except ValueError:
                continue  # missing/non-numeric taxid
            acc = entry.split("_", 1)[0]
            vname = fields[3] if len(fields) > 3 else ""
            yield entry, acc, vtax, vname


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--bfvd-taxids", type=Path,
        default=Path("/fast/sunny/bfvd/2023_02_v2/bfvd_taxid_rank_scientificname_lineage.tsv"),
        help="BFVD taxid TSV (entry, taxid, rank, name, lineage; no header)")
    ap.add_argument(
        "--virushostdb", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/virushostdb/virushostdb.tsv"))
    ap.add_argument(
        "--taxdump", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/taxdump"),
        help="Directory containing nodes.dmp and merged.dmp")
    ap.add_argument("-o", "--output", type=Path, default=None,
                    help="Output TSV (default: stdout)")
    ap.add_argument("--host-mapped-only", action="store_true",
                    help="Emit only entries with a VHDB record (drop no_host rows)")
    ap.add_argument("--query-list", type=Path, default=None,
                    help="Also write the unique GROUPABLE BFVD entry names here "
                         "(>=1 real host division; one per line) for optional "
                         "foldseek query pre-filtering")
    args = ap.parse_args()

    for p in (args.bfvd_taxids, args.virushostdb,
              args.taxdump / "nodes.dmp", args.taxdump / "merged.dmp"):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")

    print(f"Loading {args.taxdump}/nodes.dmp ...", file=sys.stderr)
    nodes = load_nodes(args.taxdump / "nodes.dmp")
    print(f"Loading {args.taxdump}/merged.dmp ...", file=sys.stderr)
    merged = load_merged(args.taxdump / "merged.dmp")
    print(f"Loading {args.virushostdb} ...", file=sys.stderr)
    vhdb = load_virushostdb(args.virushostdb)
    print(f"  {len(vhdb):,} viruses with >=1 host in Virus-Host DB", file=sys.stderr)

    out = open(args.output, "w") if args.output else sys.stdout
    out.write("bfvd_entry\taccession\tvirus_taxid\tvirus_name\t"
              "host_taxid\thost_name\thost_bucket\n")
    qlist = open(args.query_list, "w") if args.query_list else None

    n_entries = n_record = n_groupable = 0
    entries_by_bucket = {}   # bucket -> set(entry)
    viruses_by_bucket = {}   # bucket -> set(virus_taxid)
    no_host_viruses = set()

    for entry, acc, vtax, vname in parse_bfvd_taxids(args.bfvd_taxids):
        n_entries += 1
        # Exact virus-taxid join, with a merged.dmp rescue for stale BFVD taxids.
        hosts = vhdb.get(vtax)
        if hosts is None:
            hosts = vhdb.get(merged.get(vtax, vtax))
        if not hosts:
            no_host_viruses.add(vtax)
            if not args.host_mapped_only:
                out.write(f"{entry}\t{acc}\t{vtax}\t{vname}\t-\t-\tno_host\n")
            continue
        n_record += 1
        # One BFVD entry -> one virus, so all its host rows are produced here
        # together; decide groupability per entry in this same pass.
        ent_buckets = set()
        for (htid, _vn, hname) in hosts:
            bucket = "unspecified" if htid == ROOT_TAXID else classify(str(htid), nodes, merged)
            out.write(f"{entry}\t{acc}\t{vtax}\t{vname}\t{htid}\t{hname}\t{bucket}\n")
            ent_buckets.add(bucket)
            entries_by_bucket.setdefault(bucket, set()).add(entry)
            viruses_by_bucket.setdefault(bucket, set()).add(vtax)
        if ent_buckets & REAL_DIVISIONS:
            n_groupable += 1
            if qlist:
                qlist.write(entry + "\n")

    if args.output:
        out.close()
    if qlist:
        qlist.close()

    print("\n=== Coverage (BFVD entries) ===", file=sys.stderr)
    print(f"total:                              {n_entries:,}", file=sys.stderr)
    print(f"has a Virus-Host DB record:         {n_record:,} "
          f"({100.0*n_record/n_entries:.1f}%)", file=sys.stderr)
    print(f"  GROUPABLE (>=1 real host div):    {n_groupable:,} "
          f"({100.0*n_groupable/n_entries:.1f}%)   <- Fig 1C query set", file=sys.stderr)
    print(f"  record but ungroupable (root/NA): {n_record - n_groupable:,}", file=sys.stderr)
    print(f"no Virus-Host DB record:            {n_entries - n_record:,} "
          f"({100.0*(n_entries-n_record)/n_entries:.1f}%)", file=sys.stderr)
    print(f"(virus taxids with no record:       {len(no_host_viruses):,})", file=sys.stderr)

    print("\n=== Host-division groups (multi-host viruses counted in each) ===",
          file=sys.stderr)
    print(f"{'bucket':<14}{'viruses':>10}{'entries':>12}", file=sys.stderr)
    for b in sorted(entries_by_bucket, key=lambda x: -len(entries_by_bucket[x])):
        print(f"{b:<14}{len(viruses_by_bucket[b]):>10,}{len(entries_by_bucket[b]):>12,}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
