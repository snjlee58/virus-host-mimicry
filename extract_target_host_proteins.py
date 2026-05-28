#!/usr/bin/env python3
"""Extract all viral proteins encoded by viruses infecting a given host.

Workflow:
  1. Read virushostdb.tsv -> find viruses whose host_tax_id matches the input.
  2. Collect all source genome accessions for those viruses (column 4 in
     virushostdb.tsv, comma-separated for multipartite viruses).
  3. Stream virushostdb.formatted.cds.faa.gz -> emit FASTA records whose
     source genome accession (header field 6) is in the set.

Output: a FASTA file with the matching viral proteins, ready for foldseek.

Example:
    python3 extract_target_host_proteins.py 562 \\
        -o /fast/sunny/virus-host-mimicry/queries/ecoli_phages.faa
"""
import argparse
import gzip
import sys
from pathlib import Path

# virushostdb.tsv column indices (0-based)
VHDB_VIRUS_TAXID = 0
VHDB_REFSEQ_IDS = 3       # comma-separated for multipartite viruses
VHDB_HOST_TAXID = 7

# FASTA header pipe-field index for the source genome accession (0-based)
FASTA_GENOME_FIELD = 5


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


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("host_taxid", type=int,
                    help="NCBI taxid of the host (e.g., 562 for E. coli)")
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
    ap.add_argument("-o", "--output", type=Path, required=True,
                    help="Output FASTA path")
    args = ap.parse_args()

    for p in (args.virushostdb, args.fasta):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")

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
