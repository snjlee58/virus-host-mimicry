#!/usr/bin/env python3
"""Split a (multi-)FASTA into fixed-size chunks for SLURM-array Foldseek.

Writes <outdir>/chunk_0000.faa, chunk_0001.faa, ... each holding up to
--per-chunk sequence records. Prints the number of chunks created (K), so you
can submit the array as: sbatch --array=0-(K-1) run_foldseek_array.sh ...

Example:
    python3 split_fasta.py tests/host_divisions/bacteria.faa \\
        --outdir tests/host_divisions/bacteria_chunks --per-chunk 100000
"""
import argparse
import sys
from pathlib import Path


def split(fasta_path, outdir, per_chunk, prefix):
    outdir.mkdir(parents=True, exist_ok=True)
    idx = -1
    n_in_chunk = per_chunk  # force a new chunk on the first record
    n_records = 0
    out = None
    try:
        with open(fasta_path) as fin:
            for line in fin:
                if line.startswith(">"):
                    if n_in_chunk >= per_chunk:
                        if out:
                            out.close()
                        idx += 1
                        out = open(outdir / f"{prefix}{idx:04d}.faa", "w")
                        n_in_chunk = 0
                    n_in_chunk += 1
                    n_records += 1
                if out is not None:
                    out.write(line)
    finally:
        if out:
            out.close()
    return idx + 1, n_records


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", type=Path, help="Input (multi-)FASTA")
    ap.add_argument("--outdir", type=Path, required=True, help="Chunk output directory")
    ap.add_argument("--per-chunk", type=int, default=100_000,
                    help="Max records per chunk (default: %(default)s)")
    ap.add_argument("--prefix", default="chunk_", help="Chunk filename prefix")
    args = ap.parse_args()

    if not args.fasta.exists():
        sys.exit(f"ERROR: file not found: {args.fasta}")

    n_chunks, n_records = split(args.fasta, args.outdir, args.per_chunk, args.prefix)
    print(f"{n_records:,} records -> {n_chunks} chunks in {args.outdir}", file=sys.stderr)
    if n_chunks:
        print(f"Submit with: --array=0-{n_chunks - 1}", file=sys.stderr)
    # Also print the bare count to stdout for scripting.
    print(n_chunks)


if __name__ == "__main__":
    main()
