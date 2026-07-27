#!/bin/bash
# Prepare BFVD2 query chunks for the SLURM-array Foldseek run.
#
# Lists every PDB in the BFVD2 folder and splits that list into fixed-size chunk
# manifests: chunks/list_0000, list_0001, ... (each is a list of PDB paths, one
# per line). Each array task later symlinks its chunk's PDBs and searches them.
#
# Run this ONCE, on the login shell (it's just listing + splitting text — no heavy
# compute). It prints K (the number of chunks) and the exact sbatch command.
#
# Usage: bash prepare_chunks.sh [bfvd2_dir] [output_dir] [chunk_size]

set -euo pipefail

BFVD2="${1:-/fast/sunny/bfvd2/bfvd2_pdb_only}"
OUT="${2:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_proteome}"
CHUNK_SIZE="${3:-20000}"   # PDBs per chunk; aim for K ~ 2-3x your concurrency

[[ -d "$BFVD2" ]] || { echo "ERROR: BFVD2 dir not found: $BFVD2" >&2; exit 1; }

mkdir -p "$OUT/chunks" "$OUT/hits" "$OUT/logs"

echo "Listing PDBs in $BFVD2 ..."
find "$BFVD2" -maxdepth 1 -name '*.pdb' | sort > "$OUT/all_pdbs.txt"
N=$(wc -l < "$OUT/all_pdbs.txt")
[[ "$N" -gt 0 ]] || { echo "ERROR: no *.pdb files found in $BFVD2 (gzipped? adjust the -name glob)" >&2; exit 1; }

rm -f "$OUT"/chunks/list_*
split -d -a 4 -l "$CHUNK_SIZE" "$OUT/all_pdbs.txt" "$OUT/chunks/list_"
K=$(find "$OUT/chunks" -name 'list_*' | wc -l | tr -d ' ')

echo ""
echo "  PDBs (N):     $N"
echo "  chunk size:   $CHUNK_SIZE"
echo "  chunks (K):   $K"
echo ""
echo "Now submit the array (from this folder):"
echo "  sbatch --array=0-$((K - 1))%10 run_foldseek_array.sh"
echo ""
echo "(%10 = at most 10 chunks at once. Raise it if the cluster is quiet, lower if busy.)"
echo "After all tasks finish:  bash merge_hits.sh"
