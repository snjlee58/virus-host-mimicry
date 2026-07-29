#!/bin/bash
# Prepare BFVD2 query chunks for the SLURM-array Foldseek run vs FULL AFDB.
#
# Same idea as the afdb_proteome version: list every PDB in the BFVD2 folder and
# split that list into fixed-size chunk manifests (chunks/list_0000, ...). Chunks
# are LARGE here (250k default -> ~5 chunks) because each chunk re-scans the whole
# 241M-entry target, so fewer/larger chunks avoid redundant full-target scans.
#
# Run ONCE on the login shell (just listing + splitting text). Prints K + the
# exact sbatch command.
#
# Usage: bash prepare_chunks.sh [bfvd2_dir] [output_dir] [chunk_size]

set -euo pipefail

BFVD2="${1:-/fast/sunny/bfvd2/bfvd2_pdb_only}"
OUT="${2:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb}"
CHUNK_SIZE="${3:-250000}"  # FEW large chunks: each chunk re-scans the whole 241M target, so minimize chunk count

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
echo "PILOT ONE CHUNK FIRST (full AFDB is a different beast — measure before committing):"
echo "  sbatch --array=0 run_foldseek_array.sh"
echo "then, once its runtime/memory look OK:"
echo "  sbatch --array=1-$((K - 1))%3 run_foldseek_array.sh"
echo ""
echo "After all tasks finish:  bash merge_hits.sh"
