#!/bin/bash
#SBATCH --job-name=fs_array
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1                 # 1 GPU per chunk -> many chunks run concurrently
#SBATCH --nodes=1
#SBATCH --nodelist=devlss001,devlss002,devbox002  # healthy GPU nodes; nodes=1 -> each task lands on one (skips broken devbox001)
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00              # one ~50k chunk on 1 GPU is ~2h; headroom for slower A5000 nodes
#SBATCH --output=logs/fs_array_%A_%a.log

# SLURM-array Foldseek: each array task searches ONE chunk produced by split_fasta.py.
# Each task asks for gpu:1, so with the three nodelist nodes (6+8+8 = 22 GPUs)
# SLURM packs up to 22 chunks across them at once. devbox001 is excluded (broken driver).
#
# Usage:
#   1. split:   python3 split_fasta.py bacteria.faa --outdir bacteria_chunks --per-chunk 50000
#               (prints K = number of chunks)
#   2. submit:  sbatch --array=0-<K-1> run_foldseek_array.sh bacteria_chunks bacteria_foldseek
#               (the nodelist already bounds concurrency to 22 GPUs; add %N to cap lower if sharing)
#   3. combine: cat bacteria_foldseek/hits_*.m8 > bacteria_foldseek/hits.m8
#
# Optional env vars: TARGET_DB, PROSTT5_WEIGHTS (same defaults as run_foldseek.sh).

set -euo pipefail

CHUNK_DIR="${1:?Usage: sbatch --array=0-N run_foldseek_array.sh <chunk_dir> <output_dir>}"
OUTPUT="${2:?Usage: sbatch --array=0-N run_foldseek_array.sh <chunk_dir> <output_dir>}"
TARGET_DB="${TARGET_DB:-/fast2/yewon1/AFCDB_analysis_data/foldseek_search_PDBe/foldseek_pdb_db/gpu_pdb}"
PROSTT5_WEIGHTS="${PROSTT5_WEIGHTS:-/fast/sunny/virus-host-mimicry/prostt5}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: not an array job. Submit with: sbatch --array=0-<K-1> $0 <chunk_dir> <output_dir>" >&2
  exit 1
fi

CHUNK=$(printf "%s/chunk_%04d.faa" "$CHUNK_DIR" "$SLURM_ARRAY_TASK_ID")
if [[ ! -f "$CHUNK" ]]; then
  echo "ERROR: chunk not found: $CHUNK" >&2
  exit 1
fi

mkdir -p "$OUTPUT"
# Per-task temp dir so concurrent array tasks never collide.
TMP="${SCRATCH:-/tmp}/fs_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

OUT_M8="$OUTPUT/hits_$(printf '%04d' "$SLURM_ARRAY_TASK_ID").m8"

start=$(date +%s)
srun foldseek easy-search \
  "$CHUNK" \
  "$TARGET_DB" \
  "$OUT_M8" \
  "$TMP" \
  --prostt5-model "$PROSTT5_WEIGHTS" \
  --alignment-type 0 \
  --gpu 1 \
  --prefilter-mode 1 \
  --format-output "query,target,evalue,bits,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK}"
echo "chunk $SLURM_ARRAY_TASK_ID done in $(( $(date +%s) - start )) s -> $OUT_M8"
