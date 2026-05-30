#!/bin/bash
#SBATCH --job-name=foldseek
#SBATCH --partition=gpu                 # GPU partition (needed for ProstT5 inference)
#SBATCH --gres=gpu:6                    # All 6 L40S on devlss001; foldseek spreads the target DB across them
#SBATCH --nodelist=devlss001            # devbox001 has a broken NVML driver, pin to devlss001
#SBATCH --cpus-per-task=64              # All 64 cores on devlss001; helps the post-prefilter alignment stage
#SBATCH --mem=128G                      # Plenty of headroom; bump to 256G for AFDB-scale searches
#SBATCH --time=04:00:00
#SBATCH --output=logs/foldseek_%j.log

# Usage:
#   sbatch run_foldseek.sh <query> <output_dir>
#
# <query> may be either:
#   - a PDB/CIF structure (.pdb, .cif, .ent)         → uses 3Di from structure
#   - an amino-acid FASTA (.fa, .faa, .fasta)        → uses ProstT5 to predict 3Di from sequence
#
# Optional env vars:
#   TARGET_DB         path to Foldseek database (default: /fast/databases/foldseek/pdb/pdb100)
#   PROSTT5_WEIGHTS   path to ProstT5 weights dir   (default: /fast/sunny/virus-host-mimicry/prostt5)
#                     only used when <query> is a FASTA
#
# Examples:
#   sbatch run_foldseek.sh virprot1_output/top.pdb virprot1_output/foldseek
#   TARGET_DB=/fast/databases/foldseek/afdb sbatch run_foldseek.sh queries.fasta fs_afdb

set -euo pipefail

# --- Args ---
QUERY="${1:?Usage: sbatch run_foldseek.sh <query> <output_dir>}"
OUTPUT="${2:?Usage: sbatch run_foldseek.sh <query> <output_dir>}"
TARGET_DB="${TARGET_DB:-/fast2/yewon1/AFCDB_analysis_data/foldseek_search_PDBe/foldseek_pdb_db/gpu_pdb}"
PROSTT5_WEIGHTS="${PROSTT5_WEIGHTS:-/fast/sunny/virus-host-mimicry/prostt5}"

# --- Detect input type by extension ---
case "$QUERY" in
  *.fa|*.faa|*.fasta|*.fa.gz|*.faa.gz|*.fasta.gz) USE_PROSTT5=1 ;;
  *)                                              USE_PROSTT5=0 ;;
esac

# --- Sanity checks (fail early with clear messages) ---
if [[ ! -f "$QUERY" ]]; then
  echo "ERROR: query file not found: $QUERY" >&2
  exit 1
fi
if [[ ! -f "${TARGET_DB}" && ! -f "${TARGET_DB}.dbtype" ]]; then
  echo "ERROR: Foldseek DB not found at: $TARGET_DB" >&2
  echo "       Build it once with: foldseek databases PDB $TARGET_DB tmp" >&2
  exit 1
fi
if [[ "$USE_PROSTT5" == "1" && ! -d "$PROSTT5_WEIGHTS" ]]; then
  echo "ERROR: ProstT5 weights dir not found: $PROSTT5_WEIGHTS" >&2
  echo "       Download once with: foldseek databases ProstT5 $PROSTT5_WEIGHTS tmp" >&2
  exit 1
fi

## TODO
# --- Setup ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek

mkdir -p "$OUTPUT"

# Flags that only apply when the query is a FASTA going through ProstT5.
# - --alignment-type 0: 3Di+AA scoring only. ProstT5-built DBs have no Cα, so
#   the foldseek default (TMalign-based, needs Cα) fails at convertalis.
# - --gpu 1: GPU-accelerated ProstT5 inference. With SLURM allocating multiple GPUs
#   (--gres=gpu:N), CUDA_VISIBLE_DEVICES is auto-set to those indices and foldseek
#   splits work across them. No manual export needed.
EXTRA_ARGS=()
if [[ "$USE_PROSTT5" == "1" ]]; then
  EXTRA_ARGS+=(--prostt5-model "$PROSTT5_WEIGHTS")
  EXTRA_ARGS+=(--alignment-type 0)
  EXTRA_ARGS+=(--gpu 1)
  EXTRA_ARGS+=(--prefilter-mode 1)
fi

# echo "Query:      $QUERY  (prostt5=$USE_PROSTT5)"
# echo "Database:   $TARGET_DB"
# echo "Output dir: $OUTPUT"
# echo "Threads:    ${SLURM_CPUS_PER_TASK}"

# --- Run Foldseek ---
srun foldseek easy-search \
  "$QUERY" \
  "$TARGET_DB" \
  "$OUTPUT/hits.m8" \
  "$SCRATCH" \
  "${EXTRA_ARGS[@]}" \
  --format-output "query,target,evalue,bits,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK}"

# echo "Done. Top hits:"
# head -n 10 "$OUTPUT/hits"