#!/bin/bash
#SBATCH --job-name=foldseek
#SBATCH --partition=gpu                 # GPU partition (needed for ProstT5 inference)
#SBATCH --gres=gpu:1                    # Request 1 GPU; PDB-only runs leave it idle
#SBATCH --cpus-per-task=16              # Foldseek search step still uses CPU threads
#SBATCH --mem=32G                       # Bump higher (e.g. 64G) if searching AFDB
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
TARGET_DB="${TARGET_DB:-/fast/databases/foldseek/pdb/pdb100}"
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
# - --gpu 1: enables GPU-accelerated ProstT5 inference (~100-1000x vs CPU).
EXTRA_ARGS=()
if [[ "$USE_PROSTT5" == "1" ]]; then
  EXTRA_ARGS+=(--prostt5-model "$PROSTT5_WEIGHTS")
  EXTRA_ARGS+=(--alignment-type 0)
  # EXTRA_ARGS+=(--gpu 1)   # re-enable once GPU access is sorted out
fi

# echo "Query:      $QUERY  (prostt5=$USE_PROSTT5)"
# echo "Database:   $TARGET_DB"
# echo "Output dir: $OUTPUT"
# echo "Threads:    ${SLURM_CPUS_PER_TASK}"

# --- Run ---
# Use $SCRATCH (node-local NVMe at /mnt/scratch/$USER/$JOB_ID) for foldseek's
# tmp dir — MMseqs2-style intermediate I/O on /fast saturates the network FS.
# $SCRATCH is auto-created by SLURM and auto-deleted at job end.
srun foldseek easy-search \
  "$QUERY" \
  "$TARGET_DB" \
  "$OUTPUT/hits.m8" \
  "$SCRATCH" \
  "${EXTRA_ARGS[@]}" \
  --format-output "query,target,evalue,prob,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK}"

# echo "Done. Top hits:"
# head -n 10 "$OUTPUT/hits"