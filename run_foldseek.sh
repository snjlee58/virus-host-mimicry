#!/bin/bash
#SBATCH --job-name=foldseek
#SBATCH --partition=compute             # CPU partition — no GPU needed
#SBATCH --cpus-per-task=16              # Foldseek scales well with threads
#SBATCH --mem=16G                       # Bump to 32–64G if searching AFDB
#SBATCH --time=01:00:00
#SBATCH --output=logs/foldseek_%j.log

# Usage:
#   sbatch run_foldseek.sh <query.pdb> <output_dir>
#
# Optional env var:
#   TARGET_DB   path to Foldseek database (default: /fast/sunny/foldseek_dbs/pdb)
#
# Examples:
#   sbatch run_foldseek.sh virprot1_output/top.pdb virprot1_output/foldseek
#   TARGET_DB=/fast/sunny/foldseek_dbs/afdb sbatch run_foldseek.sh top.pdb fs_afdb

set -euo pipefail

# --- Args ---
QUERY_PDB="${1:?Usage: sbatch run_foldseek.sh <query.pdb> <output_dir>}"
OUTPUT="${2:?Usage: sbatch run_foldseek.sh <query.pdb> <output_dir>}"
TARGET_DB="${TARGET_DB:-/fast/databases/foldseek/pdb/pdb100}"

# --- Sanity checks (fail early with clear messages) ---
if [[ ! -f "$QUERY_PDB" ]]; then
  echo "ERROR: query PDB not found: $QUERY_PDB" >&2
  exit 1
fi
if [[ ! -f "${TARGET_DB}" && ! -f "${TARGET_DB}.dbtype" ]]; then
  echo "ERROR: Foldseek DB not found at: $TARGET_DB" >&2
  echo "       Build it once with: foldseek databases PDB $TARGET_DB tmp" >&2
  exit 1
fi

## TODO
# --- Setup ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek

mkdir -p "$OUTPUT"

# echo "Query PDB:  $QUERY_PDB"
# echo "Database:   $TARGET_DB"
# echo "Output dir: $OUTPUT"
# echo "Threads:    ${SLURM_CPUS_PER_TASK}"

# --- Run ---
srun foldseek easy-search \
  "$QUERY_PDB" \
  "$TARGET_DB" \
  "$OUTPUT/hits.m8" \
  "$OUTPUT/tmp" \
  --format-output "query,target,evalue,prob,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK}" ## TODO

# echo "Done. Top hits:"
# head -n 10 "$OUTPUT/hits"