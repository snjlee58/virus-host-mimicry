#!/bin/bash
#SBATCH --job-name=foldseek
#SBATCH --partition=compute
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --time=48:00:00                 # BFVD x AFDB/Proteome TMalign is large; adjust to your partition's limit
# NOTE: SLURM opens these before the job starts, so this dir must already exist
# (mkdir -p it once) and must match $OUTPUT below.
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_reproduction_bfvd/foldseek_%j.log
#SBATCH --error=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_reproduction_bfvd/foldseek_%j.err

# Foldseek structure-vs-structure search for the Fig 1C reproduction
# (Lasso et al. 2021), detecting viral mimicry of host proteins.
#
#   Query  = every viral protein structure in BFVD (a folder of PDB files)
#   Target = AlphaFold/Proteome Foldseek DB (host proteins)
#
# Both sides are real 3D structures, so we align with TMalign
# (--alignment-type 1). 
#
# Usage:
#   sbatch run_foldseek.sh [query] [target_db] [output_dir]
#
# The defaults below are the Fig 1C test, so a bare
#   sbatch run_foldseek.sh
# runs BFVD -> AFDB/Proteome. Positional args override them.

set -euo pipefail

# --- Args (defaults reproduce the Fig 1C test) ---
QUERY="${1:-/fast/sunny/bfvd/2023_02_v2/bfvd}"                     # dir of viral PDB files
TARGET_DB="${2:-/fast/databases/foldseek/afdb_v6/afdb_proteome}"  # prebuilt Foldseek DB
OUTPUT="${3:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_reproduction_bfvd}"

# --- Sanity checks (fail early with clear messages) ---
if [[ ! -e "$QUERY" ]]; then
  echo "ERROR: query not found: $QUERY" >&2
  exit 1
fi
if [[ ! -f "${TARGET_DB}" && ! -f "${TARGET_DB}.dbtype" ]]; then
  echo "ERROR: Foldseek target DB not found at: $TARGET_DB" >&2
  exit 1
fi

# --- Setup (uncomment whichever your cluster needs to expose `foldseek`) ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek
# module load foldseek

mkdir -p "$OUTPUT"
TMP="${SCRATCH:-$OUTPUT/tmp}"
mkdir -p "$TMP"

echo "Query:     $QUERY"
echo "Target DB: $TARGET_DB"
echo "Output:    $OUTPUT/hits.m8"
echo "Threads:   ${SLURM_CPUS_PER_TASK:-64}"

# --- Run Foldseek ---
# --alignment-type 1 : TMalign structural superposition (query & target are both structures).
# taxid/taxname/taxlineage feed the downstream host-division classification; they require the
# target DB to carry taxonomy (the AFDB Foldseek DBs do).
srun foldseek easy-search \
  "$QUERY" \
  "$TARGET_DB" \
  "$OUTPUT/hits.m8" \
  "$TMP" \
  --alignment-type 1 \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,qcov,tcov,qlen,tlen,alnlen,fident,lddt,prob,evalue,bits,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK:-64}"

echo "Done. Top hits:"
head -n 10 "$OUTPUT/hits.m8"
