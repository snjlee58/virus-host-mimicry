#!/bin/bash
#SBATCH --job-name=fs_bfvd2_afdb
#SBATCH --partition=compute
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --time=120:00:00                # extended BFVD2 x afdb_proteome TMalign is large; `compute` has no wall limit, adjust as needed
# NOTE: SLURM opens these before the job starts, so this dir must already exist (mkdir -p it once).
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_proteome/foldseek_%j.log
#SBATCH --error=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_proteome/foldseek_%j.err

# Fig 1C reproduction (Lasso et al. 2021): BFVD2 viral structures vs AFDB/Proteome.
#
#   Query  = every viral protein structure in BFVD2 (a folder of PDB files)
#   Target = AlphaFold/Proteome Foldseek DB (candidate host proteins)
#
# Both sides are real structures -> TMalign (--alignment-type 1); no ProstT5/GPU.
#   --tmscore-threshold 0.5 : keep only TM-score >= 0.5 hits = "significant
#                             structural relationship" (the paper's SAS < 2.5 A analog)
#   --cov-mode 2 -c 0.7     : require the alignment to cover >= 70% of the viral query
# taxid/taxname/taxlineage give each AFDB target's taxonomy for the downstream
# enrichment (afdb_proteome carries taxonomy, so this is free).
#
# Usage:
#   sbatch run_foldseek.sh [query] [target_db] [output_dir]

set -euo pipefail

QUERY="${1:-/fast/sunny/bfvd2/bfvd2_pdb_only}"                       # dir of BFVD2 viral PDBs
TARGET_DB="${2:-/fast/databases/foldseek/afdb_v6/afdb_proteome}"     # prebuilt Foldseek DB
OUTPUT="${3:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_proteome}"

# --- Setup (uncomment whichever your cluster needs to expose `foldseek`) ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek
# module load foldseek

# --- Sanity checks ---
[[ -e "$QUERY" ]] || { echo "ERROR: query not found: $QUERY" >&2; exit 1; }
if [[ ! -f "${TARGET_DB}" && ! -f "${TARGET_DB}.dbtype" ]]; then
  echo "ERROR: Foldseek target DB not found: $TARGET_DB" >&2; exit 1
fi

mkdir -p "$OUTPUT"
TMP="${SCRATCH:-$OUTPUT/tmp}"
mkdir -p "$TMP"

echo "Query:     $QUERY"
echo "Target DB: $TARGET_DB"
echo "Output:    $OUTPUT/hits.m8"
echo "Threads:   ${SLURM_CPUS_PER_TASK:-64}"

# --- Run Foldseek ---
srun foldseek easy-search \
  "$QUERY" \
  "$TARGET_DB" \
  "$OUTPUT/hits.m8" \
  "$TMP" \
  --alignment-type 1 --tmscore-threshold 0.5 \
  --cov-mode 2 -c 0.7 \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,qcov,tcov,qlen,tlen,alnlen,fident,lddt,prob,evalue,bits,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK:-64}"

echo "Done. Top hits:"
head -n 10 "$OUTPUT/hits.m8"
