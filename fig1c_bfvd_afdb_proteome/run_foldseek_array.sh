#!/bin/bash
#SBATCH --job-name=fs_bfvd2_arr
#SBATCH --partition=compute
#SBATCH --exclude=hulk                 # keep hulk (the 500-bldg login/gateway) responsive; drop this to use its 128 cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32             # a 64-core node runs 2 tasks; super011 (96c) runs 3
#SBATCH --mem=64G                      # mostly the afdb_proteome index; bump if a task OOMs
#SBATCH --time=12:00:00                # per task; one ~20k chunk is ~1-3h, this is generous headroom
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_proteome/logs/fs_%A_%a.log
#SBATCH --error=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_proteome/logs/fs_%A_%a.err

# SLURM-array Foldseek: each task searches ONE chunk of BFVD2 structures
# (from prepare_chunks.sh) against afdb_proteome with TMalign.
#
# Concurrency is set at submit time via --array=...%N (N tasks at once). With
# --cpus-per-task=32 and hulk excluded, N=10 uses ~320 of the ~480 non-hulk
# cores, leaving headroom for other users and the running 351k job.
#
# Workflow:
#   bash prepare_chunks.sh                      # prints K
#   sbatch --array=0-<K-1>%10 run_foldseek_array.sh
#   bash merge_hits.sh                          # after all tasks finish

set -euo pipefail

OUT="${1:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_proteome}"
TARGET_DB="${2:-/fast/databases/foldseek/afdb_v6/afdb_proteome}"

# --- Setup (uncomment whichever your cluster needs to expose `foldseek`) ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek
# module load foldseek

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: submit as an array: sbatch --array=0-<K-1>%10 $0" >&2; exit 1; }

TASK=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
LIST="$OUT/chunks/list_$TASK"
[[ -f "$LIST" ]] || { echo "ERROR: chunk list not found: $LIST" >&2; exit 1; }
if [[ ! -f "${TARGET_DB}" && ! -f "${TARGET_DB}.dbtype" ]]; then
  echo "ERROR: Foldseek target DB not found: $TARGET_DB" >&2; exit 1
fi

# Per-task scratch (node-local); auto-cleaned on exit so concurrent tasks never collide.
TMP="${SCRATCH:-/tmp}/fs_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
QDIR="$TMP/query"
mkdir -p "$QDIR" "$TMP/tmp"
trap 'rm -rf "$TMP"' EXIT

# Symlink this chunk's PDBs into the per-task query dir (targets are absolute /fast paths).
while IFS= read -r pdb; do
  [[ -n "$pdb" ]] && ln -s "$pdb" "$QDIR/$(basename "$pdb")"
done < "$LIST"
n=$(find "$QDIR" -maxdepth 1 -name '*.pdb' | wc -l | tr -d ' ')

OUT_M8="$OUT/hits/hits_${TASK}.m8"
echo "task $TASK: $n structures on $(hostname), ${SLURM_CPUS_PER_TASK} threads -> $OUT_M8"

start=$(date +%s)
srun foldseek easy-search \
  "$QDIR" \
  "$TARGET_DB" \
  "$OUT_M8" \
  "$TMP/tmp" \
  --alignment-type 1 --tmscore-threshold 0.5 \
  --cov-mode 2 -c 0.7 \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,qcov,tcov,qlen,tlen,alnlen,fident,lddt,prob,evalue,bits,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK}"
echo "task $TASK done in $(( $(date +%s) - start )) s"
