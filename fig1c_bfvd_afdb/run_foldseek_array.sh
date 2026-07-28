#!/bin/bash
#SBATCH --job-name=fs_bfvd2_afdb
#SBATCH --partition=compute
#SBATCH --exclude=hulk                 # keep the login node responsive; drop to use its 128 cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G                     # generous for the afdb50 cluster DB; check the pilot's real peak and lower for better packing
#SBATCH --time=24:00:00                # expansion search is fast, but PILOT chunk 0 to confirm runtime/output size
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb/logs/fs_%A_%a.log
#SBATCH --error=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb/logs/fs_%A_%a.err

# SLURM-array Foldseek EXPANSION (cluster) search: each task searches ONE chunk of
# BFVD2 structures against the AFDB50 cluster DB. --cluster-search 1 aligns the
# query to cluster representatives and then expands each hit to all its members
# (the fast, intended way to search the full AFDB).
#
# Uses the low-level createdb -> search -> convertalis flow (NOT easy-search),
# because the target is a cluster-search DB and `search` needs a real query DB.
# The align params are the same ones from the proteome run and all apply to
# `foldseek search`:  --alignment-type 1 (TMalign), --tmscore-threshold 0.5,
# --cov-mode 2 -c 0.7.
#
# Workflow:
#   bash prepare_chunks.sh                       # prints K
#   sbatch --array=0 run_foldseek_array.sh       # PILOT one chunk (runtime/mem/output size)
#   sbatch --array=1-<K-1>%8 run_foldseek_array.sh
#   bash merge_hits.sh

set -euo pipefail

OUT="${1:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb}"
TARGET_DB="${2:-/fast/databases/foldseek/afdb_v6/afdb50}"       # cluster-search REP DB; foldseek auto-finds <TARGET_DB>_seq (members) + _clu
MEMBER_DB="${TARGET_DB}_seq"                                     # expanded hits are member-space -> convertalis resolves names/taxonomy here

# --- Setup (uncomment whichever your cluster needs to expose `foldseek`) ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek
# module load foldseek

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: submit as an array: sbatch --array=0-<K-1>%8 $0" >&2; exit 1; }

TASK=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
LIST="$OUT/chunks/list_$TASK"
[[ -f "$LIST" ]] || { echo "ERROR: chunk list not found: $LIST" >&2; exit 1; }
[[ -f "${TARGET_DB}.dbtype" ]] || { echo "ERROR: cluster-search rep DB not found: ${TARGET_DB}.dbtype" >&2; exit 1; }
[[ -f "${MEMBER_DB}.dbtype" ]] || { echo "ERROR: member DB not found: ${MEMBER_DB}.dbtype (needed for expansion)" >&2; exit 1; }

# Per-task scratch (node-local); auto-cleaned so concurrent tasks never collide.
TMP="${SCRATCH:-/tmp}/fs_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
QDIR="$TMP/query_pdbs"
mkdir -p "$QDIR" "$TMP/searchtmp"
trap 'rm -rf "$TMP"' EXIT

while IFS= read -r pdb; do
  [[ -n "$pdb" ]] && ln -s "$pdb" "$QDIR/$(basename "$pdb")"
done < "$LIST"
n=$(find "$QDIR" -maxdepth 1 -name '*.pdb' | wc -l | tr -d ' ')

QUERYDB="$TMP/queryDB"
RESULT="$TMP/result"
OUT_M8="$OUT/hits/hits_${TASK}.m8"
echo "task $TASK: $n structures on $(hostname), ${SLURM_CPUS_PER_TASK} threads -> $OUT_M8"

start=$(date +%s)

# 1) build the query DB from this chunk's structures
foldseek createdb "$QDIR" "$QUERYDB"

# 2) expansion (cluster) search against the AFDB50 cluster DB
foldseek search "$QUERYDB" "$TARGET_DB" "$RESULT" "$TMP/searchtmp" \
  --cluster-search 1 \
  --alignment-type 1 --tmscore-threshold 0.5 \
  --cov-mode 2 -c 0.7 \
  --threads "${SLURM_CPUS_PER_TASK}"

# 3) format the results against the MEMBER DB (expanded hits are member-space).
#    afdb50 carries taxonomy (afdb50_taxonomy / _mapping), so the taxid columns fill in.
foldseek convertalis "$QUERYDB" "$MEMBER_DB" "$RESULT" "$OUT_M8" \
  --format-output "query,target,alntmscore,qcov,tcov,evalue,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK}"

echo "task $TASK done in $(( $(date +%s) - start )) s ($(wc -l < "$OUT_M8") rows)"
