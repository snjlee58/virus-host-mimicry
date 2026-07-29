#!/bin/bash
#SBATCH --job-name=fs_bfvd2_afdb
#SBATCH --partition=compute
#SBATCH --exclude=hulk                 # keep the login node responsive; drop to use its 128 cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread           # bind to PHYSICAL cores -> avoids the 2x oversubscription seen earlier
#SBATCH --mem=256G                     # TMalign needs the Ca coords over a 241M target; tune to your colleague's numbers
#SBATCH --time=48:00:00                # flat TMalign vs 241M per big chunk is heavy; PILOT chunk 0 first
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb/logs/fs_%A_%a.log
#SBATCH --error=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb/logs/fs_%A_%a.err

# SLURM-array Foldseek FLAT search: each task TMaligns ONE chunk of BFVD2 structures
# against the FULL, un-clustered AFDB (afdb50_seq = 241M members) -- NO --cluster-search,
# so every hit is a direct query-vs-member alignment carrying that member's OWN taxonomy
# (afdb50_seq's mapping covers all 241M members). This is Gorka's un-clustered design;
# it also replaces the separate afdb_proteome run.
#
#   --alignment-type 1        TMalign (mimicry-faithful structural superposition)
#   --tmscore-threshold 0.5   keep only "significant structural relationships"
#   --cov-mode 2 -c 0.7       require >=70% query coverage
#   --max-seqs 4000           don't truncate the neighbour set at the default 1000 -- a common
#                             viral fold can exceed 1000 homologs in 241M, which would bias the
#                             enrichment counts. Confirm from the pilot's hits-per-query spread.
#   --split-memory-limit      cap target-split RAM: foldseek sizes splits from the NODE's total
#                             RAM, not the SLURM cgroup, so without this it can exceed --mem and OOM.
#
# NOTE on chunking: each chunk re-scans the full 241M target, so keep chunks LARGE and FEW
# (see prepare_chunks.sh CHUNK_SIZE ~250k -> ~5 chunks) to avoid redundant target scans.
#
# Workflow:
#   bash prepare_chunks.sh                        # prints K (few large chunks)
#   sbatch --array=0 run_foldseek_array.sh        # PILOT one chunk: time / mem / hits-per-query / taxids
#   sbatch --array=1-<K-1>%3 run_foldseek_array.sh
#   bash merge_hits.sh

set -euo pipefail

OUT="${1:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb}"
TARGET_DB="${2:-/fast/databases/foldseek/afdb_v6/afdb50_seq}"   # FULL un-clustered AFDB (241M), member-level taxonomy
SPLIT_MEM="${SPLIT_MEM:-220G}"                                  # keep < --mem

# --- Setup (uncomment whichever your cluster needs to expose `foldseek`) ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek
# module load foldseek

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: submit as an array: sbatch --array=0-<K-1>%3 $0" >&2; exit 1; }

TASK=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
LIST="$OUT/chunks/list_$TASK"
[[ -f "$LIST" ]] || { echo "ERROR: chunk list not found: $LIST" >&2; exit 1; }
[[ -f "${TARGET_DB}.dbtype" ]] || { echo "ERROR: target DB not found: ${TARGET_DB}.dbtype" >&2; exit 1; }

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
echo "task $TASK: $n structures on $(hostname), ${SLURM_CPUS_PER_TASK} threads, split-mem $SPLIT_MEM -> $OUT_M8"

start=$(date +%s)

# 1) build the query DB from this chunk's structures
foldseek createdb "$QDIR" "$QUERYDB"

# 2) flat TMalign search against the full un-clustered AFDB (no --cluster-search)
foldseek search "$QUERYDB" "$TARGET_DB" "$RESULT" "$TMP/searchtmp" \
  --alignment-type 1 --tmscore-threshold 0.5 \
  --cov-mode 2 -c 0.7 \
  --max-seqs 4000 \
  --split-memory-limit "$SPLIT_MEM" \
  --threads "${SLURM_CPUS_PER_TASK}"

# 3) format the results. afdb50_seq carries member-level taxonomy, so taxid/taxname/
#    taxlineage populate for every hit.
foldseek convertalis "$QUERYDB" "$TARGET_DB" "$RESULT" "$OUT_M8" \
  --format-output "query,target,alntmscore,qcov,tcov,evalue,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK}"

echo "task $TASK done in $(( $(date +%s) - start )) s ($(wc -l < "$OUT_M8") rows)"
