#!/bin/bash
#SBATCH --job-name=fs_bfvd2_afdb
#SBATCH --partition=compute
#SBATCH --exclude=hulk                 # keep the login node responsive; drop to use its 128 cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread           # bind to PHYSICAL cores -> avoids the 2x oversubscription seen earlier
#SBATCH --mem=256G                     # TMalign needs the Ca coords over a 241M target; tune to your colleague's numbers
#SBATCH --time=0                       # no limit (compute MaxTime is 'infinite') -> won't be walltime-killed
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb/logs/fs_%A_%a.log
#SBATCH --error=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb/logs/fs_%A_%a.err

# SLURM-array Foldseek FLAT search: each task 3Di+AA-searches ONE chunk of BFVD2
# structures against the AFDB50 cluster REPS (afdb50 = 66.7M) -- NO --cluster-search.
# Hits are query-vs-rep; expand rep->cluster-members downstream (afdb50_clu + afdb50_mapping)
# to assign member taxonomy. NOTE: this is the CLUSTERED-reps route -- a fast, finishable
# first pass; the un-clustered afdb50_seq (241M) is what Gorka ultimately wants, revisit later.
#
#   --alignment-type 2        3Di+AA (fast default; TM-score still reported as qtmscore/
#                             ttmscore, so apply the TM>=0.5 cutoff downstream)
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
TARGET_DB="${2:-/fast/databases/foldseek/afdb_v6/afdb50}"       # AFDB50 cluster reps (66.7M); expand hits -> members downstream

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
echo "task $TASK: $n structures on $(hostname), ${SLURM_CPUS_PER_TASK} threads -> $OUT_M8"

start=$(date +%s)

# 1) build the query DB from this chunk's structures
foldseek createdb "$QDIR" "$QUERYDB"

# 2) flat 3Di+AA search against the AFDB50 cluster reps (no --cluster-search)
foldseek search "$QUERYDB" "$TARGET_DB" "$RESULT" "$TMP/searchtmp" \
  --alignment-type 2 \
  --cov-mode 2 -c 0.7 \
  --threads "${SLURM_CPUS_PER_TASK}"

# 3) format the results (structural columns; TM-score as qtmscore/ttmscore). Taxonomy is
#    NOT here -- assign it in the downstream rep->member expansion via the target IDs.
foldseek convertalis "$QUERYDB" "$TARGET_DB" "$RESULT" "$OUT_M8" \
  --format-output "query,target,fident,evalue,qlen,tlen,qtmscore,ttmscore,qcov,tcov" \
  --threads "${SLURM_CPUS_PER_TASK}"

echo "task $TASK done in $(( $(date +%s) - start )) s ($(wc -l < "$OUT_M8") rows)"
