#!/bin/bash
#SBATCH --job-name=fs_family_gpu
#SBATCH --partition=gpu
#SBATCH --exclude=hulk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --hint=nomultithread
#SBATCH --mem=90G
#SBATCH --gres=gpu:1
#SBATCH --time=36:00:00
#SBATCH --array=0-15
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_family_subset/logs/gpu_%A_%a.out
#
# fig1c: ALL of BFVD2 (1,218,808 structures) vs the per-family AFDB subset,
# split 16 ways on GPU.
#
# TARGET
# /fast/sunny/databases/afdb_v6_family_subset/family_pad -- 22,971,764 structures,
# one reference proteome per taxonomic family (1,912 families), sliced out of
# afdb50_seq by build_reference_db/build_subdb.sh and then GPU-padded.
# Padding verified not to disturb _h: target names resolve as AF-<acc>-F1-model_v6
# and 99.89% of hits map back to a family via _meta/acc2label.tsv.
#
# MEASURED (1,000-query pilot, job 550525, one L40S, 16 threads)
#   wall 16m26s   -- prefilter 10m41s, structurealign 5m42s, convertalis 1m44s
#   RSS  98.9 GB  -- mmap; --alignment-type 2 never touches the Ca data
#   rows 2,203,208 unfiltered
# => ~0.0164 min/query. 1,218,808 / 16 = 76,176 queries per task ~= 21 h.
#    --time=36:00:00 leaves headroom for the createsubdb slice and a slow card.
#
# For comparison the same 1,000 queries under CPU TMalign (--alignment-type 1)
# took 30m41s and 27.4 CPU-hours; the full run that way is ~65 h on 512 cores.
# This arm uses --alignment-type 2 with no score thresholds, matching
# foldseek_afdb50.sh and foldseek_afdb50_seq.sh so all three arms are comparable.
# Filter post-hoc on qtmscore/qcov. NOTE: TMalign is more sensitive -- on the
# pilot, 99 of 438 queries it called significant were missed here.
#
# SIZING
#   --cpus-per-task=16  the GPU prefilter dominates; the pilot used 16 and the
#                       alignment stage was only a third of the wall clock.
#   --mem=90G           pilot peaked at 98.9 GB RSS, but that is mmap page cache
#                       (file-backed, reclaimable) -- anonymous use is far lower.
#   --db-load-mode 2    mmap, so tasks sharing a node share one page-cache copy
#                       of the padded _ss instead of each reading its own.
#   devbox/devrtx NOT excluded: the padded _ss for 23M entries is ~6.4 GB
#                       (afdb50's 66.7M reps need ~18.5 GB), which fits a 24 GB
#                       card comfortably. More cards = less queueing.
#
# PREREQUISITES
#   bash  build_reference_db/build_subdb.sh ...            # the target subset
#   foldseek makepaddedseqdb <subset>/family <subset>/family_pad
#   sbatch build_bfvd2_querydb.sh                          # the query DB
#
# USAGE
#   sbatch --array=0    foldseek_family_gpu_array.sh   # pilot ONE chunk first
#   sbatch --array=1-15 foldseek_family_gpu_array.sh   # then the rest
#   bash   merge_family_gpu.sh
set -euo pipefail

export PATH=$HOME/bin/foldseek/bin:$PATH

OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_family_subset
QDB=$OUT/bfvd2DB
TARGET=/fast/sunny/databases/afdb_v6_family_subset/family_pad

NCHUNKS=16         # must match the --array range in the header
MAXSEQS=10000      # matches the afdb50 arms; the default 1000 saturated for ~10%
EVALUE=10          # of pilot queries and would bias enrichment counts downward

THREADS="${SLURM_CPUS_PER_TASK:-16}"
TASK=$(printf '%02d' "${SLURM_ARRAY_TASK_ID:?submit with sbatch --array=..., not bash}")

[ "$SLURM_ARRAY_TASK_ID" -lt "$NCHUNKS" ] || {
  echo "ERROR: task $SLURM_ARRAY_TASK_ID >= NCHUNKS=$NCHUNKS" >&2; exit 1; }

# --- preflight: every one of these has bitten a previous attempt -------------
[ -e "${QDB}.dbtype" ] || {
  echo "ERROR: query DB missing: $QDB. Run: sbatch build_bfvd2_querydb.sh" >&2; exit 1; }
[ -e "${TARGET}_ss.dbtype" ] || {
  echo "ERROR: padded target missing: $TARGET" >&2
  echo "       foldseek makepaddedseqdb ${TARGET%_pad} $TARGET" >&2; exit 1; }
PADFLAG=$(od -An -tu1 "${TARGET}_ss.dbtype" | tr -s ' ')
[ "$PADFLAG" = " 0 0 8 0" ] || {
  echo "ERROR: $TARGET is not GPU-padded (_ss.dbtype =$PADFLAG, want ' 0 0 8 0')." >&2; exit 1; }

mkdir -p "$OUT/hits" "$OUT/logs"

# --- node-local scratch -----------------------------------------------------
# /fast sits at ~98% full and local NVMe beats it for prefilter spill anyway.
TMP="/tmp/fs_fam_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMP/search"
trap 'rm -rf "$TMP"' EXIT

# --- this task's slice of the query DB --------------------------------------
# Slice by DB key rather than by file path: 76k paths on a command line exceed
# ARG_MAX. Round-robin (key % NCHUNKS) rather than contiguous blocks, because the
# DB was built from a sorted directory listing and long structures cluster.
cut -f1 "${QDB}.index" | awk -v n="$NCHUNKS" -v t="$SLURM_ARRAY_TASK_ID" \
  'NR % n == t' > "$TMP/keys.txt"
NQ=$(wc -l < "$TMP/keys.txt")
[ "$NQ" -gt 0 ] || { echo "ERROR: empty chunk $TASK" >&2; exit 1; }

QCHUNK="$TMP/qdb"
for suf in "" _ss _ca _h; do
  [ -e "${QDB}${suf}.index" ] || continue
  foldseek createsubdb "$TMP/keys.txt" "${QDB}${suf}" "${QCHUNK}${suf}" --subdb-mode 1
  [ -s "${QCHUNK}${suf}.index" ] || { echo "ERROR: query slice ${suf:-base} empty" >&2; exit 1; }
done

echo "task $TASK  node $(hostname)  $NQ queries  $THREADS threads"
nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv
echo "target: $TARGET ($(wc -l < "${TARGET}_ss.index") structures)"
echo

FMT="query,target,fident,evalue,qlen,tlen,qtmscore,ttmscore,qcov,tcov"
M8="$OUT/hits/family_gpu_${TASK}.m8"

/usr/bin/time -v -o "$OUT/logs/time_${TASK}.txt" \
foldseek search "$QCHUNK" "$TARGET" "$TMP/res" "$TMP/search" \
  --gpu 1 \
  --alignment-type 2 \
  -e "$EVALUE" \
  --max-seqs "$MAXSEQS" \
  --db-load-mode 2 \
  --threads "$THREADS" \
  -a

foldseek convertalis "$QCHUNK" "$TARGET" "$TMP/res" "$M8" \
  --db-load-mode 2 --threads "$THREADS" --format-output "$FMT"

echo
echo "task $TASK done: $(wc -l < "$M8") rows -> $M8"
grep -E 'Elapsed \(wall|Maximum resident' "$OUT/logs/time_${TASK}.txt"
