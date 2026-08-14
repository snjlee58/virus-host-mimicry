#!/bin/bash
#SBATCH --job-name=bfvd2_qdb
#SBATCH --partition=compute
#SBATCH --exclude=hulk
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_family_subset/logs/bfvd2_qdb_%j.out
#
# ONE-TIME: build a Foldseek query DB over ALL of BFVD2 (1,218,808 structures).
#
# WHY A WHOLE-DB BUILD AND NOT PER-TASK createdb
# The array tasks each need ~76k of these structures. Passing 76k paths on the
# command line the way build_bfvd1000_querydb.sh does for 1,000 blows past
# ARG_MAX (~76k x 80 chars = 6MB vs a 2MB limit), and symlinking chunks into
# per-task directories costs 1.2M inodes on a /fast that sits at ~98% full.
# Building once and slicing with createsubdb avoids both: each task's slice is
# an .index over data this DB already holds.
#
# Usage:  sbatch build_bfvd2_querydb.sh      (then foldseek_family_gpu_array.sh)
set -euo pipefail

export PATH=$HOME/bin/foldseek/bin:$PATH

PDBDIR=/fast/sunny/bfvd2/bfvd2_pdb_only
OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_family_subset
QDB=$OUT/bfvd2DB
THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "$OUT/logs"

if [ -e "${QDB}.dbtype" ]; then
  echo "Query DB already exists: $QDB  ($(wc -l < "${QDB}.lookup") entries)"
  echo "Delete bfvd2DB* first if you really want to rebuild."
  exit 0
fi

[ -d "$PDBDIR" ] || { echo "ERROR: $PDBDIR not found" >&2; exit 1; }

echo "node $(hostname)  threads $THREADS"
echo "building query DB from $PDBDIR"
/usr/bin/time -v -o "$OUT/logs/bfvd2_qdb_time.txt" \
  foldseek createdb "$PDBDIR" "$QDB" --threads "$THREADS"

N=$(wc -l < "${QDB}.lookup")
echo
echo "query DB entries: $N"
echo "components:"; ls -la "${QDB}"* | awk '{print "  " $5 "\t" $9}'
grep -E 'Elapsed \(wall|Maximum resident' "$OUT/logs/bfvd2_qdb_time.txt"
echo
echo "Next: sbatch --array=0 foldseek_family_gpu_array.sh   # pilot one chunk"
