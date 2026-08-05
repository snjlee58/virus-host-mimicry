#!/bin/bash
#SBATCH --job-name=fs_afdb50
#SBATCH --partition=gpu
#SBATCH --nodelist=devlss001,devlss002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --mem=200G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50/slurm_afdb50_%j.out

set -euo pipefail

# ---------------------------------------------------------------------------
# fig1c arm 1: 54M AFDB50 representatives + cluster expansion.
#
# Sized to fit alongside other jobs on a partially-allocated node:
#   1 GPU  (~40GB padded DB, needs a 48GB L40S -- hence the nodelist;
#            devbox/devrtx have 24GB cards and cannot hold it)
#   32 physical cores of 64
#   200GB  (AFDB50 with Ca coordinates is ~151GB)
#
# This is the arm that gives you the production timing estimate.
# ---------------------------------------------------------------------------

DBDIR=/fast/sunny/databases/afdb_v6
OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50

QDB=$OUT/bfvd1000DB
REPS=$DBDIR/afdb50_pad

THREADS=32

cd "$OUT"

if [ ! -e "${QDB}.dbtype" ]; then
  echo "ERROR: query DB missing. Run build_bfvd1000_querydb.sh first." >&2
  exit 1
fi

echo "node: $(hostname)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv
echo

FMT="query,target,fident,evalue,qlen,tlen,qtmscore,ttmscore,qcov,tcov"

/usr/bin/time -v -o afdb50.time \
foldseek search "$QDB" "$REPS" res_afdb50 tmp_afdb50 \
  --gpu 1 \
  --cluster-search 1 \
  --alignment-type 2 \
  -e 10 \
  --max-seqs 10000 \
  --threads "$THREADS" \
  -a

foldseek convertalis "$QDB" "$REPS" res_afdb50 afdb50.m8 \
  --threads "$THREADS" --format-output "$FMT"

echo
echo "rows: $(wc -l < afdb50.m8)"
echo "elapsed: $(grep -m1 'Elapsed' afdb50.time | sed 's/.*): //')"
echo
echo "CHECK qtmscore/ttmscore (cols 7,8) are populated:"
head -3 afdb50.m8 | column -t
