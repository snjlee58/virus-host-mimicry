#!/bin/bash
#SBATCH --job-name=fs_afdb50seq
#SBATCH --partition=gpu
#SBATCH --nodelist=devlss001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --mem=480G
#SBATCH --gres=gpu:6
#SBATCH --time=24:00:00
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50/slurm_afdb50seq_%j.out

set -euo pipefail

# ---------------------------------------------------------------------------
# fig1c arm 2: all 241M AFDB members, searched flat.
#
# This one genuinely needs most of the node (~170GB padded across 6 cards),
# so expect to queue until devlss001 drains. Submit it and let it wait --
# arm 1 is not blocked by it.
#
# --db-load-mode 2 mmaps the target instead of loading it: the Ca data for
# 241M members would need ~680GB of host RAM, and the node has 515GB.
# ---------------------------------------------------------------------------

DBDIR=/fast/sunny/databases/afdb_v6
OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50

QDB=$OUT/bfvd1000DB
MEMBERS=$DBDIR/afdb50_seq_pad

THREADS=64

cd "$OUT"

if [ ! -e "${QDB}.dbtype" ]; then
  echo "ERROR: query DB missing. Run build_bfvd1000_querydb.sh first." >&2
  exit 1
fi

echo "node: $(hostname)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv
echo

FMT="query,target,fident,evalue,qlen,tlen,qtmscore,ttmscore,qcov,tcov"

# --max-seqs is 4.5x arm 1's so both reach comparable cluster depth rather
# than this arm burning its budget on near-duplicate members.

/usr/bin/time -v -o afdb50_seq.time \
foldseek search "$QDB" "$MEMBERS" res_afdb50_seq tmp_afdb50_seq \
  --gpu 1 \
  --alignment-type 2 \
  -e 10 \
  --max-seqs 45000 \
  --db-load-mode 2 \
  --threads "$THREADS" \
  -a

foldseek convertalis "$QDB" "$MEMBERS" res_afdb50_seq afdb50_seq.m8 \
  --threads "$THREADS" --format-output "$FMT"

echo
echo "rows: $(wc -l < afdb50_seq.m8)"
echo "elapsed: $(grep -m1 'Elapsed' afdb50_seq.time | sed 's/.*): //')"
