#!/bin/bash
#SBATCH --job-name=fs_bench_1k
#SBATCH --partition=gpu
#SBATCH --nodelist=devlss001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --mem=480G
#SBATCH --gres=gpu:6
#SBATCH --time=24:00:00
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50/slurm_%j.out

set -euo pipefail

# ---------------------------------------------------------------------------
# fig1c: two ways to search BFVD against AlphaFold DB
#
#   afdb50      = 54M representatives + cluster expansion. ~40GB padded,
#                 fits one 48GB L40S.
#   afdb50_seq  = all 241M members directly. ~170GB padded, needs 6 cards.
#                 Uses --db-load-mode 2 because the Ca data needed for
#                 qtmscore/ttmscore would otherwise exceed 503GB host RAM.
#
# devlss001: 128 CPUs (64 physical cores), 515GB, 6x L40S, 7TB scratch.
# --hint=nomultithread + 64 cpus-per-task = all 64 physical cores, one
# thread each. Hyperthreads add little to Foldseek's SIMD alignment work.
# Alternative is devlss002 (8x L40S, fewer cores) -- more GPUs, but since
# --gpu only accelerates the prefilter, cores drive the stage that
# actually dominates. devlss001 is the better trade.
#
# Prerequisite: run build_bfvd1000_querydb.sh first.
# ---------------------------------------------------------------------------

DBDIR=/fast/sunny/databases/afdb_v6
OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50

QDB=$OUT/bfvd1000DB               # prebuilt, persistent
REPS=$DBDIR/afdb50_pad            # padded AFDB50 representatives
MEMBERS=$DBDIR/afdb50_seq_pad     # padded full member DB

THREADS=64

cd "$OUT"

if [ ! -e "${QDB}.dbtype" ]; then
  echo "ERROR: query DB missing. Run build_bfvd1000_querydb.sh first." >&2
  exit 1
fi

nvidia-smi --query-gpu=index,name,memory.total --format=csv
echo

FMT="query,target,fident,evalue,qlen,tlen,qtmscore,ttmscore,qcov,tcov"

# --- afdb50: representatives, then cluster expansion ------------------------
# One GPU. In production this is the arm you replicate 6x across the node.

export CUDA_VISIBLE_DEVICES=0

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

# --- afdb50_seq: flat search against all members ----------------------------
# All 6 GPUs. --max-seqs is 4.5x the afdb50 run's so both reach comparable
# cluster depth rather than afdb50_seq burning its budget on near-duplicates.

export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5

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

# --- compare ----------------------------------------------------------------

cut -f1,2 afdb50.m8     | sort -u > afdb50.pairs
cut -f1,2 afdb50_seq.m8 | sort -u > afdb50_seq.pairs

# hits afdb50_seq found that the faster afdb50 route did not
comm -13 afdb50.pairs afdb50_seq.pairs > missed_by_afdb50.pairs
# hits only afdb50 found; expect ~0, a large number means misconfiguration
comm -23 afdb50.pairs afdb50_seq.pairs > only_in_afdb50.pairs

{
  echo "node:                   $(hostname)"
  echo "queries in DB:          $(wc -l < "${QDB}.lookup")"
  echo
  echo "afdb50 pairs:           $(wc -l < afdb50.pairs)"
  echo "afdb50_seq pairs:       $(wc -l < afdb50_seq.pairs)"
  echo "missed by afdb50:       $(wc -l < missed_by_afdb50.pairs)"
  echo "only in afdb50:         $(wc -l < only_in_afdb50.pairs)"
  echo "queries affected:       $(cut -f1 missed_by_afdb50.pairs | sort -u | wc -l)"
  echo
  echo "afdb50 elapsed:         $(grep -m1 'Elapsed' afdb50.time     | sed 's/.*): //' || true)"
  echo "afdb50_seq elapsed:     $(grep -m1 'Elapsed' afdb50_seq.time | sed 's/.*): //' || true)"
  echo "afdb50 peak RSS:        $(grep -m1 'Maximum resident' afdb50.time     | sed 's/.*: //' || true)"
  echo "afdb50_seq peak RSS:    $(grep -m1 'Maximum resident' afdb50_seq.time | sed 's/.*: //' || true)"
} > summary.txt

cat summary.txt

# Uncomment once you trust the run -- these get large:
# rm -rf tmp_afdb50 tmp_afdb50_seq