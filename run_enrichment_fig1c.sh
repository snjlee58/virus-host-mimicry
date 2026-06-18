#!/bin/bash
#SBATCH --job-name=enrich1c
#SBATCH --partition=compute          # CPU-only; enrichment needs no GPU
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/enrich1c_%j.log

# Fig 1C taxonomic-enrichment matrix from classified Foldseek hits.
# Mirrors the Fig 1B run, but builds the 5x5 division matrix across every
# <div>_foldseek/hits.classified.tsv found under --hits-root (auto-discovers
# whichever columns you've classified so far).
#
# Submit from your repo dir (so `python enrichment_fig1c.py` resolves):
#   cd ~/virus-host-mimicry
#   sbatch run_enrichment_fig1c.sh                # default E-value cutoff 1e-10
#   sbatch run_enrichment_fig1c.sh 1e-6           # optional: override the cutoff
#
# Outputs (under $ROOT):
#   fig1c_matrix.tsv   -log10(Bonferroni p), rows=hosts x cols=virus groups
#   fig1c_detail.tsv   full k/K/n/N/obs/exp/p per cell
#   fig1c_heatmap.png  the colored matrix

set -euo pipefail

EVALUE="${1:-1e-10}"
ROOT=/fast/sunny/virus-host-mimicry/tests/host_divisions

python enrichment_fig1c.py \
  --hits-root "$ROOT" \
  --background "$ROOT/background_counts.tsv" \
  --evalue "$EVALUE" \
  -o "$ROOT/fig1c"
