#!/bin/bash
# Merge the per-chunk Foldseek outputs into one hits.m8, after the array finishes.
# Chunks hold disjoint queries, so a plain concatenation is the full result.
#
# Usage: bash merge_hits.sh [output_dir]

set -euo pipefail

OUT="${1:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb}"

expected=$(find "$OUT/chunks" -name 'list_*' | wc -l | tr -d ' ')
produced=$(find "$OUT/hits" -name 'hits_*.m8' | wc -l | tr -d ' ')
echo "chunk lists: $expected | hits files produced: $produced"
if [[ "$produced" -ne "$expected" ]]; then
  echo "WARNING: $((expected - produced)) chunk(s) have no hits file — check logs/ for failed/timed-out/OOM tasks" >&2
  echo "         re-run a missing chunk with: sbatch --array=<N> run_foldseek_array.sh" >&2
fi

# columns: query target alntmscore qcov tcov evalue taxid taxname taxlineage
cat "$OUT"/hits/hits_*.m8 > "$OUT/hits.m8"
echo "merged -> $OUT/hits.m8 ($(wc -l < "$OUT/hits.m8") rows)"
