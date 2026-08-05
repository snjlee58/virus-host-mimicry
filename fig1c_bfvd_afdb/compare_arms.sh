#!/bin/bash
# Run on the login node after BOTH search jobs finish.
set -euo pipefail

OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50
QDB=$OUT/bfvd1000DB
cd "$OUT"

for f in afdb50.m8 afdb50_seq.m8; do
  [ -e "$f" ] || { echo "missing $f -- has that job finished?" >&2; exit 1; }
done

cut -f1,2 afdb50.m8     | sort -u > afdb50.pairs
cut -f1,2 afdb50_seq.m8 | sort -u > afdb50_seq.pairs

comm -13 afdb50.pairs afdb50_seq.pairs > missed_by_afdb50.pairs
comm -23 afdb50.pairs afdb50_seq.pairs > only_in_afdb50.pairs

{
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
} > summary.txt

cat summary.txt
