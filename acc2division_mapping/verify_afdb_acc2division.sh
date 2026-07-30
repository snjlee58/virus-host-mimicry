#!/usr/bin/env bash
# Post-run checks for afdb_acc2division.py output.
# Usage: bash verify_afdb_acc2division.sh afdb_acc2division.tsv.gz /fast/databases/foldseek/afdb_v6/afdb50_seq /fast/sunny/ncbi_taxdump/taxdump
set -euo pipefail

OUT=${1:-afdb_acc2division.tsv.gz}
DB=${2:-/fast/databases/foldseek/afdb_v6/afdb50_seq}
TAX=${3:-/fast/sunny/ncbi_taxdump/taxdump}

echo "=== 1. row count: output must equal .lookup lines (+1 header) ==="
LK=$(wc -l < "${DB}.lookup")
GZ=$(( $(zcat -f "$OUT" | wc -l) - 1 ))
echo "  .lookup lines : $LK"
echo "  output rows   : $GZ"
[ "$LK" -eq "$GZ" ] && echo "  MATCH" || echo "  *** MISMATCH ***"

echo
echo "=== 2. division breakdown ==="
column -t -s$'\t' afdb_division_counts.tsv 2>/dev/null \
  || zcat -f "$OUT" | awk -F'\t' 'NR>1{c[$3]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' | sort -k2 -rn

echo
echo "=== 3. the 382 unresolved: deleted, or would a newer taxdump fix them? ==="
zcat -f "$OUT" | awk -F'\t' 'NR>1 && $3=="NA" && $2!=0 {print $2}' | sort -u > /tmp/unresolved.txt
echo "  distinct unresolved taxids: $(wc -l < /tmp/unresolved.txt)"

if [ -f "$TAX/delnodes.dmp" ]; then
  sed 's/[^0-9]//g' "$TAX/delnodes.dmp" | sort -u > /tmp/del_ids.txt
  echo "  of which in delnodes.dmp (permanently DELETED by NCBI): \
$(comm -12 /tmp/unresolved.txt /tmp/del_ids.txt | wc -l)"
  echo "  -> deleted taxids will NEVER resolve; a newer taxdump cannot recover them."
else
  echo "  delnodes.dmp not found in $TAX — re-extract the taxdump to run this check."
fi

awk -F'\t\\|\t' '{print $1}' "$TAX/merged.dmp" | sort -u > /tmp/merged_ids.txt
echo "  of which in merged.dmp (should be 0 — those are handled already): \
$(comm -12 /tmp/unresolved.txt /tmp/merged_ids.txt | wc -l)"

echo
echo "=== 4. spot-check reference accessions ==="
echo "  expect: P69905 PRI | P0A7B8 BCT | P00330 PLN | Q58232 BCT | P00476 PHG"
zcat -f "$OUT" | awk -F'\t' '$1=="P69905"||$1=="P0A7B8"||$1=="P00330"||$1=="Q58232"||$1=="P00476"'
