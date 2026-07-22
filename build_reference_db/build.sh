#!/usr/bin/env bash
#
# build.sh — select one representative reference proteome per taxonomic ORDER.
#
# Pipeline:
#   1. Input: UniProt reference-proteomes TSV(.gz) downloaded from the UniProt
#      Proteomes website ("Download all", TSV, columns:
#        Proteome Id, Organism, Organism Id, Protein count, BUSCO, CPD, Taxonomic lineage).
#   2. taxonkit reformat -f "{d}\t{o}"  ->  map each Organism Id (taxid) to its
#      DOMAIN and ORDER (both rank-accurate, from the same lookup).
#   3. select_representatives.py   ->  group proteomes by order, keep the most
#      complete one (BUSCO %C primary, CPD-not-Outlier, then protein count),
#      and emit the download path for its gene2acc file.
#
# Both domain and order come from taxonkit. taxonkit writes the domain capitalized
# (Bacteria/Archaea/Eukaryota) which is exactly the FTP folder the proteome lives
# under. Viruses return an empty domain (dropped); no-order taxa are dropped too.
#
# Requires: taxonkit (v0.20+), python3, an NCBI taxdump.
# Usage:    ./build.sh [REFERENCE_PROTEOMES_TSV_GZ]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PROJ="$(cd "$HERE/.." && pwd)"
TD="$PROJ/taxdump"

SRC="${1:-$HERE/proteomes_proteome_type_REFERENCE_2026_07_22.tsv.gz}"
[[ -f "$SRC" ]] || { echo "ERROR: reference-proteomes TSV not found: $SRC" >&2; exit 1; }
[[ -f "$TD/nodes.dmp" ]] || { echo "ERROR: taxdump not found at $TD" >&2; exit 1; }

echo "[1/2] taxid -> domain + order via taxonkit"
# skip header, take Organism Id (col 3), unique, resolve domain + order rank
gzcat "$SRC" | tail -n +2 | cut -f3 | sort -un \
  | taxonkit reformat --data-dir "$TD" -I 1 -f "{d}\t{o}" 2>/dev/null \
  > "$HERE/taxid2taxonomy.tsv"
echo "      $(wc -l < "$HERE/taxid2taxonomy.tsv") taxids mapped"

echo "[2/2] select most-complete proteome per order"
python3 "$HERE/select_representatives.py" \
  --tsv "$SRC" \
  --taxid2taxonomy "$HERE/taxid2taxonomy.tsv" \
  --orders-all "$HERE/orders_tree_of_life.tsv"

echo
echo "Done. Outputs in $HERE :"
echo "  selected_proteomes.order.tsv  (order -> chosen proteome + gene2acc URL)"
echo "  coverage_report.txt"
