#!/usr/bin/env bash
#
# derive_orders_tree_of_life.sh
# -----------------------------
# Derive the list of all taxonomic ORDERS in the tree of life (cellular life
# only, viruses excluded) from an NCBI taxonomy dump, using taxonkit.
#
# Output: orders_tree_of_life.tsv with columns:
#   1. taxid       NCBI taxonomy id of the order
#   2. domain      Bacteria | Archaea | Eukaryota
#   3. order       scientific name of the order
#   4. lineage     full name lineage (root -> order), ';'-separated
#
# Rationale: at "order" host-range resolution you take one representative
# proteome per order, so the order count is the (upper-bound) proteome budget.
# Note: order count != usable proteomes — not every order has a sequenced,
# well-annotated reference proteome (esp. Candidatus/environmental lineages).
#
# Requires: taxonkit (tested with v0.20.0)
# Usage:    ./derive_orders_tree_of_life.sh [TAXDUMP_DIR]

set -euo pipefail

# Directory containing nodes.dmp, names.dmp, etc. (NCBI taxdump)
# (taxdump lives in the project root, one level up from this build folder)
TD="${1:-$(dirname "$0")/../taxdump}"
OUT="$(dirname "$0")/orders_tree_of_life.tsv"

if [[ ! -f "$TD/nodes.dmp" ]]; then
  echo "ERROR: nodes.dmp not found in '$TD'" >&2
  echo "Download with: taxonkit or 'wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'" >&2
  exit 1
fi

# 1. Pull every tax_id whose rank is exactly 'order' straight from nodes.dmp.
#    nodes.dmp is delimited by '\t|\t'; col1=taxid, col3=rank.
# 2. taxonkit reformat -f "{d}"  -> append the DOMAIN (Bacteria/Archaea/Eukaryota;
#    empty for viruses, which have no NCBI 'domain').
# 3. taxonkit lineage            -> append the full name lineage.
# 4. awk                         -> order name = last ';'-separated lineage element.
# 5. keep only cellular domains  -> drops viruses & unclassified. Sort by domain,order.
awk -F'\t\\|\t' '$3=="order"{print $1}' "$TD/nodes.dmp" \
  | taxonkit reformat --data-dir "$TD" -I 1 -f "{d}" 2>/dev/null \
  | taxonkit lineage  --data-dir "$TD" -i 1        2>/dev/null \
  | awk -F'\t' 'BEGIN{OFS="\t"}{
        n = split($3, a, ";"); ordername = a[n];
        print $1, $2, ordername, $3
      }' \
  | awk -F'\t' '$2=="Bacteria" || $2=="Archaea" || $2=="Eukaryota"' \
  | sort -t$'\t' -k2,2 -k3,3 \
  > "$OUT"

echo "Wrote $(wc -l < "$OUT") non-viral orders -> $OUT"
echo
echo "Orders by domain:"
cut -f2 "$OUT" | sort | uniq -c | sort -rn
