#!/usr/bin/env bash
#
# build.sh [rank] — select one representative reference proteome per taxon at RANK.
#          rank = family (default) | order | genus | class | phylum
#
# Pipeline:
#   1. Input: UniProt reference-proteomes TSV(.gz) downloaded from the UniProt
#      Proteomes website ("Download all", TSV, columns:
#        Proteome Id, Organism, Organism Id, Protein count, BUSCO, CPD, Taxonomic lineage).
#   2. taxonkit reformat -f "{d}\t{RANK}"  ->  map each Organism Id (taxid) to its
#      DOMAIN and its taxon at RANK (both rank-accurate, from the same lookup).
#   3. select_representatives.py --rank RANK ->  group proteomes by the taxon at RANK,
#      keep the most complete one (BUSCO %C primary, CPD-not-Outlier, then protein
#      count), and emit the download path for its gene2acc file.
#
# Both domain and rank value come from taxonkit. taxonkit writes the domain
# capitalized (Bacteria/Archaea/Eukaryota) = the FTP folder the proteome lives under.
# Viruses return an empty domain (dropped); taxa with no value at RANK are dropped.
#
# Requires: taxonkit (v0.20+), python3, an NCBI taxdump at ../taxdump.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PROJ="$(cd "$HERE/.." && pwd)"
TD="$PROJ/taxdump"

RANK="${1:-family}"
SRC="$HERE/proteomes_proteome_type_REFERENCE_2026_07_22.tsv.gz"
[[ -f "$SRC" ]] || { echo "ERROR: reference-proteomes TSV not found: $SRC" >&2; exit 1; }
[[ -f "$TD/nodes.dmp" ]] || { echo "ERROR: taxdump not found at $TD" >&2; exit 1; }

# rank -> taxonkit reformat placeholder
case "$RANK" in
  phylum) PH="{p}";; class) PH="{c}";; order) PH="{o}";;
  family) PH="{f}";; genus) PH="{g}";;
  *) echo "ERROR: unsupported rank '$RANK' (use order|family|genus|class|phylum)" >&2; exit 1;;
esac

# reader for the source TSV (.gz or plain)
read_src() { case "$SRC" in *.gz) gzcat "$SRC";; *) cat "$SRC";; esac; }

echo "[1/2] taxid -> domain + $RANK via taxonkit"
read_src | tail -n +2 | cut -f3 | sort -un \
  | taxonkit reformat --data-dir "$TD" -I 1 -f "{d}\t$PH" 2>/dev/null \
  > "$HERE/taxid2taxonomy.$RANK.tsv"
echo "      $(wc -l < "$HERE/taxid2taxonomy.$RANK.tsv") taxids mapped"

echo "[2/2] select most-complete proteome per $RANK"
python3 "$HERE/select_representatives.py" \
  --tsv "$SRC" \
  --taxid2taxonomy "$HERE/taxid2taxonomy.$RANK.tsv" \
  --rank "$RANK"

echo
echo "Done ($RANK). Outputs in $HERE :"
echo "  selected_proteomes.$RANK.tsv   ($RANK -> chosen proteome + gene2acc URL)"
echo "  coverage_report.$RANK.txt"
echo
echo "Next, to turn the selection into a Foldseek target DB:"
echo "  bash build_subdb.sh selected_proteomes.$RANK.tsv <outdir>"
