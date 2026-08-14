#!/bin/bash
# Merge the 16 array chunks into one .m8 and label every hit with its proteome
# and family. Run on the login node after all tasks finish -- it is I/O, not compute.
#
# Usage: bash merge_family_gpu.sh
set -euo pipefail

OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd_afdb_family_subset
D=/fast/sunny/databases/afdb_v6_family_subset
SEL=${SEL:-$HOME/virus-host-mimicry/build_reference_db/selected_proteomes.family.tsv}
NCHUNKS=${NCHUNKS:-16}

cd "$OUT"
[ -f "$SEL" ] || { echo "ERROR: selection not found: $SEL (set SEL=...)" >&2; exit 1; }

# --- 1. completeness check before merging ------------------------------------
MISSING=0
for i in $(seq 0 $((NCHUNKS-1))); do
  f=$(printf 'hits/family_gpu_%02d.m8' "$i")
  [ -s "$f" ] || { echo "MISSING or empty: $f"; MISSING=$((MISSING+1)); }
done
[ "$MISSING" -eq 0 ] || {
  echo "ERROR: $MISSING/$NCHUNKS chunks absent -- resubmit those array tasks before merging." >&2
  exit 1; }

cat hits/family_gpu_*.m8 > hits/family_all.m8
echo "merged rows:    $(wc -l < hits/family_all.m8)"
echo "queries hit:    $(cut -f1 hits/family_all.m8 | sort -u | wc -l)"

# --- 2. accession -> proteome label -> family --------------------------------
# Strip the AF- wrapper, the -F<n> fragment tag AND any -<n> isoform suffix:
# AFDB carries isoform entries (AF-A0A0B4KGY6-10-F1-model_v6) which build_subdb.sh
# matched to their canonical accession, so the isoform tail must go too or those
# targets look unmappable.
cut -f2 hits/family_all.m8 \
  | sed -E 's/^AF-//; s/-F[0-9]+.*$//; s/-[0-9]+$//' | sort -u > hit_accs.txt
awk -F'\t' -v OFS='\t' 'NR==FNR{w[$1]=1;next} ($1 in w){print $1,$2}' \
  hit_accs.txt "$D/_meta/acc2label.tsv" | sort -u > hit_acc2label.tsv
awk -F'\t' -v OFS='\t' 'NR>1{print $3"_"$4, $1}' "$SEL" > label2family.tsv

echo "hit accessions: $(wc -l < hit_accs.txt)   mapped: $(wc -l < hit_acc2label.tsv)"

# appends accession(11), label(12), family(13)
awk -F'\t' -v OFS='\t' '
  FILENAME==ARGV[1] { lab[$1]=$2; next }
  FILENAME==ARGV[2] { fam[$1]=$2; next }
  { acc=$2; sub(/^AF-/,"",acc); sub(/-F[0-9]+.*$/,"",acc); sub(/-[0-9]+$/,"",acc)
    l=(acc in lab?lab[acc]:"NA")
    print $0, acc, l, (l in fam?fam[l]:"NA") }
' hit_acc2label.tsv label2family.tsv hits/family_all.m8 > hits/family_all.labeled.m8

NA=$(awk -F'\t' '$13=="NA"' hits/family_all.labeled.m8 | wc -l)
echo "unmapped rows:  $NA  (expect <0.2%)"

# --- 3. per-family counts, raw and normalised by proteome size ---------------
# Raw counts track proteome size: the pilot's top families were a 133,953-protein
# rust fungus and a BUSCO-77.5% swallow, whose hits piled onto ~15% of their
# proteins. Neither ranking is an enrichment -- feed these to compute_background.py
# and enrichment_fig1c.py, do not read them directly.
cut -f13 hits/family_all.labeled.m8 | sort | uniq -c \
  | sed -E 's/^ *([0-9]+) (.*)$/\2\t\1/' | sort > fam_hits.tsv
awk -F'\t' -v OFS='\t' '
  NR==FNR { h[$1]=$2; next }
  FNR>1 && ($1 in h) { print $1, h[$1], $8, h[$1]*1000/$8 }
' fam_hits.tsv "$SEL" | sort -t$'\t' -k4,4gr > fam_hits.normalised.tsv

echo "families hit:   $(wc -l < fam_hits.tsv) / $(($(wc -l < "$SEL")-1))"
echo
echo "top 10 by hits per 1000 proteins (family, hits, proteins, ratio):"
head -10 fam_hits.normalised.tsv | sed 's/^/  /'
echo
echo "wrote hits/family_all.m8, hits/family_all.labeled.m8, fam_hits.normalised.tsv"
