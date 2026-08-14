#!/usr/bin/env bash
# build_subdb.sh -- selected proteomes -> gene2acc -> accessions -> ONE Foldseek target DB
#
#   bash build_subdb.sh selected_proteomes.family.tsv /fast/sunny/databases/afdb_v6_family_subset
#
# Reads the proteome selection, fetches each proteome's gene2acc (cached), collects
# the UniProt accessions, and slices those entries out of an existing AlphaFold
# Foldseek DB into a single pooled target DB under <outdir>.
#
# args:  1 = selection TSV (needs a header with a gene2acc_url column)
#        2 = output directory
#        3 = DB name (default: the selection's rank column, e.g. "family")
#
# env:   SRC           source Foldseek DB (default: afdb50_seq, see below)
#        SUBDB_MODE    1 = symlink parent data (default), 0 = hard copy
#        PER_PROTEOME  1 = also build one DB per proteome under proteomes/<label>/
#        WORKERS       parallel gene2acc downloads (default 16)
#        MIRROR        ebi (default) | uniprot
#        FOLDSEEK_BIN  foldseek binary
#
# WHICH SOURCE DB -- use the MEMBERS side, not the representatives. afdb50 is a
# 50%-identity clustering: its ~66.7M reps are cluster centroids only, so most
# proteins of our family representatives are not reps and would simply be absent.
# afdb50_seq (~241M entries) is effectively all of AFDB. See
# fig1c_bfvd_afdb/build_afdb50_gpu_db.sh for how reps/_seq/_clu relate here.
#
# Prefer an UNPADDED source: subsetting a GPU-padded DB does not regenerate its
# .gpu_mapping sidecars, so the result is not reliably GPU-ready. Slice the plain
# DB, then pad the (much smaller) subset:  foldseek makepaddedseqdb OUT OUT_pad
#
# PARTIAL COVERAGE IS NORMAL -- not every UniProt accession has an AFDB model.
# Unmatched accessions land in missing.txt; feed that to fetch_structures.py.
set -euo pipefail

MANIFEST=${1:?usage: bash build_subdb.sh <selected_proteomes.tsv> <outdir> [name]}
OUTDIR=${2:?usage: bash build_subdb.sh <selected_proteomes.tsv> <outdir> [name]}

SRC=${SRC:-/fast/databases/foldseek/afdb_v6/afdb50_seq}
SUBDB_MODE=${SUBDB_MODE:-1}
PER_PROTEOME=${PER_PROTEOME:-0}
WORKERS=${WORKERS:-16}
MIRROR=${MIRROR:-ebi}
FOLDSEEK=${FOLDSEEK_BIN:-foldseek}

[ -f "$MANIFEST" ] || { echo "ERROR: selection not found: $MANIFEST" >&2; exit 1; }
command -v "$FOLDSEEK" >/dev/null || { echo "ERROR: foldseek not on PATH (set FOLDSEEK_BIN)" >&2; exit 1; }
command -v curl >/dev/null || { echo "ERROR: curl not on PATH" >&2; exit 1; }
[ -f "${SRC}.dbtype" ] || { echo "ERROR: not a Foldseek DB: $SRC (no ${SRC}.dbtype)" >&2; exit 1; }

GENE2ACC=$OUTDIR/gene2acc
PROTEOMES=$OUTDIR/proteomes
META=$OUTDIR/_meta
mkdir -p "$GENE2ACC" "$PROTEOMES" "$META"
LOG=$OUTDIR/build.log; : > "$LOG"

# rank = first column header; also the default DB name
RANK=$(head -1 "$MANIFEST" | cut -f1)
NAME=${3:-$RANK}
OUT=$OUTDIR/$NAME

echo "selection : $MANIFEST"
echo "source DB : $SRC"
echo "target DB : $OUT"
echo "mirror    : $MIRROR   workers: $WORKERS   subdb-mode: $SUBDB_MODE"

# Padded sources are searchable but cannot be sliced into a GPU-ready DB.
if [ -f "${SRC}_ss.dbtype" ] && od -An -tu1 "${SRC}_ss.dbtype" 2>/dev/null | grep -qE '0 +0 +8 +0'; then
  echo "WARNING: ${SRC}_ss is GPU-padded; .gpu_mapping sidecars will NOT be regenerated."
  echo "         Prefer the unpadded DB and pad the subset instead."
fi
[ -f "${SRC}_ca.dbtype" ] || echo "NOTE: ${SRC}_ca absent -> --alignment-type 1 (TMalign) unavailable; type 2 still works."

# ---- 1. selection -> gene2acc URLs (column resolved by header name) ----------
awk -F'\t' -v mirror="$MIRROR" '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i
          if (!("gene2acc_url" in h)) { print "NO_COL" > "/dev/stderr"; exit 1 }
          next }
  $(h["gene2acc_url"]) != "" {
    u = $(h["gene2acc_url"])
    if (mirror == "ebi") sub(/^https:\/\/ftp\.uniprot\.org\/pub\/databases\/uniprot\//,
                             "https://ftp.ebi.ac.uk/pub/databases/uniprot/", u)
    print u
  }
' "$MANIFEST" 2>"$META/awk.err" | sort -u > "$META/urls.txt" || {
  grep -q NO_COL "$META/awk.err" && echo "ERROR: $MANIFEST has no gene2acc_url column" >&2
  exit 1
}
N_PROT=$(grep -c . "$META/urls.txt" || true)
echo
echo "[1/5] proteomes in selection: $N_PROT"
[ "$N_PROT" -gt 0 ] || { echo "ERROR: no gene2acc URLs found" >&2; exit 1; }

# ---- 2. fetch gene2acc, cached and parallel ----------------------------------
# The EBI mirror serves byte-identical files far faster than ftp.uniprot.org
# (~1 s vs ~53 s per file measured), which over ~1900 proteomes is the difference
# between half an hour and most of a day.
: > "$META/failed_urls.txt"
export GENE2ACC META
fetch_one() {
  url=$1
  gz="$GENE2ACC/$(basename "$url")"
  [ -s "$gz" ] && exit 0
  if curl -sSfL --retry 3 --retry-delay 3 --max-time 300 -o "$gz.part" "$url"; then
    mv -f "$gz.part" "$gz"
  else
    rm -f "$gz.part"
    printf '%s\n' "$url" >> "$META/failed_urls.txt"
  fi
}
export -f fetch_one
echo "[2/5] fetching gene2acc (cached; re-runs are free)"
xargs -P "$WORKERS" -I{} bash -c 'fetch_one "$@"' _ {} < "$META/urls.txt" || true
N_FAIL=$(grep -c . "$META/failed_urls.txt" || true)
N_HAVE=$(find "$GENE2ACC" -name '*.gene2acc.gz' -size +0 | wc -l | tr -d ' ')
echo "      cached: $N_HAVE / $N_PROT   failed: $N_FAIL"
[ "$N_FAIL" -eq 0 ] || echo "      (failed URLs in $META/failed_urls.txt -- re-run to retry)"
[ "$N_HAVE" -gt 0 ] || { echo "ERROR: no gene2acc files fetched" >&2; exit 1; }

# ---- 3. gene2acc -> per-proteome accessions ----------------------------------
# gene2acc column 2 is the UniProtKB accession and REPEATS when two gene symbols
# share a translation (CHBEV_001 and CHBEV_334 both -> A0A916NY62), so the dedup
# is required, not cosmetic. Isoform suffixes (-2, -3) are stripped because AFDB
# models the canonical sequence. sed -E, not \+, so this works on BSD sed too.
: > "$META/labels.txt"
while IFS= read -r url; do
  base=$(basename "$url"); label=${base%.gene2acc.gz}
  gz=$GENE2ACC/$base
  [ -s "$gz" ] || continue
  d=$PROTEOMES/$label; mkdir -p "$d"
  if [ ! -s "$d/accessions.txt" ]; then
    gzip -dc "$gz" | cut -f2 | sed -E 's/-[0-9]+$//' \
      | grep -E '^[A-Z0-9]{6,10}$' | sort -u > "$d/accessions.txt"
  fi
  printf '%s\n' "$label" >> "$META/labels.txt"
done < "$META/urls.txt"
sort -u -o "$META/labels.txt" "$META/labels.txt"

# accession -> label, one flat file (keeps the awk join to two inputs and avoids
# putting ~1900 filenames on a command line)
awk 'BEGIN{OFS="\t"}
  { n=split(FILENAME,p,"/"); print $1, p[n-1] }
' $(sed "s|^|$PROTEOMES/|; s|$|/accessions.txt|" "$META/labels.txt") > "$META/acc2label.tsv"
cut -f1 "$META/acc2label.tsv" | sort -u > "$OUTDIR/accessions.txt"
N_ACC=$(grep -c . "$OUTDIR/accessions.txt" || true)
echo "[3/5] proteomes with accessions: $(grep -c . "$META/labels.txt")   unique accessions: $N_ACC"

# ---- 4. ONE pass over the source entry names for every proteome at once ------
# Names come from <src>.lookup, or from <src>_h via prefixid when there is no
# usable .lookup (padded DBs often have none; target names always resolve via _h).
# Matching is fragment-agnostic: AF-<acc>-F1-, -F2-, -F3- all map back to <acc>,
# so large multi-fragment proteins are not silently dropped.
#
# The cache is keyed to the label list AND the accession count, so adding or
# removing proteomes correctly forces a re-join instead of silently reporting the
# new ones as "0 hits in AFDB".
STAMP=$( { cat "$META/labels.txt"; echo "$N_ACC"; } | cksum | awk '{print $1"-"$2}' )
if [ ! -s "$META/key_map.tsv" ] || [ "$(cat "$META/stamp" 2>/dev/null)" != "$STAMP" ]; then
  if [ -f "${SRC}.lookup" ]; then
    echo "[4/5] joining against $(basename "$SRC").lookup (single pass, minutes on afdb50_seq)"
    NAMES="${SRC}.lookup"
  else
    echo "[4/5] no .lookup -- dumping entry names from $(basename "$SRC")_h via prefixid"
    [ -f "${SRC}_h.dbtype" ] || { echo "ERROR: neither ${SRC}.lookup nor ${SRC}_h exists" >&2; exit 1; }
    "$FOLDSEEK" prefixid "${SRC}_h" "$META/key2name.tsv" --tsv >>"$LOG" 2>&1 \
      || { echo "ERROR: prefixid failed, see $LOG" >&2; tail -10 "$LOG" >&2; exit 1; }
    NAMES="$META/key2name.tsv"
  fi
  echo "      sample source name: $(head -1 "$NAMES" | cut -f2)"
  awk -F'\t' '
    NR==FNR { want[$1] = ($1 in want ? want[$1] " " $2 : $2); next }
    {
      n = $2
      if (n ~ /^AF-/) { split(n, a, "-"); acc = a[2] }
      else            { acc = n; sub(/-F[0-9]+.*$/, "", acc) }
      if (acc in want) {
        k = split(want[acc], labs, " ")
        for (i=1; i<=k; i++) print labs[i] "\t" $1 "\t" acc
      }
    }
  ' "$META/acc2label.tsv" "$NAMES" > "$META/key_map.tsv.part"
  mv -f "$META/key_map.tsv.part" "$META/key_map.tsv"
  printf '%s' "$STAMP" > "$META/stamp"
else
  echo "[4/5] reusing cached key_map.tsv (selection unchanged)"
fi
N_MATCH=$(grep -c . "$META/key_map.tsv" || true)
echo "      matched entries (incl. accessions shared by >1 proteome): $N_MATCH"
[ "$N_MATCH" -gt 0 ] || {
  echo "ERROR: nothing matched. Check the accession format against the sample name above." >&2
  exit 1
}

# per-proteome keys/matched/missing. Sorted by label so only one fd is open at a time.
sort -k1,1 -k2,2n "$META/key_map.tsv" \
  | awk -F'\t' -v P="$PROTEOMES" '
      $1 != prev { if (prev != "") { close(kf); close(mf) }
                   prev=$1; kf=P"/"$1"/keys.txt"; mf=P"/"$1"/matched.txt"
                   printf "" > kf; printf "" > mf }
      { print $2 > kf; print $3 > mf }'
while IFS= read -r label; do
  d=$PROTEOMES/$label
  [ -s "$d/matched.txt" ] || { : > "$d/matched.txt"; }
  sort -u "$d/matched.txt" -o "$d/matched.txt"
  comm -23 "$d/accessions.txt" "$d/matched.txt" > "$d/missing.txt"
done < "$META/labels.txt"

# ---- 5. pooled target DB ----------------------------------------------------
cut -f2 "$META/key_map.tsv" | sort -un > "$OUTDIR/keys.txt"
N_KEY=$(grep -c . "$OUTDIR/keys.txt")
cut -f3 "$META/key_map.tsv" | sort -u > "$META/found_acc.txt"
comm -23 "$OUTDIR/accessions.txt" "$META/found_acc.txt" > "$OUTDIR/missing.txt"
N_MISS=$(grep -c . "$OUTDIR/missing.txt" || true)

echo "[5/5] createsubdb -> $OUT   ($N_KEY unique entries)"
for suf in "" _ss _ca _h; do
  [ -f "${SRC}${suf}.index" ] || [ -f "${SRC}${suf}" ] || continue
  # Never redirect into ${OUT}.lookup: with --subdb-mode 1 foldseek symlinks it to
  # ${SRC}.lookup, and a redirect would truncate the shared source DB's lookup.
  "$FOLDSEEK" createsubdb "$OUTDIR/keys.txt" "${SRC}${suf}" "${OUT}${suf}" \
    --subdb-mode "$SUBDB_MODE" >>"$LOG" 2>&1
  [ -s "${OUT}${suf}.index" ] || { echo "ERROR: ${suf:-base} subset is empty, see $LOG" >&2; exit 1; }
  echo "      ${suf:-base}: $(grep -c . "${OUT}${suf}.index") entries"
done
for extra in _mapping _taxonomy; do
  if [ -f "${SRC}${extra}" ] && [ ! -L "${OUT}${extra}" ]; then
    rm -f "${OUT}${extra}"; cp "${SRC}${extra}" "${OUT}${extra}"
    echo "      carried over ${extra}"
  fi
done

# optional: one DB per proteome, for per-host searches. Nearly free under
# --subdb-mode 1 (symlinked data + a small .index each).
if [ "$PER_PROTEOME" = "1" ]; then
  echo "      building per-proteome DBs"
  while IFS= read -r label; do
    d=$PROTEOMES/$label; o=$d/$label
    [ -s "$d/keys.txt" ] || { echo "        $label: 0 hits, skipped"; continue; }
    if [ -s "$o.index" ] && [ "$(grep -c . "$o.index")" -eq "$(grep -c . "$d/keys.txt")" ]; then
      continue
    fi
    for suf in "" _ss _ca _h; do
      [ -f "${SRC}${suf}.index" ] || [ -f "${SRC}${suf}" ] || continue
      "$FOLDSEEK" createsubdb "$d/keys.txt" "${SRC}${suf}" "${o}${suf}" \
        --subdb-mode "$SUBDB_MODE" >>"$LOG" 2>&1
    done
  done < "$META/labels.txt"
fi

# ---- report -----------------------------------------------------------------
{
  echo "target DB:          $OUT"
  echo "source DB:          $SRC"
  echo "selection:          $MANIFEST   ($N_PROT proteomes, $N_HAVE fetched)"
  echo "unique accessions:  $N_ACC"
  echo "structures in DB:   $N_KEY"
  echo "no AFDB model:      $N_MISS   (-> missing.txt)"
  [ "$N_FAIL" -eq 0 ] || echo "gene2acc failures:  $N_FAIL   (-> _meta/failed_urls.txt)"
  awk -v f="$N_KEY" -v n="$N_ACC" 'BEGIN{printf "coverage:           %.2f%%\n", f*100.0/n}'
} > "$OUTDIR/coverage.txt"
echo
cat "$OUTDIR/coverage.txt"
echo
echo "target DB ready: $OUT"
echo "  foldseek easy-search <queries> $OUT hits.m8 tmp --alignment-type 2"
echo "  GPU: foldseek makepaddedseqdb $OUT ${OUT}_pad"
