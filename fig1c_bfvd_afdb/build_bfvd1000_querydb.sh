#!/bin/bash
# ---------------------------------------------------------------------------
# Run ONCE. Samples N BFVD structures and builds a persistent Foldseek
# query DB from them, saved alongside the fig1c results.
#
# Usage:  bash build_bfvd1000_querydb.sh
# ---------------------------------------------------------------------------

set -euo pipefail

PDBDIR=/fast/sunny/bfvd2/bfvd2_pdb_only
OUT=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c_bfvd1000_afdb50
QDB=$OUT/bfvd1000DB

N=1000
THREADS=16

mkdir -p "$OUT"
cd "$OUT"

if [ -e "${QDB}.dbtype" ]; then
  echo "Query DB already exists: $QDB"
  echo "Delete bfvd1000DB* first if you really want to rebuild it."
  exit 0
fi

# --- inventory and sample ---------------------------------------------------

if [ ! -e sample.txt ]; then
  find "$PDBDIR" \( -name '*.pdb' -o -name '*.pdb.gz' \) | sort > all_pdbs.txt
  echo "found $(wc -l < all_pdbs.txt) structures in $PDBDIR"
  shuf -n "$N" all_pdbs.txt | sort > sample.txt
  echo "sampled $(wc -l < sample.txt) -> sample.txt"
else
  echo "reusing existing sample.txt ($(wc -l < sample.txt) paths)"
fi

# --- build the DB -----------------------------------------------------------
# Paths passed as an array rather than symlinked into a flat directory,
# so nested subdirs with colliding basenames can't silently collide.

mapfile -t FILES < sample.txt
foldseek createdb "${FILES[@]}" "$QDB" --threads "$THREADS"

# --- verify -----------------------------------------------------------------

echo
echo "DB entries: $(wc -l < "${QDB}.lookup")  (expected ~$N)"
echo "components:"
ls -la "${QDB}"* | awk '{print "  " $5 "\t" $9}'
echo
echo "Done. Reusable query DB: $QDB"
