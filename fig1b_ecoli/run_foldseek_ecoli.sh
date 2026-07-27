#!/bin/bash
#SBATCH --job-name=foldseek_ecoli
#SBATCH --partition=compute
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --time=12:00:00
# NOTE: SLURM opens these before the job starts, so this dir must already exist
# (mkdir -p it once) and must match $OUTPUT below.
#SBATCH --output=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1b_reproduction_ecoli/foldseek_%j.log
#SBATCH --error=/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1b_reproduction_ecoli/foldseek_%j.err

# Foldseek search for the Fig 1B reproduction (Lasso et al. 2021):
# E. coli-infecting phage structures (from BFVD) vs host proteomes.
#
# Fig 1B compares the structure space of E. coli phages against the E. coli,
# S. cerevisiae, and human proteomes and finds significant structural overlap
# ONLY with E. coli (their natural host). We search the full afdb_proteome DB
# (which contains all three organisms) and restrict to them downstream. Both
# sides are real structures -> TMalign (--alignment-type 1); no ProstT5 / GPU.
#
# Query set = the E. coli-infecting viral proteins, given as UniProt ACCESSIONS
# (one per line) in the accession list. BFVD names each structure by accession:
# <acc>.pdb, or <acc>_1.pdb / <acc>_2.pdb ... for long proteins split into
# fragments. So we link every file matching <acc>.pdb AND <acc>_*.pdb.
#
# Usage (run from INSIDE fig1b_ecoli/ so the default list path resolves):
#   sbatch run_foldseek_ecoli.sh [accession_list] [bfvd_src] [target_db] [output_dir]

set -euo pipefail

ACC_LIST="${1:-ecoli_infecting_accessions.txt}"
BFVD_SRC="${2:-/fast/sunny/bfvd/2023_02_v2/bfvd}"
TARGET_DB="${3:-/fast/databases/foldseek/afdb_v6/afdb_proteome}"
OUTPUT="${4:-/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1b_reproduction_ecoli}"

# --- Setup (uncomment whichever your cluster needs to expose `foldseek`) ---
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate foldseek
# module load foldseek

# --- Sanity checks ---
[[ -f "$ACC_LIST" ]] || { echo "ERROR: accession list not found: $ACC_LIST" >&2; exit 1; }
[[ -d "$BFVD_SRC" ]] || { echo "ERROR: BFVD source dir not found: $BFVD_SRC" >&2; exit 1; }
if [[ ! -f "${TARGET_DB}" && ! -f "${TARGET_DB}.dbtype" ]]; then
  echo "ERROR: Foldseek target DB not found: $TARGET_DB" >&2; exit 1
fi

mkdir -p "$OUTPUT"
QUERY_DIR="$OUTPUT/query_pdbs"
mkdir -p "$QUERY_DIR"

# --- Build the query subset ---
# BFVD files are named by accession (<acc>.pdb), with optional fragment suffixes
# (<acc>_1.pdb, <acc>_2.pdb ...). For each accession, link every matching file.
# The literal "." / "_" right after the accession anchors the match so it can't
# spill onto a longer accession that merely shares this prefix.
: > "$OUTPUT/missing_accessions.txt"
n_acc=0; linked=0; missing=0
while IFS= read -r acc; do
  [[ -z "$acc" ]] && continue
  n_acc=$((n_acc + 1))
  found=0
  for pdb in "$BFVD_SRC/$acc".pdb "$BFVD_SRC/$acc"_*.pdb; do
    [[ -e "$pdb" ]] || continue
    ln -sf "$pdb" "$QUERY_DIR/$(basename "$pdb")"
    linked=$((linked + 1)); found=1
  done
  [[ "$found" -eq 1 ]] || { echo "$acc" >> "$OUTPUT/missing_accessions.txt"; missing=$((missing + 1)); }
done < "$ACC_LIST"
echo "Accessions read: $n_acc | PDB files linked: $linked | accessions with no file: $missing (see missing_accessions.txt)"
[[ "$linked" -gt 0 ]] || { echo "ERROR: no query PDBs were linked" >&2; exit 1; }

TMP="${SCRATCH:-$OUTPUT/tmp}"
mkdir -p "$TMP"

echo "Query dir: $QUERY_DIR ($linked structures)"
echo "Target DB: $TARGET_DB"
echo "Output:    $OUTPUT/hits.m8"

# --- Run Foldseek (flags mirror run_foldseek.sh) ---
# --alignment-type 1     : TMalign structural superposition
# --tmscore-threshold 0.5: only keep TM-score >= 0.5 hits (the "significant
#                          structural relationship" cutoff; SAS < 2.5 Angstrom analog)
# --cov-mode 2 -c 0.7    : require the alignment to cover >= 70% of the query
srun foldseek easy-search \
  "$QUERY_DIR" \
  "$TARGET_DB" \
  "$OUTPUT/hits.m8" \
  "$TMP" \
  --alignment-type 1 --tmscore-threshold 0.5 \
  --cov-mode 2 -c 0.7 \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,qcov,tcov,qlen,tlen,alnlen,fident,lddt,prob,evalue,bits,taxid,taxname,taxlineage" \
  --threads "${SLURM_CPUS_PER_TASK:-64}"

echo "Done. Top hits:"
head -n 10 "$OUTPUT/hits.m8"
