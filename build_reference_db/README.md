# build_reference_db

Build a Foldseek **target database** of host protein structures: one representative
reference proteome per taxonomic family, structures from AlphaFold DB.

Used as the target side of the virus→host structural mimicry search (queries = BFVD
viral structures).

## Workflow

```
proteomes_proteome_type_REFERENCE_2026_07_22.tsv.gz     UniProt Reference Proteomes
        │                                               (manual download, see below)
        │  build.sh  [1/2]  taxonkit reformat "{d}\t{f}"
        ├────────────────────────────────────────────►  taxid2taxonomy.family.tsv
        │
        │  build.sh  [2/2]  select_representatives.py
        │    · drop proteomes with no documented family
        │    · drop non-cellular (viruses)
        │    · group by family, keep best BUSCO %C -> CPD-not-Outlier -> protein count
        ├────────────────────────────────────────────►  selected_proteomes.family.tsv
        └────────────────────────────────────────────►  coverage_report.family.txt
                             │
                             │  fetch_accessions.py   (gene2acc per proteome)
                             ├───────────────────────►  accessions.txt
                             └───────────────────────►  accessions.family.tsv  (provenance)
                             │
                             │  build_subdb.sh        ← preferred
                             │    slice afdb50_seq by accession
                             └───────────────────────►  family{,_ss,_ca,_h}       ONE target DB
                             │
                             │  fetch_structures.py   ← fallback / pilot only
                             └───────────────────────►  structures/<upid>_<taxid>/
                                                        db{,_ss,_ca,_h}           ONE target DB
```

## Running it

```bash
# 1-2. select one proteome per family  (family is the default rank)
./build.sh                      # or: ./build.sh order|genus|class|phylum

# 3. gene2acc -> accessions -> ONE pooled target DB
bash build_subdb.sh selected_proteomes.family.tsv /fast/sunny/databases/afdb_v6_family_subset
# DB lands at <outdir>/family   (third arg overrides the name)
```

`build_subdb.sh` fetches the gene2acc files itself (cached under `<outdir>/gene2acc`,
so re-runs are free), so `fetch_accessions.py` is **optional** — it exists only to
produce the combined `accessions.<rank>.tsv` provenance table.

Useful env overrides:

| var | default | |
|---|---|---|
| `SRC` | `/fast/databases/foldseek/afdb_v6/afdb50_seq` | source Foldseek DB |
| `SUBDB_MODE` | `1` | 1 = symlink parent data, 0 = hard copy |
| `PER_PROTEOME` | `0` | 1 = also build one DB per proteome |
| `WORKERS` | `16` | parallel gene2acc downloads |
| `MIRROR` | `ebi` | `uniprot` to use the URLs as written |
| `DL` | auto | downloader: `curl`, `wget` or `python3`. Auto-detected in that order; only needed when something isn't already cached, so a warm `<outdir>/gene2acc` needs none. Compute nodes here have no `curl`. |

Then search it:

```bash
foldseek easy-search <bfvd_queries> \
  /fast/sunny/databases/afdb_v6_family_subset/family hits.m8 tmp --alignment-type 2
```

Requires `taxonkit` (v0.20+), `python3`, `foldseek`, and an NCBI taxdump at `../taxdump`.

## Why step 4 slices instead of downloading

AlphaFold DB ships bulk per-proteome tars for only ~46 model organisms — **9 of our
1,912 family representatives**. The selection spans ~28.3M proteins, so fetching
structures one HTTP request at a time is roughly two months of wall clock and ~1.4 TB
of `.cif` files. Slicing an on-disk AFDB Foldseek DB takes minutes instead.

`fetch_structures.py` remains for pilots (`--limit`) and for filling the gaps listed in
the subset's `missing.txt`.

### Which source DB

Use the **members** side, `afdb50_seq` (~241M entries ≈ all of AFDB) — *not* the
representatives. `afdb50` is a 50%-identity clustering whose ~66.7M reps are cluster
centroids only, so most proteins of our family representatives are not reps and would
simply be absent. See [../fig1c_bfvd_afdb/build_afdb50_gpu_db.sh](../fig1c_bfvd_afdb/build_afdb50_gpu_db.sh)
for how the reps / `_seq` / `_clu` pieces relate on this cluster.

Prefer an **unpadded** source. Padding (`makepaddedseqdb`) is what makes a DB
GPU-searchable, and subsetting a padded DB does not regenerate its `.gpu_mapping`
sidecars. Slice the plain DB, then pad the much smaller subset:

```bash
foldseek makepaddedseqdb <out_db> <out_db>_pad
```

### Mechanics, verified against foldseek 10.941cd33

- **Fragment-agnostic matching.** AFDB names are `AF-<acc>-F1-model_v6`, and large
  proteins carry `-F2-`, `-F3-` … fragments. Matching splits the accession out of the
  name, so multi-fragment proteins are not silently dropped. Bare-accession names work
  too.
- **`gene2acc` column 2 repeats** when two gene symbols share a translation
  (`CHBEV_001` and `CHBEV_334` → `A0A916NY62`), so the dedup is required, not cosmetic.
  Isoform suffixes (`-2`, `-3`) are stripped first, since AFDB models the canonical
  sequence.
- **Matching is by numeric DB key**, for base / `_ss` / `_ca` / `_h` alike. This avoids
  `createsubdb --id-mode 1`, which needs a `.lookup` beside every DB it touches (`_ss`
  and `_ca` have none) and exits non-zero on a structure DB even after succeeding.
- **Names come from `<src>.lookup`, or from `<src>_h` via `prefixid --tsv`** when the
  source has no usable `.lookup` — as padded DBs often don't. Target names always
  resolve through `_h`.
- **`--subdb-mode 1` symlinks parent data** instead of copying it, which matters with
  `/fast` near capacity. It also symlinks `<out>.lookup` → `<src>.lookup`, so the script
  must never redirect into `<out>.lookup`: that would write through the symlink and
  **truncate the shared source DB's lookup**. It doesn't need to — `createsubdb`
  preserves keys, so the parent lookup resolves subset names correctly.
- Partial coverage is normal; unmatched accessions land in `missing.txt`.

## Input file

`proteomes_proteome_type_REFERENCE_2026_07_22.tsv.gz` — UniProt Proteomes site,
"Download all" as TSV:

| Proteome Id | Organism | Organism Id | Protein count | BUSCO | CPD | Taxonomic lineage |
|---|---|---|---|---|---|---|
| UP000000625 | Escherichia coli (strain K12) | 83333 | 4403 | `C:100.0%[S:99.1%,D:0.9%],F:0.0%,M:0.0%,n:440` | Unknown | cellular organisms, Bacteria, …, Escherichia coli |

The `Taxonomic lineage` column is deliberately **not** used for classification — family
and domain are re-derived with taxonkit so they are rank-accurate rather than positional.

`uniprot_reference_proteomes/` holds the matching UniProt `STATS` and `README`. The
README is the authoritative spec for the `gene2acc` format that `fetch_accessions.py`
parses (column 2 repeats when two genes have identical translations, hence the dedup);
`STATS` gives per-proteome gene2acc counts and assembly accessions.

## Current numbers (family rank)

| | |
|---|---|
| proteomes read | 36,465 |
| dropped, no family rank | 2,668 |
| dropped, non-cellular | 11,979 |
| **families with a representative** | **1,912** (Archaea 57 / Bacteria 759 / Eukaryota 1,096) |
| proteins across the selection | 28,292,926 |
| gene2acc accessions (actual download count) | 28,668,868 |

## Known caveats

- **Viruses are excluded** beyond your stated "family must be documented" rule. Viruses
  *do* have families (`Birnaviridae`), so they are dropped by the non-cellular filter at
  `select_representatives.py:84-86`, not the family filter. That is correct for a *host*
  target DB — but it is an extra filter, so it is called out here.
- **295 cellular taxids are silently lost** for having no family node in NCBI (183
  Bacteria, 63 Eukaryota, 49 Archaea) — mostly `incertae sedis` and `Candidatus`
  lineages. Their proteomes exist in UniProt but can never enter a family-grouped DB.
- **Eukaryota dominates**: 1,096 families / 25.2M proteins vs 2.9M bacterial. Cap with
  `--domain` or `--limit` for a first pass.
- **`gene2acc_url` is synthesized, not verified** — it interpolates
  `{domain}/{upid}/{upid}_{taxid}` against `current_release`, while the proteomes TSV is a
  July 2026 snapshot and local `STATS`/`README` are release 2026_02. Retired proteomes
  404 at fetch time and land in `fetch_accessions.failed.tsv`. Re-download all three
  from the same release when refreshing.
