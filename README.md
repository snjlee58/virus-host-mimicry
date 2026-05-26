# virus-host-mimicry

Structure-based detection of viral protein mimicry of host proteins,
following Lasso et al. (Cell Systems 2021).

## Pipeline

1. `run_localcolabfold.sh` — predict structure of a viral protein
2. `run_foldseek.sh` — search predicted structure against PDB100
3. `classify_hits.py` — bucket each hit's host by NCBI division

## Reference data (not in repo)

Reference databases are kept on `/fast` and gitignored. To set up:

### NCBI taxdump
- Source: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
- Location on cluster: `/fast/sunny/virus-host-mimicry/taxdump/`
- Version pinned: **2026-05-26** (MD5: `784fe01e3dd0e231e46ea64c12aa0fd1`)
- Download / refresh:
  ```bash
  curl -O https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  mkdir -p /fast/sunny/virus-host-mimicry/taxdump
  tar -xzf taxdump.tar.gz -C /fast/sunny/virus-host-mimicry/taxdump/
  ```

### Virus-Host DB
Virus → host taxid mapping (Mihara et al. 2016). Used to assign each virus
to its host's taxonomic division.

- Source: https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv
- Location on cluster: `/fast/sunny/virus-host-mimicry/virushostdb/`
- Version pinned: snapshot dated **2026-02-14** (MD5: `66a10783690148dce327d4b90b4338f3`)
- Download / refresh:
  ```bash
  mkdir -p /fast/sunny/virus-host-mimicry/virushostdb
  cd /fast/sunny/virus-host-mimicry/virushostdb
  curl -O https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv
  curl -O https://www.genome.jp/ftp/db/virushostdb/README
  ```
- Verify: `md5sum virushostdb.tsv` (linux) or `md5 virushostdb.tsv` (macOS).
- Note: the source FTP does not publish per-file MD5s; the hash above is
  computed locally from the downloaded snapshot.
