#!/usr/bin/env python3
"""
fetch_accessions.py
-------------------
Stage 3 of the target-DB build: expand the one-proteome-per-family selection into
the flat list of UniProt accessions whose AFDB structures make up the Foldseek
target DB.

Input:
  --selected  selected_proteomes.<rank>.tsv from select_representatives.py.
              Its `gene2acc_url` column points at UniProt's gene2acc.gz for the
              chosen proteome; column 1 is named after the rank (family/order/...).

For each selected proteome: fetch gene2acc.gz and take the UNIQUE column-2
accessions. Column 2 is the UniProtKB accession and it repeats when two or more
genes have identical translations (see UniProt_Reference_Proteomes_README), so
deduping is required, not cosmetic.

Outputs (next to --selected unless --outdir given):
  accessions.txt                 one accession per line, deduped across all proteomes.
                                 (build_subdb.sh does its own fetching; this
                                 script is for the provenance table.)
  accessions.<rank>.tsv          accession -> rank value, domain, upid, taxid.
                                 Provenance: which family/proteome each structure
                                 came from, for grouping hits downstream.
  fetch_accessions.failed.tsv    proteomes whose gene2acc could not be fetched.

Mirror choice: ftp.uniprot.org served a gene2acc.gz in ~53 s where ftp.ebi.ac.uk
served the identical bytes in ~1 s. Over ~1900 proteomes that is the difference
between most of a day and about half an hour, so the EBI mirror is the default.
Pass --mirror uniprot to use the URL exactly as written in the selection TSV.

Fetched files are cached under --cache-dir, so a re-run costs nothing and an
interrupted run resumes.

Examples:
  ./fetch_accessions.py                                  # all families, EBI mirror
  ./fetch_accessions.py --limit 20                       # quick pilot
  ./fetch_accessions.py --domain Bacteria --workers 24
"""
import argparse, csv, gzip, io, os, sys, time, urllib.error, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

HEADERS = {"User-Agent": "virus-host-mimicry-builder/2.0"}
UNIPROT_PREFIX = "https://ftp.uniprot.org/pub/databases/uniprot/"
EBI_PREFIX     = "https://ftp.ebi.ac.uk/pub/databases/uniprot/"


def to_mirror(url, mirror):
    """Rewrite a UniProt FTP URL to the chosen mirror (byte-identical content)."""
    if mirror == "ebi" and url.startswith(UNIPROT_PREFIX):
        return EBI_PREFIX + url[len(UNIPROT_PREFIX):]
    return url


def http_get(url, timeout, retries=3):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:                      # noqa: BLE001 - network is flaky by nature
            last = e
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    raise last


def accessions_for(row, mirror, cache_dir, timeout):
    """Unique column-2 accessions of one proteome's gene2acc, order preserved."""
    tag = f'{row["upid"]}_{row["taxid"]}'
    cached = os.path.join(cache_dir, f"{tag}.gene2acc.gz") if cache_dir else None
    data = None
    if cached and os.path.exists(cached) and os.path.getsize(cached) > 0:
        with open(cached, "rb") as fh:
            data = fh.read()
    if data is None:
        data = http_get(to_mirror(row["gene2acc_url"], mirror), timeout)
        if cached:
            tmp = cached + ".part"
            with open(tmp, "wb") as fh:
                fh.write(data)
            os.replace(tmp, cached)

    seen, accs = set(), []
    with gzip.open(io.BytesIO(data), "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc = parts[1].strip()
            if acc and acc not in seen:
                seen.add(acc)
                accs.append(acc)
    return accs


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--selected", default=os.path.join(here, "selected_proteomes.family.tsv"))
    ap.add_argument("--outdir", default=None, help="default: alongside --selected")
    ap.add_argument("--cache-dir", default=None,
                    help="default: <outdir>/gene2acc_cache; '' disables caching")
    ap.add_argument("--mirror", choices=["ebi", "uniprot"], default="ebi")
    ap.add_argument("--domain", default="", help="filter: Bacteria|Archaea|Eukaryota")
    ap.add_argument("--limit", type=int, default=0, help="first N proteomes; 0 = all")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--timeout", type=int, default=120)
    args = ap.parse_args()

    if not os.path.exists(args.selected):
        sys.exit(f"ERROR: selection not found: {args.selected}\n"
                 f"Run build.sh first.")
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.selected))
    os.makedirs(outdir, exist_ok=True)
    cache_dir = args.cache_dir
    if cache_dir is None:
        cache_dir = os.path.join(outdir, "gene2acc_cache")
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)

    # column 1 is named after the rank the selection was built at
    with open(args.selected) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        rank = r.fieldnames[0]
        rows = [row for row in r if not args.domain or row["domain"] == args.domain]
    if args.limit > 0:
        rows = rows[:args.limit]
    if not rows:
        sys.exit("ERROR: no proteomes selected (check --domain).")

    print(f"rank={rank}  proteomes={len(rows)}  mirror={args.mirror}  "
          f"workers={args.workers}  cache={'off' if not cache_dir else cache_dir}")

    results, failed = {}, []
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(accessions_for, row, args.mirror, cache_dir, args.timeout): row
                for row in rows}
        for f in as_completed(futs):
            row = futs[f]
            done += 1
            try:
                results[row["upid"]] = f.result()
            except Exception as e:                  # noqa: BLE001
                failed.append((row[rank], row["upid"], row["taxid"], repr(e)))
            if done % 100 == 0 or done == len(rows):
                got = sum(len(v) for v in results.values())
                print(f"  {done}/{len(rows)} proteomes  {got} accessions  "
                      f"{len(failed)} failed", flush=True)

    # provenance table + deduped accession list
    prov_path = os.path.join(outdir, f"accessions.{rank}.tsv")
    acc_path  = os.path.join(outdir, "accessions.txt")
    seen = set()
    n_dup = 0
    with open(prov_path, "w", newline="") as fh, open(acc_path, "w") as afh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", rank, "domain", "upid", "taxid"])
        for row in rows:
            for acc in results.get(row["upid"], ()):
                w.writerow([acc, row[rank], row["domain"], row["upid"], row["taxid"]])
                if acc in seen:
                    n_dup += 1
                    continue
                seen.add(acc)
                afh.write(acc + "\n")

    if failed:
        fail_path = os.path.join(outdir, "fetch_accessions.failed.tsv")
        with open(fail_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([rank, "upid", "taxid", "error"])
            w.writerows(failed)

    print(f"\nproteomes fetched : {len(results)}/{len(rows)}")
    print(f"unique accessions : {len(seen)}")
    if n_dup:
        print(f"  ({n_dup} accession rows were duplicates across proteomes, "
              f"kept once in accessions.txt)")
    if failed:
        print(f"FAILED proteomes  : {len(failed)}  -> fetch_accessions.failed.tsv")
        print("  (re-run to retry; successful fetches are cached)")
    print(f"\nwrote {acc_path}\nwrote {prov_path}")
    print("\nNote: build_subdb.sh fetches gene2acc itself; this script is only\n      needed for the combined accession->family provenance table.")


if __name__ == "__main__":
    main()
