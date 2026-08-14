#!/usr/bin/env python3
"""
fetch_structures.py
-------------------
Download AFDB structures for the selected reference proteomes into per-proteome
folders, then build ONE Foldseek target DB over the whole tree.

Prefer build_subdb.sh when the structures already exist inside a full AFDB
Foldseek DB (afdb50_seq) - slicing that takes minutes, whereas this script issues
one HTTP request per structure. At the full family selection (~28.7M accessions)
that is roughly two months of wall clock and ~1.4 TB. Use this script for:
  * pilots (--limit) to check the pipeline end to end, and
  * filling gaps, i.e. the accessions in build_subdb.sh's missing.txt.

Input is auto-detected:
  * accessions.<rank>.tsv from fetch_accessions.py (header starting "accession")
    -> structures are filed per proteome, as UP<id>_<taxid>/
  * a plain accession list (one per line), e.g. missing.txt
    -> structures go in a single unassigned/ folder, since a bare list carries no
       proteome provenance

Layout under --outdir:
  structures/UP000000346_666510/AF-*.cif    one folder per proteome
  db, db_ss, db_ca, db_h                    ONE Foldseek target DB
  missing.txt                               accessions with no AFDB model
  downloaded.txt                            accessions present locally

foldseek createdb recurses into subdirectories, so the per-proteome folders are
purely organisational - they still yield a single search target, not one DB per
proteome (which would need one search per proteome).

Per accession:
    primary   https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_{ver}.cif
    fallback  resolve the real URL via https://alphafold.ebi.ac.uk/api/prediction/{acc}
              (some expansion-release entries use non-standard filenames)
Accessions with no AFDB model are recorded, not fatal. Existing files are skipped,
so an interrupted run resumes.

Examples:
  ./fetch_structures.py --input accessions.family.tsv --limit 500 --outdir pilot
  ./fetch_structures.py --input /fast/sunny/databases/afdb_v6_family_subset/family/missing.txt \\
                        --outdir gapfill
"""
import argparse, csv, json, os, subprocess, sys, time, urllib.error, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

AFDB_FILE = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_{ver}.cif"
AFDB_API  = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
HEADERS   = {"User-Agent": "virus-host-mimicry-builder/2.0"}
UNASSIGNED = "unassigned"


def http_get(url, timeout=60, retries=2):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise                       # genuinely absent, do not retry
            last = e
        except Exception as e:              # noqa: BLE001
            last = e
        if attempt < retries:
            time.sleep(2 ** attempt)
    raise last


def read_input(path):
    """-> list of (accession, proteome_folder). Auto-detects provenance TSV vs list."""
    with open(path) as fh:
        first = fh.readline()
    is_prov = first.startswith("accession\t")
    items = []
    with open(path) as fh:
        if is_prov:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                acc = (row.get("accession") or "").strip()
                if not acc:
                    continue
                upid, taxid = (row.get("upid") or "").strip(), (row.get("taxid") or "").strip()
                items.append((acc, f"{upid}_{taxid}" if upid and taxid else UNASSIGNED))
        else:
            for ln in fh:
                acc = ln.strip()
                if acc:
                    items.append((acc, UNASSIGNED))
    return items, is_prov


def download_one(acc, folder, ver, root, timeout):
    """-> (acc, status) with status in ok|ok-api|cached|missing|error."""
    sdir = os.path.join(root, folder)
    dest = os.path.join(sdir, f"AF-{acc}-F1-model_{ver}.cif")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return acc, "cached"
    try:
        data = http_get(AFDB_FILE.format(acc=acc, ver=ver), timeout)
        status = "ok"
    except urllib.error.HTTPError as e:
        if e.code != 404:
            return acc, "error"
        try:
            meta = json.loads(http_get(AFDB_API.format(acc=acc), timeout))
            url = (meta[0].get("cifUrl") or meta[0].get("pdbUrl")) if meta else None
            if not url:
                return acc, "missing"
            data = http_get(url, timeout)
            status = "ok-api"
        except Exception:                   # noqa: BLE001
            return acc, "missing"
    except Exception:                       # noqa: BLE001
        return acc, "error"

    os.makedirs(sdir, exist_ok=True)
    tmp = dest + ".part"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, dest)
    return acc, status


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--input", default=os.path.join(here, "accessions.family.tsv"),
                    help="accessions.<rank>.tsv (per-proteome folders) or a plain accession list")
    ap.add_argument("--outdir", default=os.path.join(here, "afdb_structures"))
    ap.add_argument("--limit", type=int, default=0, help="first N accessions; 0 = all")
    ap.add_argument("--model-version", default="v6")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=60)
    ap.add_argument("--no-foldseek", action="store_true", help="download only")
    ap.add_argument("--foldseek-bin", default="foldseek")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: input not found: {args.input}\n"
                 f"Run fetch_accessions.py first, or point --input at a missing.txt.")
    items, is_prov = read_input(args.input)
    if args.limit > 0:
        items = items[:args.limit]
    if not items:
        sys.exit("ERROR: no accessions in input.")

    # dedupe, first proteome wins (an accession belongs to one reference proteome)
    seen, uniq = set(), []
    for acc, folder in items:
        if acc not in seen:
            seen.add(acc)
            uniq.append((acc, folder))
    items = uniq

    root = os.path.join(args.outdir, "structures")
    os.makedirs(root, exist_ok=True)
    n_folders = len({f for _, f in items})
    print(f"input      : {args.input}")
    print(f"layout     : {'per-proteome folders' if is_prov else 'single unassigned/ folder'}"
          f"  ({n_folders} folder{'s' if n_folders != 1 else ''})")
    print(f"accessions : {len(items)}   workers: {args.workers}")

    if args.limit == 0 and len(items) > 100_000:
        print(f"\nNOTE: {len(items)} accessions at one request each. This is the slow path -\n"
              f"      prefer build_subdb.sh if the structures are already in afdb50_seq.\n")

    stats = {"ok": 0, "ok-api": 0, "cached": 0, "missing": 0, "error": 0}
    missing, present = [], []
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(download_one, acc, folder, args.model_version, root, args.timeout): acc
                for acc, folder in items}
        for f in as_completed(futs):
            acc, status = f.result()
            stats[status] = stats.get(status, 0) + 1
            if status == "missing":
                missing.append(acc)
            elif status != "error":
                present.append(acc)
            done += 1
            if done % 1000 == 0 or done == len(items):
                print(f"  {done}/{len(items)}  ok={stats['ok']} api={stats['ok-api']} "
                      f"cached={stats['cached']} missing={stats['missing']} "
                      f"error={stats['error']}", flush=True)

    with open(os.path.join(args.outdir, "missing.txt"), "w") as fh:
        fh.write("".join(a + "\n" for a in sorted(missing)))
    with open(os.path.join(args.outdir, "downloaded.txt"), "w") as fh:
        fh.write("".join(a + "\n" for a in sorted(present)))

    got = stats["ok"] + stats["ok-api"] + stats["cached"]
    print(f"\nstructures on disk: {got}/{len(items)}")
    print(f"no AFDB model:      {stats['missing']}  -> missing.txt")
    if stats["error"]:
        print(f"transient errors:   {stats['error']}  (re-run to retry; downloads resume)")

    if args.no_foldseek:
        return
    if got == 0:
        sys.exit("ERROR: nothing downloaded, skipping foldseek createdb.")

    # ONE database over the whole tree - createdb recurses into the per-proteome dirs.
    db = os.path.join(args.outdir, "db")
    print(f"\nbuilding Foldseek DB over {root} -> {db}")
    try:
        subprocess.run([args.foldseek_bin, "createdb", root, db],
                       check=True, capture_output=True, text=True)
    except FileNotFoundError:
        sys.exit(f"ERROR: '{args.foldseek_bin}' not found on PATH.")
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR: foldseek createdb failed:\n{e.stderr}")
    n_db = sum(1 for _ in open(db + ".index"))
    print(f"DB entries: {n_db}")
    print("Search it with:")
    print(f"  foldseek easy-search <queries> {db} hits.m8 tmp --alignment-type 1")


if __name__ == "__main__":
    main()
