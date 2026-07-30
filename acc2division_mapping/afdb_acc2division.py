#!/usr/bin/env python3
"""
AFDB accession -> NCBI taxonomic division.

Inputs (all local; the taxdump is the only download):
    <db>.lookup     key <TAB> AF-<acc>-F<n>-model_v<v> <TAB> file#
    <db>_mapping    key <TAB> ncbi_taxid
    taxdump/        nodes.dmp, merged.dmp, division.dmp
                    ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

Output: TSV, one row per AFDB entry (241,070,489 for full AFDB v6)

    uniprot_acc  taxid  div  div_name

    div      NCBI division code, read from division.dmp via field 5 of nodes.dmp
    div_name its long name; omit with --no-name

The 12 NCBI divisions:
    BCT Bacteria                INV Invertebrates   MAM Mammals
    PHG Phages                  PLN Plants and Fungi PRI Primates
    ROD Rodents                 SYN Synthetic and Chimeric
    UNA Unassigned              VRL Viruses         VRT Vertebrates
    ENV Environmental samples
    NA  taxid missing from _mapping, or absent from the taxdump

Quirks of this scheme, all inherent to NCBI (see README):
    - archaea are labelled BCT; there is no archaeal division
    - fungi are inside PLN ("Plants and Fungi")
    - MAM excludes primates and rodents; VRT excludes all mammals
    - protists (e.g. Plasmodium) are in INV
    - ENV means "sampled from the environment" and cuts across clades

Usage:
    python3 afdb_acc2division.py \
        --db /fast/databases/foldseek/afdb_v6/afdb50_seq \
        --taxdump /fast/sunny/ncbi_taxdump/taxdump \
        --out afdb_acc2division.tsv.gz \
        --summary afdb_division_counts.tsv

Requires numpy. ~5 min plain / ~13 min gzipped for 241M rows; ~4 GB RAM.
"""
import argparse, gzip, os, re, sys
import numpy as np


def read_dmp(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            yield line.rstrip("\t|\n").split("\t|\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="DB prefix, e.g. /data/afdb50_seq")
    ap.add_argument("--taxdump", required=True, help="dir containing nodes.dmp etc.")
    ap.add_argument("--out", required=True, help=".tsv or .tsv.gz")
    ap.add_argument("--summary", help="optional per-division count table")
    ap.add_argument("--no-name", action="store_true",
                    help="emit only the 3-letter code, drop div_name")
    a = ap.parse_args()

    # ---- 1. taxid -> division, as an id-indexed array (~3 MB) ----------------
    sys.stderr.write("reading taxdump...\n")
    rows = [(int(f[0]), int(f[4])) for f in read_dmp(os.path.join(a.taxdump, "nodes.dmp"))]
    maxt = max(r[0] for r in rows)
    divnum = np.full(maxt + 1, -1, dtype=np.int8)
    for t, d in rows:
        divnum[t] = d
    del rows

    # merged.dmp: AFDB was built against an older taxonomy snapshot, so some
    # taxids in _mapping have since been retired and merged into another id.
    # Without this they would silently come out NA.
    merged = np.arange(maxt + 1, dtype=np.int32)
    mpath = os.path.join(a.taxdump, "merged.dmp")
    if os.path.exists(mpath):
        for f in read_dmp(mpath):
            old, new = int(f[0]), int(f[1])
            if old <= maxt:
                merged[old] = new

    divcode = {int(f[0]): f[1] for f in read_dmp(os.path.join(a.taxdump, "division.dmp"))}
    divlong = {int(f[0]): f[2] for f in read_dmp(os.path.join(a.taxdump, "division.dmp"))}

    # resolve merges once for every taxid, then look up its division
    allt = np.arange(maxt + 1, dtype=np.int64)
    res = np.where(allt > 0, merged[allt], 0).astype(np.int64)
    res[divnum[res] < 0] = 0                      # merge target must exist
    div = np.full(maxt + 1, -1, dtype=np.int8)
    ok = res > 0
    div[ok] = divnum[res[ok]]

    # string tables indexed by (division + 1), so the row loop does no dict work
    code_tab = np.array(["NA"] + [divcode.get(i, "NA") for i in range(0, 16)], dtype=object)
    name_tab = np.array(["NA"] + [divlong.get(i, "NA") for i in range(0, 16)], dtype=object)

    # ---- 2. key -> taxid as a flat array ------------------------------------
    sys.stderr.write("reading _mapping...\n")
    kparts, tparts = [], []
    with open(a.db + "_mapping", "rb") as fh:
        buf = b""
        while True:
            chunk = fh.read(1 << 26)
            if not chunk:
                break
            buf += chunk
            cut = buf.rfind(b"\n")
            if cut < 0:
                continue
            v = np.array(buf[:cut].split(), dtype=np.int64).reshape(-1, 2)
            kparts.append(v[:, 0].copy())
            tparts.append(v[:, 1].astype(np.int32))
            buf = buf[cut + 1:]
        if buf.strip():
            v = np.array(buf.split(), dtype=np.int64).reshape(-1, 2)
            kparts.append(v[:, 0].copy())
            tparts.append(v[:, 1].astype(np.int32))
    keys = np.concatenate(kparts); del kparts
    taxs = np.concatenate(tparts); del tparts
    kmax = int(keys.max())
    key2tax = np.zeros(kmax + 1, dtype=np.int32)      # 0 == absent from _mapping
    key2tax[keys] = taxs
    del keys, taxs

    # ---- 3. stream .lookup, write output ------------------------------------
    sys.stderr.write("writing rows...\n")
    pat = re.compile(r"AF-([A-Za-z0-9]+)-F\d+")
    with_name = not a.no_name
    counts = {}
    n = no_map = no_tax = 0
    opener = gzip.open if a.out.endswith(".gz") else open

    with open(a.db + ".lookup") as fh, opener(a.out, "wt") as out:
        out.write("uniprot_acc\ttaxid\tdiv\tdiv_name\n" if with_name
                  else "uniprot_acc\ttaxid\tdiv\n")
        for line in fh:
            p = line.split("\t")
            if len(p) < 2:
                continue
            m = pat.search(p[1])
            acc = m.group(1) if m else p[1].split("-")[0]

            key = int(p[0])
            tax = int(key2tax[key]) if key <= kmax else 0
            if tax == 0:
                no_map += 1
                d = -1
            else:
                d = int(div[tax]) if tax <= maxt else -1
                if d < 0:
                    no_tax += 1

            code = code_tab[d + 1]
            if with_name:
                out.write(f"{acc}\t{tax}\t{code}\t{name_tab[d + 1]}\n")
            else:
                out.write(f"{acc}\t{tax}\t{code}\n")

            counts[code] = counts.get(code, 0) + 1
            n += 1
            if n % 20_000_000 == 0:
                sys.stderr.write(f"  {n:,}\n"); sys.stderr.flush()

    sys.stderr.write(f"done: rows={n:,} not_in_mapping={no_map:,} "
                     f"taxid_not_in_taxdump={no_tax:,}\n")

    if a.summary:
        with open(a.summary, "w") as fh:
            fh.write("div\tdiv_name\tn_proteins\tpct\n")
            rev = {v: k for k, v in divcode.items()}
            for c, v in sorted(counts.items(), key=lambda x: -x[1]):
                fh.write(f"{c}\t{divlong.get(rev.get(c, -1), 'NA')}\t{v}\t{100*v/n:.4f}\n")


if __name__ == "__main__":
    main()
