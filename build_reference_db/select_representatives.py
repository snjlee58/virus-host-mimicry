#!/usr/bin/env python3
"""
select_representatives.py
-------------------------
Pick ONE representative reference proteome per taxon at a chosen RANK
(order / family / genus / ...), to seed a structure search database.

Inputs:
  --tsv           UniProt reference-proteomes TSV(.gz) downloaded from the UniProt
                  Proteomes site. Columns:
                    Proteome Id | Organism | Organism Id | Protein count | BUSCO | CPD | Taxonomic lineage
  --taxid2taxonomy  Three-column TSV: taxid <TAB> domain <TAB> rankvalue
                    (from `taxonkit reformat -f "{d}\\t{o|f|g}"`)
  --rank          the rank name (family|order|genus|...); sets the first output
                  column header and the output filenames. Default: family.

Both domain and rankvalue come from taxonkit, so classification is rank-accurate
and consistent. taxonkit writes the domain capitalized (Bacteria/Archaea/Eukaryota)
which is exactly the FTP folder the proteome lives under. Viruses come back with an
empty domain, so they drop out naturally; taxa with no value at this rank are dropped.

Selection rule (per taxon at the rank):
  rank candidates by  BUSCO %Complete (desc)   <- primary; CPD is mostly "Unknown"
                 then CPD-not-Outlier (prefer)  <- avoids redundant/incomplete flagged proteomes
                 then Protein count (desc)       <- final tie-break
  keep the top one.

Outputs (written next to this script):
  selected_proteomes.<rank>.tsv  one row per taxon-at-rank with chosen proteome + gene2acc URL
  coverage_report.<rank>.txt     how many proteomes were read/dropped and how many
                                 taxa at this rank ended up with a representative
"""
import argparse, csv, gzip, os, re

CELLULAR = {"Bacteria", "Archaea", "Eukaryota"}
FTP_BASE = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/reference_proteomes")

def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def busco_complete(busco):
    # "C:99.6%[S:97.4%,D:2.2%],F:0.1%,M:0.3%,n:2137" -> 99.6  (missing -> -1.0)
    m = re.match(r"\s*C:([\d.]+)%", busco or "")
    return float(m.group(1)) if m else -1.0

def cpd_is_outlier(cpd):
    return "outlier" in (cpd or "").lower()

def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--taxid2taxonomy", required=True)
    ap.add_argument("--rank", default="family")
    ap.add_argument("--outdir", default=here)
    args = ap.parse_args()

    rank = args.rank

    # taxid -> (domain, rankvalue) from taxonkit
    taxid2tax = {}
    with open(args.taxid2taxonomy) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 3:
                taxid2tax[row[0].strip()] = (row[1].strip(), row[2].strip())

    # walk proteomes, keep best per taxon-at-rank
    best = {}   # rankvalue -> chosen record (dict)
    counts = {} # rankvalue -> number of candidates considered
    n_total = n_no_rank = n_no_domain = 0
    with opener(args.tsv) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)  # header
        for row in r:
            if len(row) < 7:
                continue
            n_total += 1
            upid, organism, taxid, pcount, busco, cpd, _lineage = row[:7]
            domain, value = taxid2tax.get(taxid.strip(), ("", ""))
            if not value:
                n_no_rank += 1
                continue
            if domain not in CELLULAR:   # empty domain == virus / acellular
                n_no_domain += 1
                continue
            try:
                pcount_i = int(pcount)
            except ValueError:
                pcount_i = 0
            rec = {
                "value": value, "domain": domain, "upid": upid,
                "taxid": taxid, "organism": organism,
                "busco_C": busco_complete(busco), "cpd": cpd or "",
                "protein_count": pcount_i,
                "gene2acc_url": f"{FTP_BASE}/{domain}/{upid}/{upid}_{taxid}.gene2acc.gz",
            }
            counts[value] = counts.get(value, 0) + 1
            # sort key: higher is better
            key = (rec["busco_C"], 0 if cpd_is_outlier(cpd) else 1, rec["protein_count"])
            if value not in best or key > best[value]["_key"]:
                rec["_key"] = key
                best[value] = rec

    # write selection
    out_sel = os.path.join(args.outdir, f"selected_proteomes.{rank}.tsv")
    cols = [rank, "domain", "upid", "taxid", "organism",
            "busco_C", "cpd", "protein_count", "n_candidates", "gene2acc_url"]
    with open(out_sel, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for value in sorted(best, key=lambda v: (best[v]["domain"], v)):
            rec = best[value]
            w.writerow([rec["value"], rec["domain"], rec["upid"], rec["taxid"],
                        rec["organism"], f'{rec["busco_C"]:.1f}', rec["cpd"],
                        rec["protein_count"], counts[value], rec["gene2acc_url"]])

    # summary: how many taxa at this rank got a proteome?
    lines = []
    lines.append(f"rank: {rank}")
    lines.append(f"proteomes read (non-header): {n_total}")
    lines.append(f"  dropped, no {rank} rank:   {n_no_rank}")
    lines.append(f"  dropped, non-cellular:     {n_no_domain}")
    lines.append(f"{rank}s with a representative: {len(best)}")
    by_dom = {}
    for v in best:
        by_dom[best[v]["domain"]] = by_dom.get(best[v]["domain"], 0) + 1
    for d in sorted(by_dom):
        lines.append(f"  {d:<10} {by_dom[d]}")

    report = "\n".join(lines)
    with open(os.path.join(args.outdir, f"coverage_report.{rank}.txt"), "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nwrote {out_sel}")

if __name__ == "__main__":
    main()
