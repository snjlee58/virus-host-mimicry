#!/usr/bin/env python3
"""
select_representatives.py
-------------------------
Pick ONE representative reference proteome per taxonomic ORDER, to seed an
order-level structure search database.

Inputs:
  --tsv           UniProt reference-proteomes TSV(.gz) downloaded from the UniProt
                  Proteomes site. Columns:
                    Proteome Id | Organism | Organism Id | Protein count | BUSCO | CPD | Taxonomic lineage
  --taxid2taxonomy  Three-column TSV: taxid <TAB> domain <TAB> order
                    (from `taxonkit reformat -f "{d}\\t{o}"`)

Both domain (superregnum) and order come from taxonkit, so classification is
rank-accurate and consistent. taxonkit writes the domain capitalized
(Bacteria/Archaea/Eukaryota) which is exactly the FTP folder the proteome lives
under. Viruses come back with an empty domain, so they drop out naturally;
taxa with no order rank are also dropped.

Selection rule (per order):
  rank candidates by  BUSCO %Complete (desc)   <- primary; CPD is mostly "Unknown"
                 then CPD-not-Outlier (prefer)  <- avoids redundant/incomplete flagged proteomes
                 then Protein count (desc)       <- final tie-break
  keep the top one.

Outputs (written next to this script):
  selected_proteomes.order.tsv  one row per order with the chosen proteome + gene2acc URL
  coverage_report.txt           how many of the tree-of-life orders got a proteome
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
    ap.add_argument("--orders-all",
                    default=os.path.join(here, "orders_tree_of_life.tsv"),
                    help="full order list, for the coverage denominator")
    ap.add_argument("--outdir", default=here)
    args = ap.parse_args()

    # taxid -> (domain, order) from taxonkit
    taxid2tax = {}
    with open(args.taxid2taxonomy) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 3:
                taxid2tax[row[0].strip()] = (row[1].strip(), row[2].strip())

    # walk proteomes, keep best per order
    best = {}   # order -> chosen record (dict)
    counts = {} # order -> number of candidates considered
    n_total = n_no_order = n_no_domain = 0
    with opener(args.tsv) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)  # header
        for row in r:
            if len(row) < 7:
                continue
            n_total += 1
            upid, organism, taxid, pcount, busco, cpd, _lineage = row[:7]
            domain, order = taxid2tax.get(taxid.strip(), ("", ""))
            if not order:
                n_no_order += 1
                continue
            if domain not in CELLULAR:   # empty domain == virus / acellular
                n_no_domain += 1
                continue
            try:
                pcount_i = int(pcount)
            except ValueError:
                pcount_i = 0
            rec = {
                "order": order, "domain": domain, "upid": upid,
                "taxid": taxid, "organism": organism,
                "busco_C": busco_complete(busco), "cpd": cpd or "",
                "protein_count": pcount_i,
                "gene2acc_url": f"{FTP_BASE}/{domain}/{upid}/{upid}_{taxid}.gene2acc.gz",
            }
            counts[order] = counts.get(order, 0) + 1
            # sort key: higher is better
            key = (rec["busco_C"], 0 if cpd_is_outlier(cpd) else 1, rec["protein_count"])
            if order not in best or key > best[order]["_key"]:
                rec["_key"] = key
                best[order] = rec

    # write selection
    out_sel = os.path.join(args.outdir, "selected_proteomes.order.tsv")
    cols = ["order", "domain", "upid", "taxid", "organism",
            "busco_C", "cpd", "protein_count", "n_candidates", "gene2acc_url"]
    with open(out_sel, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for order in sorted(best, key=lambda o: (best[o]["domain"], o)):
            rec = best[order]
            w.writerow([rec["order"], rec["domain"], rec["upid"], rec["taxid"],
                        rec["organism"], f'{rec["busco_C"]:.1f}', rec["cpd"],
                        rec["protein_count"], counts[order], rec["gene2acc_url"]])

    # coverage: how many tree-of-life orders got a proteome?
    all_orders_by_domain = {}
    if os.path.exists(args.orders_all):
        with open(args.orders_all) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 3:
                    all_orders_by_domain.setdefault(row[1], set()).add(row[2])
    covered = set(best)
    lines = []
    lines.append(f"proteomes read (non-header): {n_total}")
    lines.append(f"  dropped, no order rank:    {n_no_order}")
    lines.append(f"  dropped, non-cellular:     {n_no_domain}")
    lines.append(f"orders with a representative: {len(best)}")
    by_dom = {}
    for o in best:
        by_dom[best[o]["domain"]] = by_dom.get(best[o]["domain"], 0) + 1
    for d in sorted(by_dom):
        lines.append(f"  {d:<10} {by_dom[d]}")
    if all_orders_by_domain:
        lines.append("")
        lines.append("coverage vs tree-of-life orders (orders_tree_of_life.tsv):")
        for d in sorted(all_orders_by_domain):
            total = len(all_orders_by_domain[d])
            got = len(covered & all_orders_by_domain[d])
            lines.append(f"  {d:<10} {got}/{total} orders covered")
        allset = set().union(*all_orders_by_domain.values())
        lines.append(f"  {'TOTAL':<10} {len(covered & allset)}/{len(allset)} orders covered")
        missing = sorted(allset - covered)
        lines.append(f"  orders with NO reference proteome: {len(missing)}")

    report = "\n".join(lines)
    with open(os.path.join(args.outdir, "coverage_report.txt"), "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nwrote {out_sel}")

if __name__ == "__main__":
    main()
