#!/usr/bin/env python3
"""Fig 1C: per-division taxonomic enrichment matrix (Foldseek-based).

For each virus group (column = the host division of the virus) and each host
division (row), compute the hypergeometric significance of OVER-representation of
that division among the group's unique structural neighbors, following Lasso
et al. 2021 (verified against scripts/taxonomy_enrichment/3_hypergeometric.pl):

    P(X >= k) = hypergeom.sf(k-1, N, K, n)
      k = unique structural neighbors (PDB targets) in the row division
      n = total unique structural neighbors for the virus group (E-value filtered)
      K = background DB entries in the row division
      N = background total DB entries
    p-values are Bonferroni-corrected across the rows of each column.

Two paper-faithful rules baked in (see project memory 'fig1c-known-divergences'):
  - VERTEBRATE row is INCLUSIVE of human: k and K use (vertebrate + human). The
    paper's vertebrate = human + nonHumanMammal + nonMammalVertebrate; our
    classify_hits buckets are disjoint, so we re-add human here. Other rows are
    disjoint.
  - n and N include ALL buckets (virus/phage/unknown too), consistently on both
    sides; only the 5 host divisions are emitted as rows.

CAVEAT: counts are PDB-chain level, not UniProt/SIFTS-deduped (see
compute_background.py TODO). Chain-level over-weights heavily-crystallized
organisms (human, E. coli); this is a known approximation, not the final figure.

Inputs:
  - background_counts.tsv from `compute_background.py --divisions`
  - <hits-root>/<div>_foldseek/hits.classified.tsv per available column
    (output of classify_hits.py; the bucket is the last column)

Example:
    python3 enrichment_fig1c.py \\
      --hits-root tests/host_divisions \\
      --background tests/host_divisions/background_counts.tsv \\
      --evalue 1e-10 \\
      -o tests/host_divisions/fig1c
"""
import argparse
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
from scipy.stats import hypergeom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Fig 1C rows (hosts) and columns (virus groups by host division), in paper order.
DIVISIONS = ["bacteria", "plantfungi", "invertebrate", "vertebrate", "human"]
DISPLAY = {
    "bacteria": "bacteria",
    "plantfungi": "plant & fungi",
    "invertebrate": "invertebrate",
    "vertebrate": "vertebrate",
    "human": "human",
}

# classified.tsv column indices; bucket is the LAST column (classify_hits appends it).
COL_TARGET = 1
COL_EVALUE = 2


def load_background(path):
    """background_counts.tsv -> (N, {bucket: K})."""
    K = {}
    N = None
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            label, count = parts[0], int(parts[1])
            if label == "_total":
                N = count
            else:
                K[label] = count
    if N is None:
        sys.exit(f"ERROR: no _total row in {path}")
    return N, K


def unique_neighbors(hits_path, evalue_max):
    """Dedup classified hits to unique targets (best/min E-value per target).

    Returns (n, counts): n = unique targets passing the E-value cutoff,
    counts = {bucket: n_targets}.
    """
    best = {}  # target -> (min_evalue, bucket)
    with open(hits_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= COL_EVALUE + 1:
                continue
            target = fields[COL_TARGET]
            try:
                ev = float(fields[COL_EVALUE])
            except ValueError:
                continue
            bucket = fields[-1]
            cur = best.get(target)
            if cur is None or ev < cur[0]:
                best[target] = (ev, bucket)
    counts = defaultdict(int)
    n = 0
    for ev, bucket in best.values():
        if ev <= evalue_max:
            counts[bucket] += 1
            n += 1
    return n, counts


def row_kK(row, counts, K):
    """(k, K) for a Fig 1C row, applying the paper's nested-vertebrate rule."""
    if row == "vertebrate":
        k = counts.get("vertebrate", 0) + counts.get("human", 0)
        Kc = K.get("vertebrate", 0) + K.get("human", 0)
    else:
        k = counts.get(row, 0)
        Kc = K.get(row, 0)
    return k, Kc


def enrich_column(n, counts, N, K):
    """Per-row enrichment for one virus group. Returns {row: stats dict}."""
    rows = {}
    for row in DIVISIONS:
        k, Kc = row_kK(row, counts, K)
        p = 1.0 if (k == 0 or Kc == 0) else float(hypergeom.sf(k - 1, N, Kc, n))
        rows[row] = dict(k=k, K=Kc, n=n, N=N, pval=p,
                         obs=(k / n if n else 0.0),
                         exp=(Kc / N if N else 0.0))
    m = len(DIVISIONS)  # Bonferroni across the rows of this column
    for row in DIVISIONS:
        pc = min(1.0, rows[row]["pval"] * m)
        rows[row]["pcor"] = pc
        rows[row]["logp"] = -np.log10(max(pc, 1e-300))
    return rows


def lasso_cmap_norm():
    """Paper-style colormap: blue (p~1) -> light -> red (very significant)."""
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "lasso", ["#2a5cad", "#e8e8e8", "#a01818"], N=256)
    norm = mcolors.TwoSlopeNorm(vmin=0.0, vcenter=-np.log10(0.05), vmax=20.0)
    return cmap, norm


def plot_matrix(matrix, cols, out_png):
    cmap, norm = lasso_cmap_norm()
    nrow, ncol = len(DIVISIONS), len(cols)
    fig, ax = plt.subplots(figsize=(1.5 * ncol + 2.0, 1.0 * nrow + 1.8))
    for j, col in enumerate(cols):
        for i, row in enumerate(DIVISIONS):
            lp = matrix[row][col]
            color = cmap(norm(lp))
            ax.add_patch(plt.Rectangle((j, nrow - 1 - i), 0.92, 0.92,
                                       facecolor=color, edgecolor="black", lw=0.6))
    ax.set_xlim(-0.1, ncol)
    ax.set_ylim(-0.1, nrow)
    ax.set_xticks([j + 0.46 for j in range(ncol)])
    ax.set_xticklabels([DISPLAY[c] for c in cols], rotation=30, ha="left", fontsize=9)
    ax.xaxis.set_ticks_position("top")
    ax.set_yticks([nrow - 1 - i + 0.46 for i in range(nrow)])
    ax.set_yticklabels([DISPLAY[r] for r in DIVISIONS], fontsize=9)
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("viruses (by host division)", fontsize=11, pad=24)
    ax.set_ylabel("hosts", fontsize=11)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", pad=0.06,
                        shrink=0.8, aspect=24)
    cbar.set_label("pval", fontsize=9)
    cbar.set_ticks([-np.log10(x) for x in (1.0, 0.05, 1e-20)])
    cbar.set_ticklabels(["1", "0.05", r"$1\times10^{-20}$"])
    cbar.ax.invert_xaxis()  # red (significant) on the left, matching the paper
    cbar.ax.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--hits-root", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/tests/host_divisions"),
        help="Dir holding <div>_foldseek/hits.classified.tsv per column",
    )
    ap.add_argument(
        "--background", type=Path,
        default=Path("/fast/sunny/virus-host-mimicry/tests/host_divisions/background_counts.tsv"),
    )
    ap.add_argument("--evalue", type=float, default=1e-10,
                    help="Max E-value for a hit to count as a structural neighbor")
    ap.add_argument("--subdir-template", default="{div}_foldseek/hits.classified.tsv",
                    help="Per-column path under --hits-root")
    ap.add_argument("-o", "--output-prefix", type=Path,
                    default=Path("/fast/sunny/virus-host-mimicry/tests/host_divisions/fig1c"))
    args = ap.parse_args()

    N, K = load_background(args.background)
    print(f"Background: N={N:,}", file=sys.stderr)

    cols, col_paths = [], {}
    for div in DIVISIONS:
        p = args.hits_root / args.subdir_template.format(div=div)
        if p.exists():
            cols.append(div)
            col_paths[div] = p
        else:
            print(f"  (skip column {div}: {p} not found)", file=sys.stderr)
    if not cols:
        sys.exit("ERROR: no <div>_foldseek/hits.classified.tsv found under --hits-root")

    matrix = {row: {} for row in DIVISIONS}
    per_col = {}
    for div in cols:
        print(f"Column {div}: {col_paths[div]}", file=sys.stderr)
        n, counts = unique_neighbors(col_paths[div], args.evalue)
        rows = enrich_column(n, counts, N, K)
        per_col[div] = (n, rows)
        for row in DIVISIONS:
            matrix[row][div] = rows[row]["logp"]
        print(f"  n (unique neighbors, E<={args.evalue:g}): {n:,}", file=sys.stderr)

    out_prefix = args.output_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    out_matrix = Path(str(out_prefix) + "_matrix.tsv")
    with open(out_matrix, "w") as f:
        f.write("host_division\t" + "\t".join(cols) + "\n")
        for row in DIVISIONS:
            f.write(row + "\t" + "\t".join(f"{matrix[row][c]:.4f}" for c in cols) + "\n")

    out_detail = Path(str(out_prefix) + "_detail.tsv")
    with open(out_detail, "w") as f:
        f.write("virus_group\thost_division\tk\tK\tn\tN\tobs\texp\tpval\tpcor\tneglog10_pcor\n")
        for div in cols:
            _, rows = per_col[div]
            for row in DIVISIONS:
                r = rows[row]
                f.write(f"{div}\t{row}\t{r['k']}\t{r['K']}\t{r['n']}\t{r['N']}\t"
                        f"{r['obs']:.4f}\t{r['exp']:.4f}\t{r['pval']:.3e}\t"
                        f"{r['pcor']:.3e}\t{r['logp']:.4f}\n")

    out_png = Path(str(out_prefix) + "_heatmap.png")
    plot_matrix(matrix, cols, out_png)

    # Echo the matrix to the terminal so you see numbers without opening the PNG.
    print("\n-log10(Bonferroni p)  [rows=hosts, cols=virus groups]", file=sys.stderr)
    print("host\\virus".ljust(14) + "".join(c[:11].rjust(13) for c in cols), file=sys.stderr)
    for row in DIVISIONS:
        print(row.ljust(14) + "".join(f"{matrix[row][c]:13.2f}" for c in cols),
              file=sys.stderr)
    print(f"\nWrote:\n  {out_matrix}\n  {out_detail}\n  {out_png}", file=sys.stderr)


if __name__ == "__main__":
    main()
