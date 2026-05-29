#!/usr/bin/env python3
"""Hypergeometric taxonomic enrichment + Figure 1B plot.

Reproduces Lasso et al. 2021 Cell Systems Figure 1B for one virus group
(e.g. E. coli phages) against three host species (E. coli, S. cerevisiae,
H. sapiens).

Pipeline:
  1. Read hits.classified.tsv (output of classify_hits.py).
  2. Drop rows with evalue > threshold (default: 1e-3).
  3. Deduplicate to unique structural neighbors at the target (PDB chain) level.
  4. Count neighbors per species by matching the taxlineage column.
  5. Run hypergeometric P(X >= k) for each species using N, K from
     compute_background.py.
  6. Save a heatmap-style figure colored by -log10(p), matching the paper.

Hits TSV columns (0-indexed) — matches classify_hits.py output when foldseek
was run with --format-output "query,target,evalue,bits,taxid,taxname,taxlineage":
   0 query, 1 target, 2 evalue, 3 bits, 4 taxid, 5 taxname, 6 taxlineage, 7 bucket

Example:
    python3 enrichment_fig1b.py tests/ecoli/foldseek_output/hits.classified.tsv \\
        --background background_counts.tsv \\
        --evalue 1e-3 \\
        --title "viruses of E. coli" \\
        -o figure_1B.png
"""
import argparse
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# (display label, taxlineage substring, background_counts.tsv key)
SPECIES = [
    ("E. coli",       "s_Escherichia coli",         "ecoli"),
    ("S. cerevisiae", "s_Saccharomyces cerevisiae", "yeast"),
    ("H. sapiens",    "s_Homo sapiens",             "human"),
]

# Column indices in hits.classified.tsv
COL_TARGET = 1
COL_EVALUE = 2
COL_LINEAGE = 6


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("hits", type=Path, help="hits.classified.tsv")
    ap.add_argument(
        "--background", type=Path, default=Path("background_counts.tsv"),
        help="TSV from compute_background.py",
    )
    ap.add_argument(
        "--evalue", type=float, default=1e-3,
        help="Max evalue for a hit to count as a structural neighbor",
    )
    ap.add_argument(
        "--title", default="viruses of E. coli",
        help="Plot title (appears above the cells)",
    )
    ap.add_argument(
        "-o", "--output", type=Path, default=Path("figure_1B.png"),
    )
    return ap.parse_args()


def load_background(path):
    counts = {}
    N = None
    with open(path) as f:
        next(f)  # header
        for line in f:
            sp, c = line.rstrip("\n").split("\t")
            c = int(c)
            if sp == "_total":
                N = c
            else:
                counts[sp] = c
    if N is None:
        sys.exit(f"ERROR: no _total row in {path}")
    return N, counts


def filter_and_dedupe(hits_path, evalue_max):
    """Return {target -> lineage} for hits passing evalue cutoff."""
    unique = {}
    n_rows = 0
    n_kept = 0
    with open(hits_path) as f:
        for line in f:
            n_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= COL_LINEAGE:
                continue
            try:
                ev = float(fields[COL_EVALUE])
            except ValueError:
                continue
            if ev > evalue_max:
                continue
            n_kept += 1
            target = fields[COL_TARGET]
            if target not in unique:
                unique[target] = fields[COL_LINEAGE]
    return unique, n_rows, n_kept


def count_per_species(unique):
    k = defaultdict(int)
    for lineage in unique.values():
        for label, tag, _ in SPECIES:
            if tag in lineage:
                k[label] += 1
                break
    return k


def hypergeom_pvals(k_counts, n, K_by_sp, N):
    """For each species: P(X >= k) under hypergeometric(N, K, n)."""
    pvals = {}
    for label, _, bg_key in SPECIES:
        k = k_counts[label]
        K = K_by_sp.get(bg_key, 0)
        if k == 0 or K == 0:
            pvals[label] = 1.0
        else:
            pvals[label] = float(hypergeom.sf(k - 1, N, K, n))
    return pvals


def plot(pvals, k_counts, n, title, output_path):
    """Stack of cells colored by -log10(p), mimicking Figure 1B."""
    labels = [l for l, _, _ in SPECIES]
    ps = [pvals[l] for l in labels]
    logp = np.array([-np.log10(max(p, 1e-20)) for p in ps])

    # Paper-style colormap: red (very significant) - light - blue (not)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "lasso", ["#2a5cad", "#e8e8e8", "#a01818"], N=256
    )
    # Center the colormap at p = 0.05 (-log10 ≈ 1.3), saturate at p = 1e-20
    norm = mcolors.TwoSlopeNorm(vmin=0.0, vcenter=-np.log10(0.05), vmax=20.0)

    fig, ax = plt.subplots(figsize=(3.2, 4.2))
    cell_w, cell_h, gap = 0.9, 0.9, 0.08

    for i, (label, p, lp) in enumerate(zip(labels, ps, logp)):
        y = (len(labels) - 1 - i) * (cell_h + gap)
        color = cmap(norm(lp))
        ax.add_patch(plt.Rectangle((0.6, y), cell_w, cell_h,
                                    facecolor=color, edgecolor="black", linewidth=0.8))
        # Species label on the left
        ax.text(0.5, y + cell_h / 2, label, va="center", ha="right",
                fontsize=11, fontstyle="italic")
        # p-value annotation inside the cell
        txt_color = "white" if lp > 5 else "black"
        ax.text(0.6 + cell_w / 2, y + cell_h / 2,
                f"p={p:.1e}\nk={k_counts[label]}/n={n}",
                va="center", ha="center", fontsize=7.5, color=txt_color)

    ax.set_xlim(-0.5, 2.0)
    ax.set_ylim(-0.3, len(labels) * (cell_h + gap))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=12, pad=12)

    # Colorbar at the bottom (mirrors the paper's pval legend)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", pad=0.08,
                        shrink=0.7, aspect=18)
    cbar.set_label("$-\\log_{10}(p)$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    # Add a tick at 0.05 to mark the conventional significance threshold
    ticks = [0, -np.log10(0.05), 5, 10, 15, 20]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.1f}" for t in ticks])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")


def main():
    args = parse_args()

    print(f"Loading background {args.background}", file=sys.stderr)
    N, K_by_sp = load_background(args.background)
    print(f"  N = {N:,}", file=sys.stderr)
    for sp_key, c in K_by_sp.items():
        print(f"  K[{sp_key}] = {c:,}", file=sys.stderr)

    print(f"\nFiltering {args.hits} (evalue <= {args.evalue})...", file=sys.stderr)
    unique, n_rows, n_kept = filter_and_dedupe(args.hits, args.evalue)
    n = len(unique)
    print(f"  rows total:           {n_rows:,}", file=sys.stderr)
    print(f"  rows passing evalue:  {n_kept:,}", file=sys.stderr)
    print(f"  unique targets (n):   {n:,}", file=sys.stderr)

    k_counts = count_per_species(unique)
    pvals = hypergeom_pvals(k_counts, n, K_by_sp, N)

    print(f"\nHypergeometric P(X >= k):", file=sys.stderr)
    print(f"{'species':<16} {'k':>7} {'K':>9} {'n':>9} {'N':>10}  p-value", file=sys.stderr)
    for label, _, bg in SPECIES:
        K = K_by_sp.get(bg, 0)
        print(f"{label:<16} {k_counts[label]:>7,} {K:>9,} {n:>9,} {N:>10,}  "
              f"{pvals[label]:.3e}", file=sys.stderr)

    plot(pvals, k_counts, n, args.title, args.output)
    print(f"\nSaved {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
