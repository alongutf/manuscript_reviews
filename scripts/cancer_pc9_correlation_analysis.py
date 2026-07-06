"""
Correlation-spectrum analysis of PC9 cancer scRNA-seq data (day 14).

Applies the exact same GMP-Cor / correlation-spectrum pipeline used for the
bacterial datasets in this project (src.analysis_functions.get_eig_dist) to the
day-14 PC9 samples in C:\\Users\\owner\\Documents\\Cancer_single_cell\\sample_data.

Pipeline
--------
1. Gene panel filtering (shared across all samples):
     - From pc9_data/marker_genes_pc9.csv take the top 400 genes by `scores`
       for group 0 and the top 400 for group 2, then take their union
       (800 genes here -- the two groups do not overlap).
2. Cell filtering (per sample):
     - Keep the top n_genes/2 cells ranked by total UMI (library size,
       summed over all genes in the raw count matrix).
     - With 800 panel genes this yields 400 cells -> a 400 x 800 matrix.
3. Correlation spectrum (per sample), identical to the paper datasets:
     - get_eig_dist(m, norm=True, log=False, norm_method='sum', norm_sum=1)
       normalizes each cell to unit total, z-scores genes, and returns the
       empirical eigenvalue spectrum (pcs) vs. the mean scrambled spectrum
       (pcs1, 10 column-permutation repeats).
4. GMP-Cor = sum(max(lambda_i - lambda_max^scrambled, 0))
       = np.sum(pcs[pcs > pcs1.max()] - pcs1.max()).
5. Outputs: a 1-CCDF (complementary CDF) log-log plot per sample of the
   original vs. scrambled spectrum, a summary CSV of GMP-Cor, and the raw
   [pcs, pcs1] arrays.

Run from the `scripts/` directory (paths are resolved from __file__, so it also
works from elsewhere).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- make `src` importable regardless of CWD ---------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
import src.analysis_functions as af  # noqa: E402

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
DATA_DIR = r"C:\Users\owner\Documents\Cancer_single_cell\sample_data"
MARKER_FILE = os.path.join(REPO_ROOT, "pc9_data", "marker_genes_pc9.csv")
OUT_DIR = os.path.join(REPO_ROOT, "results", "cancer_pc9")

SAMPLES = [
    "14_rep1_high", "14_rep1_med", "14_rep1_low",
    "14_rep2_high", "14_rep2_med", "14_rep2_low",
]

MARKER_GROUPS = [0, 2]     # gene groups to pool
TOP_N_PER_GROUP = 400      # top genes by score per group

# get_eig_dist settings -- identical to the paper's data_for_paper pipeline
NORM = True
LOG = False
NORM_METHOD = "sum"
NORM_SUM = 1

np.random.seed(0)  # reproducible scrambling


# ----------------------------------------------------------------------------
# Gene panel
# ----------------------------------------------------------------------------
def select_gene_panel():
    """Union of the top-scoring genes from the requested marker groups."""
    mg = pd.read_csv(MARKER_FILE)
    genes = pd.Index([])
    for grp in MARKER_GROUPS:
        top = (mg[mg["group"] == grp]
               .sort_values("scores", ascending=False)["names"]
               .head(TOP_N_PER_GROUP))
        genes = genes.union(pd.Index(top))
    return genes


# ----------------------------------------------------------------------------
# Per-sample filtered matrix
# ----------------------------------------------------------------------------
def build_filtered_matrix(sample, panel):
    """Return a (cells x panel-genes) count matrix for one sample.

    Cell selection: top n_genes/2 cells by total UMI (library size over ALL
    genes in the raw matrix), then restrict columns to the marker panel.
    """
    df = pd.read_csv(os.path.join(DATA_DIR, f"{sample}.csv"), index_col=0)

    panel_present = [g for g in panel if g in df.columns]
    n_genes = len(panel_present)
    n_cells = n_genes // 2

    library_size = df.sum(axis=1)                       # total UMI per cell
    keep_cells = library_size.sort_values(ascending=False).index[:n_cells]

    sub = df.loc[keep_cells, panel_present]
    return sub.to_numpy(dtype=float), keep_cells, panel_present


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def plot_ccdf(ax, pcs, pcs1, title):
    """1 - CDF (CCDF) log-log plot of original vs. scrambled eigenvalues."""
    pcs = np.sort(pcs[pcs > 0])
    pcs1 = np.sort(pcs1[pcs1 > 0])
    ccdf = lambda n: 1 - np.arange(1, n + 1) / n + 1 / n

    thr = pcs1.max()
    signal = pcs > thr

    ax.loglog(pcs[~signal], ccdf(len(pcs))[~signal], ".", ls="-",
              color="darkgray", alpha=0.7, ms=3, label="spurious")
    ax.loglog(pcs[signal], ccdf(len(pcs))[signal], ".", ls="-",
              color="#3182bd", ms=3, label="original (signal)")
    ax.loglog(pcs1, ccdf(len(pcs1)), ".", ls="-",
              color="#de2d26", alpha=0.7, ms=3, label="scrambled")

    ax.axvline(thr, color="k", ls="--", alpha=0.6)
    ax.text(thr * 1.1, 0.8, r"$\lambda_\mathrm{max}^\mathrm{scr}$",
            transform=ax.get_xaxis_transform(), fontsize=8, va="center", alpha=0.7)

    ax.set_xlim(1e-1, pcs.max() + 10)
    ax.grid(True, which="both", ls="-", alpha=0.4)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"CCDF: $P(X>\lambda)$")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc="lower left")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "ev_data"), exist_ok=True)

    panel = select_gene_panel()
    print(f"Gene panel: {len(panel)} genes "
          f"(union of top {TOP_N_PER_GROUP} in groups {MARKER_GROUPS})")

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    records = []

    for ax, sample in zip(axes.ravel(), SAMPLES):
        m, cells, genes = build_filtered_matrix(sample, panel)
        print(f"\n{sample}: filtered matrix {m.shape[0]} cells x {m.shape[1]} genes")

        pcs, pcs1, frac_nz = af.get_eig_dist(
            m, norm=NORM, log=LOG, norm_method=NORM_METHOD, norm_sum=NORM_SUM)

        thr = pcs1.max()
        gmp_cor = float(np.sum(pcs[pcs > thr] - thr))
        n_signal = int(np.sum(pcs > thr))

        plot_ccdf(ax, pcs, pcs1, f"{sample}\nGMP-Cor = {gmp_cor:.2f}")
        np.save(os.path.join(OUT_DIR, "ev_data", f"{sample}.npy"),
                np.array([pcs, pcs1]))  # shape (2, P): row 0 empirical, row 1 scrambled

        records.append({
            "sample": sample,
            "n_cells": m.shape[0],
            "n_genes": m.shape[1],
            "max_ev": float(pcs.max()),
            "max_ev_scrambled": float(thr),
            "eig_diff": float(pcs.max() - thr),
            "n_ev_above_threshold": n_signal,
            "GMP_Cor": gmp_cor,
            "fraction_non_zero": float(frac_nz),
        })
        print(f"  GMP-Cor = {gmp_cor:.3f}  (max_ev={pcs.max():.2f}, "
              f"scrambled_max={thr:.2f}, #above={n_signal})")

    fig.suptitle("PC9 day-14 correlation spectrum (original vs. scrambled)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    svg = os.path.join(OUT_DIR, "pc9_day14_ccdf.svg")
    png = os.path.join(OUT_DIR, "pc9_day14_ccdf.png")
    fig.savefig(svg)
    fig.savefig(png, dpi=200)
    print(f"\nSaved figure -> {svg}\n              -> {png}")

    summary = pd.DataFrame(records)
    csv = os.path.join(OUT_DIR, "pc9_day14_gmp_cor.csv")
    summary.to_csv(csv, index=False)
    print(f"Saved summary -> {csv}\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
