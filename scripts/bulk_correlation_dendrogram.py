"""
Four-panel Spearman correlation figure for the 12 bulk RNA-seq samples.

The panels cross the two choices that matter for how the samples group:

                      no batch correction        batch-corrected
  log fold change            A                          B
  absolute expression        C                          D

Spearman only: it answers "do the samples rank genes the same way", is invariant
to the log transform, and is not pulled around by the long tail of large fold
changes, so the four panels are directly comparable to each other.

By default each panel carries an average-linkage dendrogram on distance 1 - rho
and the heatmap is ordered by the dendrogram leaves rather than by a fixed
grouping, so which samples cluster together is read off the tree instead of being
imposed.  With --no-dendrogram the tree is dropped and every panel uses the same
fixed condition > batch > sample order, which makes the four panels directly
comparable cell for cell at the cost of no longer showing the clustering.  Either
way the two colour strips underneath re-expose the design: condition and batch.

Input:  results/bulk_pca/bulk_log_fold_changes<SUFFIX>.csv
        results/bulk_pca/bulk_normalized_counts<SUFFIX>.csv
        written by  python bulk_lfc_pca.py --bulk-control matched --shx-late
Output: results/bulk_pca/bulk_spearman_<layout><SUFFIX>.{svg,png}
        results/bulk_pca/bulk_spearman_dendrogram_clusters<SUFFIX>.csv
        (the cluster table is written only when the dendrograms are computed)

Usage:
    python bulk_correlation_dendrogram.py
    python bulk_correlation_dendrogram.py --no-dendrogram --cmap bluered
"""

import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "bulk_pca")

parser = argparse.ArgumentParser()
parser.add_argument("--no-dendrogram", action="store_true",
                    help="drop the trees and order every panel by condition > batch > "
                         "sample instead, so the four panels are cell-for-cell comparable")
parser.add_argument("--cmap", choices=["viridis", "bluered"], default="viridis",
                    help="heatmap colours; bluered is RdBu_r, blue = weakest "
                         "correlation in the panel, red = strongest")
parser.add_argument("--suffix", default="_log2_matched_late",
                    help="suffix of the input matrices; the default is the 12-sample "
                         "matched-control / late-SHX set")
args = parser.parse_args()

# every correlation here is strongly positive, so the colour scale is stretched over
# the observed range of each panel rather than centred on zero: RdBu_r is used as a
# sequential blue -> red ramp, not as a diverging map around some neutral value
CMAP = {"viridis": "viridis", "bluered": "RdBu_r"}[args.cmap]
OUT_TAG = ("blocks" if args.no_dendrogram else "dendrogram") + (
    "_bluered" if args.cmap == "bluered" else "")


# ------------------------------------------------------------------- inputs
def load(stem):
    src = os.path.join(OUT, "%s%s.csv" % (stem, args.suffix))
    if not os.path.exists(src):
        raise SystemExit("%s not found - run bulk_lfc_pca.py --bulk-control matched "
                         "--shx-late first" % src)
    return pd.read_csv(src, index_col=0)          # genes x samples


lfc = load("bulk_log_fold_changes")
# +1 keeps the handful of zero counts finite; every gene here already passes the
# mean-count filter, so the pseudocount moves nothing that carries signal
expr = np.log10(load("bulk_normalized_counts") + 1.0)
if list(lfc.columns) != list(expr.columns):
    raise SystemExit("the two matrices do not hold the same samples")
print("%d genes x %d samples" % lfc.shape)

# ------------------------------------------------- sample annotation (as in the PCA)
LABELS = {"bulk": {"Disrupted": "Dis-Arrest1", "CASP": "Reg-Arrest1"},
          "shx": "Dis-Arrest2", "casp": "Reg-Arrest2+SHX"}


def annotate(col):
    ds, s = col.split("|")
    if ds == "bulk":
        label = LABELS["bulk"]["CASP" if s.startswith("CASP") else "Disrupted"]
    else:
        label = LABELS[ds]
    return {"dataset": ds, "sample": s, "label": label,
            "condition": "Reg-Arrest" if label.startswith("Reg") else "Dis-Arrest",
            "batch": "bulk" if ds == "bulk" else "timecourse"}


meta = pd.DataFrame([annotate(c) for c in lfc.columns], index=lfc.columns)


def short(row):
    """Compact axis label: replicate number or time point."""
    s = re.sub(r"^(CASP|Disrupted)_biorep", "rep", row["sample"])
    s = re.sub(r"^CASP biorep1 (\d+)min", r"\1 min", s)
    return "%s %s" % (row["label"], s)


NAME = {c: short(r) for c, r in meta.iterrows()}


# ------------------------------------------------------------- batch removal
def remove_batch(mat):
    """limma removeBatchEffect on a genes x samples matrix: fit
    y ~ 1 + condition + batch per gene and subtract only the fitted batch term.
    `bulk` is the reference level, so the time-course samples are shifted onto the
    bulk baseline.  Both conditions appear in both batches, so the design is a
    balanced 2x2 and the batch offsets are identified without absorbing the
    biological contrast."""
    cond = pd.get_dummies(meta["condition"], drop_first=True).astype(float)
    batch = pd.get_dummies(meta["batch"], drop_first=False).astype(float)[["timecourse"]]
    keep = np.column_stack([np.ones(len(meta)), cond.values])
    X = np.column_stack([keep, batch.values])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise SystemExit("design is rank-deficient: batch and condition are confounded")
    y = mat[meta.index].T.values                       # samples x genes
    coefs, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = batch.values @ coefs[keep.shape[1]:]
    print("  mean |batch offset| = %.3f" % np.abs(coefs[keep.shape[1]:]).mean())
    return pd.DataFrame((y - fit).T, index=mat.index, columns=meta.index)


PANELS = [("A", "Log fold change", lfc, False),
          ("B", "Log fold change, batch-corrected", lfc, True),
          ("C", "Absolute expression", expr, False),
          ("D", "Absolute expression, batch-corrected", expr, True)]

# ------------------------------------------------------------------ figure
COND_COLOR = {"Reg-Arrest": "#1f5fa9", "Dis-Arrest": "#c8102e"}
BATCH_COLOR = {"bulk": "#4d4d4d", "timecourse": "#bdbdbd"}
n = len(meta)

fig = plt.figure(figsize=(16, 14))
gs = GridSpec(2, 2, figure=fig, hspace=0.34, wspace=0.30,
              left=0.13, right=0.95, top=0.90, bottom=0.16)

# fixed layout order, used when the dendrograms are switched off
COND_ORDER = ["Reg-Arrest", "Dis-Arrest"]
BATCH_ORDER = ["bulk", "timecourse"]
FIXED_ORDER = sorted(meta.index,
                     key=lambda c: (COND_ORDER.index(meta.loc[c, "condition"]),
                                    BATCH_ORDER.index(meta.loc[c, "batch"]),
                                    meta.loc[c, "sample"]))

clusters = {}
for (tag, title, mat, debatch), cell in zip(PANELS, gs):
    print("%s: %s" % (tag, title))
    m = remove_batch(mat) if debatch else mat[meta.index]
    rho = m.corr(method="spearman")

    # average linkage on 1 - rho: the standard distance for a correlation heatmap,
    # and average linkage keeps a cluster's position tied to all of its members
    d = 1.0 - rho.values
    np.fill_diagonal(d, 0.0)
    d = (d + d.T) / 2.0                                # kill float asymmetry
    Z = linkage(squareform(d, checks=False), method="average")

    if args.no_dendrogram:
        # heatmap / two-row design strip; the panel title moves onto the heatmap
        sub = cell.subgridspec(2, 1, height_ratios=[1.0, 0.09], hspace=0.03)
        ax_h, ax_a = (fig.add_subplot(sub[i]) for i in range(2))
        axes_group = [ax_h, ax_a]
        order = FIXED_ORDER
        ax_h.set_title("%s  %s" % (tag, title), fontsize=10, pad=6, loc="left")
    else:
        # dendrogram / heatmap / two-row design strip
        sub = cell.subgridspec(3, 1, height_ratios=[0.28, 1.0, 0.09], hspace=0.03)
        ax_d, ax_h, ax_a = (fig.add_subplot(sub[i]) for i in range(3))
        axes_group = [ax_d, ax_h, ax_a]
        dend = dendrogram(Z, ax=ax_d, no_labels=True, color_threshold=0,
                          link_color_func=lambda _: "0.3")
        order = [rho.index[i] for i in dend["leaves"]]
        ax_d.set_xticks([])
        ax_d.tick_params(labelsize=6)
        for side in ("top", "right", "bottom"):
            ax_d.spines[side].set_visible(False)
        ax_d.set_ylabel("1 - rho", fontsize=6.5)
        ax_d.set_title("%s  %s" % (tag, title), fontsize=10, pad=6, loc="left")

    r = rho.loc[order, order]
    off = ~np.eye(n, dtype=bool)
    vmin, vmax = np.floor(r.values[off].min() * 20) / 20, 1.0
    # aspect="auto" so the cells fill the axes box: the dendrogram above and the
    # design strip below then line up column for column with the heatmap
    im = ax_h.imshow(r.values, cmap=CMAP, vmin=vmin, vmax=vmax, aspect="auto",
                     extent=(-0.5, n - 0.5, n - 0.5, -0.5))
    for i in range(n):
        for j in range(n):
            v = r.values[i, j]
            # viridis is dark at its low end only, RdBu_r is dark at both ends
            t = (v - vmin) / (vmax - vmin)
            dark = (t < 0.2 or t > 0.8) if args.cmap == "bluered" else (t < 0.75)
            ax_h.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=5.5,
                      color="white" if dark else "0.15")
    ax_h.set_xticks(range(n))
    ax_h.set_xticklabels([])
    ax_h.set_yticks(range(n))
    ax_h.set_yticklabels([NAME[c] for c in order], fontsize=7)
    for t, c in zip(ax_h.get_yticklabels(), order):
        t.set_color(COND_COLOR[meta.loc[c, "condition"]])
    cb = fig.colorbar(im, ax=axes_group, fraction=0.032, pad=0.015)
    cb.set_label("Spearman rho", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # design strips: condition (top) and batch (bottom), in the panel order
    strip = np.zeros((2, n, 3))
    for j, c in enumerate(order):
        strip[0, j] = to_rgb(COND_COLOR[meta.loc[c, "condition"]])
        strip[1, j] = to_rgb(BATCH_COLOR[meta.loc[c, "batch"]])
    ax_a.imshow(strip, aspect="auto", extent=(-0.5, n - 0.5, 1.5, -0.5))
    ax_a.set_yticks([0, 1])
    ax_a.set_yticklabels(["condition", "batch"], fontsize=6.5)
    ax_a.set_xticks(range(n))
    ax_a.set_xticklabels([NAME[c] for c in order], rotation=90, fontsize=7)
    for t, c in zip(ax_a.get_xticklabels(), order):
        t.set_color(COND_COLOR[meta.loc[c, "condition"]])
    ax_a.tick_params(length=0)

    # two-cluster cut: does the top split of the tree recover the condition?
    cut = pd.Series(fcluster(Z, 2, criterion="maxclust"), index=rho.index)
    clusters["%s_%s" % (tag, "nobatch" if debatch else "raw")] = cut
    tab = pd.crosstab(cut, meta["condition"])
    pure = bool((tab.max(axis=1) == tab.sum(axis=1)).all())
    print("  top split vs condition: %s" % ("clean" if pure else "mixed"))
    print("    " + tab.to_string().replace("\n", "\n    "))

handles = ([plt.Line2D([], [], marker="s", ls="", ms=8, mfc=v, mec="none", label=k)
            for k, v in COND_COLOR.items()]
           + [plt.Line2D([], [], marker="s", ls="", ms=8, mfc=v, mec="none",
                         label="batch: " + k) for k, v in BATCH_COLOR.items()])
fig.legend(handles=handles, ncol=4, fontsize=8.5, frameon=False,
           loc="lower center", bbox_to_anchor=(0.5, 0.03))
fig.suptitle("Bulk RNA-seq: Spearman sample-sample correlation, %d genes x %d samples\n%s"
             % (lfc.shape[0], n,
                "samples ordered by condition > batch" if args.no_dendrogram
                else "average-linkage clustering on 1 - rho"), fontsize=12)

if not args.no_dendrogram:
    pd.DataFrame(clusters).to_csv(
        os.path.join(OUT, "bulk_spearman_dendrogram_clusters%s.csv" % args.suffix))
stem = os.path.join(OUT, "bulk_spearman_%s%s" % (OUT_TAG, args.suffix))
fig.savefig(stem + ".svg")
fig.savefig(stem + ".png", dpi=200)
print("\nwrote", stem + ".svg")
