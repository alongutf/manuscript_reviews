"""
Single-panel version of panel D of bulk_correlation_dendrogram.py: the Spearman
sample-sample correlation of the 12 bulk RNA-seq samples in *absolute expression*
space, after batch correction.

The four-panel figure crosses two choices (log fold change vs absolute expression,
raw vs batch-corrected).  This script keeps only the one panel, at a size where the
tree, the labels and the cell numbers are readable, and adds two things the
four-panel figure does not carry:

  1. validation - the plotted matrix is recomputed three independent ways
     (pandas .corr, scipy.stats.spearmanr, Pearson on explicit ranks), the batch
     correction is re-derived one gene at a time by an independent least-squares
     fit, and the raw (uncorrected) matrix is compared cell for cell against
     results/bulk_pca/bulk_expression_correlation_spearman<SUFFIX>.csv written by
     bulk_lfc_correlation.py.  Spearman is invariant to the monotone log10(x + 1),
     so the untransformed normalised counts must give the same raw matrix as well -
     that is checked too.
  2. pair scatters - a handful of sample pairs spanning the observed range of rho
     are drawn gene by gene, in value space and in rank space, with the rho of the
     plotted points annotated.  Each annotated rho is asserted to equal the
     corresponding heatmap cell, so the scatters verify the matrix rather than
     merely illustrating it.

Input:  results/bulk_pca/bulk_normalized_counts<SUFFIX>.csv   (genes x samples)
        written by  python bulk_lfc_pca.py --bulk-control matched --shx-late
Output: results/bulk_pca/bulk_expression_spearman_nobatch<SUFFIX>.{svg,png}
        results/bulk_pca/bulk_expression_spearman_nobatch<SUFFIX>.csv
        results/bulk_pca/bulk_expression_pair_scatter<SUFFIX>.{svg,png}
        results/bulk_pca/bulk_expression_spearman_nobatch_validation<SUFFIX>.txt

Usage:
    python bulk_expression_correlation_panel.py
    python bulk_expression_correlation_panel.py --no-dendrogram
    python bulk_expression_correlation_panel.py --cmap bluered --n-pairs 6
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
from scipy.stats import spearmanr, rankdata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "bulk_pca")

parser = argparse.ArgumentParser()
parser.add_argument("--no-dendrogram", action="store_true",
                    help="drop the tree and order the samples condition > batch > sample")
parser.add_argument("--cmap", choices=["viridis", "bluered"], default="viridis",
                    help="heatmap colours; bluered is RdBu_r used as a sequential ramp")
parser.add_argument("--n-pairs", type=int, default=4,
                    help="how many sample pairs to draw as validation scatters")
parser.add_argument("--suffix", default="_log2_matched_late",
                    help="suffix of the input matrix; the default is the 12-sample "
                         "matched-control / late-SHX set")
args = parser.parse_args()

CMAP = {"viridis": "viridis", "bluered": "RdBu_r"}[args.cmap]
# TAG names the figure/validation-log files (distinct from the CSV, which is keyed
# only by --suffix so bulk_correlation_dendrogram.py's reference CSV check finds it
# regardless of which colour map produced this run's plots)
TAG = args.suffix + ("_bluered" if args.cmap == "bluered" else "")
TOL = 1e-10                      # every check below is exact to float noise
log = []


def say(msg=""):
    """Print a line and also buffer it, so the same text ends up on stdout and
    in the validation .txt file written at the end of the script."""
    print(msg)
    log.append(msg)


# ------------------------------------------------------------------------- input
src = os.path.join(OUT, "bulk_normalized_counts%s.csv" % args.suffix)
if not os.path.exists(src):
    raise SystemExit("%s not found - run bulk_lfc_pca.py --bulk-control matched "
                     "--shx-late first" % src)
counts = pd.read_csv(src, index_col=0)                 # genes x samples
# +1 keeps the handful of zero counts finite; every gene here already passes the
# mean-count filter, so the pseudocount moves nothing that carries signal
expr = np.log10(counts + 1.0)
say("input   %s" % os.path.basename(src))
say("        %d genes x %d samples, log10(normalised count + 1)" % expr.shape)

# --------------------------------------------------- sample annotation (as in the PCA)
LABELS = {"bulk": {"Disrupted": "Dis-Arrest1", "CASP": "Reg-Arrest1"},
          "shx": "Dis-Arrest2", "casp": "Reg-Arrest2+SHX"}


def annotate(col):
    """Derive one metadata row (dataset, sample, publication label, condition,
    batch) from a "<dataset>|<sample>" column name of the input matrix."""
    ds, s = col.split("|")
    if ds == "bulk":
        label = LABELS["bulk"]["CASP" if s.startswith("CASP") else "Disrupted"]
    else:
        label = LABELS[ds]
    return {"dataset": ds, "sample": s, "label": label,
            "condition": "Reg-Arrest" if label.startswith("Reg") else "Dis-Arrest",
            "batch": "bulk" if ds == "bulk" else "timecourse"}


meta = pd.DataFrame([annotate(c) for c in expr.columns], index=expr.columns)


def short(row):
    """Compact axis label: replicate number or time point."""
    s = re.sub(r"^(CASP|Disrupted)_biorep", "rep", row["sample"])
    s = re.sub(r"^CASP biorep1 (\d+)min", r"\1 min", s)
    return "%s %s" % (row["label"], s)


NAME = {c: short(r) for c, r in meta.iterrows()}
n = len(meta)


# ----------------------------------------------------------------- batch removal
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
    return pd.DataFrame((y - fit).T, index=mat.index, columns=meta.index), X, coefs


corrected, X, coefs = remove_batch(expr)
offsets = coefs[2]                                     # per-gene timecourse offset
say("batch   design rank %d / %d columns; per-gene timecourse offset: "
    "mean %+.3f, mean |.| %.3f, sd %.3f"
    % (np.linalg.matrix_rank(X), X.shape[1], offsets.mean(),
       np.abs(offsets).mean(), offsets.std(ddof=0)))

rho = corrected.corr(method="spearman")
rho.to_csv(os.path.join(OUT, "bulk_expression_spearman_nobatch%s.csv" % args.suffix))

# ==================================================================== validation
say()
say("--- validation " + "-" * 63)
checks = []


def check(name, ok, detail):
    """Record one validation result and print/log it immediately; `checks`
    accumulates (name, passed) pairs so the pass/fail tally can be printed
    and the process can exit non-zero if anything failed."""
    checks.append((name, bool(ok)))
    say("  [%s] %-46s %s" % ("ok" if ok else "FAIL", name, detail))


# 1. scipy's own Spearman on the same matrix
sp = spearmanr(corrected.values)[0]
check("scipy.stats.spearmanr", np.abs(sp - rho.values).max() < TOL,
      "max |diff| = %.2e" % np.abs(sp - rho.values).max())

# 2. Spearman is Pearson on the ranks: build the ranks explicitly, use np.corrcoef
ranks = np.apply_along_axis(rankdata, 0, corrected.values)
pr = np.corrcoef(ranks, rowvar=False)
check("Pearson on explicit per-sample ranks", np.abs(pr - rho.values).max() < TOL,
      "max |diff| = %.2e" % np.abs(pr - rho.values).max())

# 3. the batch correction, re-derived one gene at a time by an independent solve
rng = np.random.default_rng(0)
probe = rng.choice(expr.shape[0], size=200, replace=False)
redo = np.empty((len(probe), n))
for k, g in enumerate(probe):
    y = expr[meta.index].values[g]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    redo[k] = y - X[:, 2] * b[2]
gap = np.abs(redo - corrected.values[probe]).max()
check("per-gene refit of the batch term (200 genes)", gap < 1e-9,
      "max |diff| = %.2e" % gap)

# 4. only the batch term was removed: the condition contrast must survive untouched.
#    The design is a balanced 2x2, so both group means shift by the same amount and
#    the per-gene Reg - Dis difference is preserved exactly.
d_before = (expr[meta.index[meta["condition"] == "Reg-Arrest"]].mean(axis=1)
            - expr[meta.index[meta["condition"] == "Dis-Arrest"]].mean(axis=1))
d_after = (corrected[meta.index[meta["condition"] == "Reg-Arrest"]].mean(axis=1)
           - corrected[meta.index[meta["condition"] == "Dis-Arrest"]].mean(axis=1))
gap = np.abs(d_before - d_after).max()
check("condition contrast preserved by the correction", gap < 1e-9,
      "max |diff| in per-gene Reg - Dis = %.2e" % gap)

# 5. the raw (uncorrected) matrix against the CSV bulk_lfc_correlation.py wrote
raw = expr[meta.index].corr(method="spearman")
ref_path = os.path.join(OUT, "bulk_expression_correlation_spearman%s.csv" % args.suffix)
if os.path.exists(ref_path):
    ref = pd.read_csv(ref_path, index_col=0)
    common = [c for c in ref.columns if c in raw.columns]
    gap = np.abs(raw.loc[common, common].values - ref.loc[common, common].values).max()
    check("raw matrix vs bulk_lfc_correlation.py CSV", gap < 5e-6,
          "%d x %d cells, max |diff| = %.2e" % (len(common), len(common), gap))
else:
    say("  [--] reference CSV not found, skipping the cross-script check")

# 6. Spearman is invariant to the monotone log10(x + 1), so the untransformed
#    normalised counts must give the same raw matrix
raw_counts = counts[meta.index].corr(method="spearman")
check("Spearman invariant to log10(x + 1)",
      np.abs(raw_counts.values - raw.values).max() < TOL,
      "max |diff| = %.2e" % np.abs(raw_counts.values - raw.values).max())

# 7. shape sanity: symmetric, unit diagonal, inside [-1, 1]
v = rho.values
off = ~np.eye(n, dtype=bool)
check("symmetric, unit diagonal, within [-1, 1]",
      np.abs(v - v.T).max() < TOL and np.abs(np.diag(v) - 1).max() < TOL
      and v.min() >= -1 and v.max() <= 1,
      "off-diagonal range %.4f - %.4f" % (v[off].min(), v[off].max()))

# what the correction actually did to the correlations
same_c = (meta["condition"].values[:, None] == meta["condition"].values[None, :])
same_b = (meta["batch"].values[:, None] == meta["batch"].values[None, :])
say()
say("  mean off-diagonal rho, before -> after batch correction")
for lab, mask in [("same condition, same batch", off & same_c & same_b),
                  ("same condition, diff batch", off & same_c & ~same_b),
                  ("diff condition, same batch", off & ~same_c & same_b),
                  ("diff condition, diff batch", off & ~same_c & ~same_b)]:
    say("    %-27s %.4f -> %.4f  (n=%d)"
        % (lab, raw.values[mask].mean(), v[mask].mean(), mask.sum() // 2))

# ======================================================================== figure
COND_COLOR = {"Reg-Arrest": "#1f5fa9", "Dis-Arrest": "#c8102e"}
BATCH_COLOR = {"bulk": "#4d4d4d", "timecourse": "#bdbdbd"}
COND_ORDER = ["Reg-Arrest", "Dis-Arrest"]
BATCH_ORDER = ["bulk", "timecourse"]
FIXED_ORDER = sorted(meta.index,
                     key=lambda c: (COND_ORDER.index(meta.loc[c, "condition"]),
                                    BATCH_ORDER.index(meta.loc[c, "batch"]),
                                    meta.loc[c, "sample"]))

# average linkage on 1 - rho: the standard distance for a correlation heatmap, and
# average linkage keeps a cluster's position tied to all of its members
d = 1.0 - v.copy()
np.fill_diagonal(d, 0.0)
d = (d + d.T) / 2.0                                    # kill float asymmetry
Z = linkage(squareform(d, checks=False), method="average")

fig = plt.figure(figsize=(8.6, 9.4))
if args.no_dendrogram:
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.0, 0.09], hspace=0.03,
                  left=0.30, right=0.88, top=0.90, bottom=0.24)
    ax_h, ax_a = (fig.add_subplot(gs[i]) for i in range(2))
    axes_group = [ax_h, ax_a]
    order = FIXED_ORDER
else:
    gs = GridSpec(3, 1, figure=fig, height_ratios=[0.28, 1.0, 0.09], hspace=0.03,
                  left=0.30, right=0.88, top=0.90, bottom=0.24)
    ax_d, ax_h, ax_a = (fig.add_subplot(gs[i]) for i in range(3))
    axes_group = [ax_d, ax_h, ax_a]
    dend = dendrogram(Z, ax=ax_d, no_labels=True, color_threshold=0,
                      link_color_func=lambda _: "0.3")
    order = [rho.index[i] for i in dend["leaves"]]
    ax_d.set_xticks([])
    ax_d.tick_params(labelsize=6)
    for side in ("top", "right", "bottom"):
        ax_d.spines[side].set_visible(False)
    ax_d.set_ylabel("1 - rho", fontsize=7)

r = rho.loc[order, order]
vmin, vmax = np.floor(r.values[off].min() * 20) / 20, 1.0
# aspect="auto" so the cells fill the axes box: the dendrogram above and the design
# strip below then line up column for column with the heatmap
im = ax_h.imshow(r.values, cmap=CMAP, vmin=vmin, vmax=vmax, aspect="auto",
                 extent=(-0.5, n - 0.5, n - 0.5, -0.5))
for i in range(n):
    for j in range(n):
        val = r.values[i, j]
        # viridis is dark at its low end only, RdBu_r is dark at both ends
        t = (val - vmin) / (vmax - vmin)
        dark = (t < 0.2 or t > 0.8) if args.cmap == "bluered" else (t < 0.75)
        ax_h.text(j, i, "%.2f" % val, ha="center", va="center", fontsize=6.5,
                  color="white" if dark else "0.15")
ax_h.set_xticks(range(n))
ax_h.set_xticklabels([])
ax_h.set_yticks(range(n))
ax_h.set_yticklabels([NAME[c] for c in order], fontsize=8)
for t, c in zip(ax_h.get_yticklabels(), order):
    t.set_color(COND_COLOR[meta.loc[c, "condition"]])
cb = fig.colorbar(im, ax=axes_group, fraction=0.045, pad=0.02)
cb.set_label("Spearman rho", fontsize=8.5)
cb.ax.tick_params(labelsize=7.5)

# design strips: condition (top) and batch (bottom), in the panel order
strip = np.zeros((2, n, 3))
for j, c in enumerate(order):
    strip[0, j] = to_rgb(COND_COLOR[meta.loc[c, "condition"]])
    strip[1, j] = to_rgb(BATCH_COLOR[meta.loc[c, "batch"]])
ax_a.imshow(strip, aspect="auto", extent=(-0.5, n - 0.5, 1.5, -0.5))
ax_a.set_yticks([0, 1])
ax_a.set_yticklabels(["condition", "batch"], fontsize=7)
ax_a.set_xticks(range(n))
ax_a.set_xticklabels([NAME[c] for c in order], rotation=90, fontsize=8)
for t, c in zip(ax_a.get_xticklabels(), order):
    t.set_color(COND_COLOR[meta.loc[c, "condition"]])
ax_a.tick_params(length=0)

handles = ([plt.Line2D([], [], marker="s", ls="", ms=8, mfc=cv, mec="none", label=k)
            for k, cv in COND_COLOR.items()]
           + [plt.Line2D([], [], marker="s", ls="", ms=8, mfc=cv, mec="none",
                         label="batch: " + k) for k, cv in BATCH_COLOR.items()])
fig.legend(handles=handles, ncol=4, fontsize=8.5, frameon=False,
           loc="lower center", bbox_to_anchor=(0.5, 0.015))
fig.suptitle("Bulk RNA-seq: Spearman sample-sample correlation\n"
             "absolute expression, batch-corrected  (%d genes x %d samples, %s)"
             % (expr.shape[0], n,
                "condition > batch order" if args.no_dendrogram
                else "average-linkage clustering on 1 - rho"), fontsize=11)

# two-cluster cut: does the top split of the tree recover the condition?
cut = pd.Series(fcluster(Z, 2, criterion="maxclust"), index=rho.index)
tab = pd.crosstab(cut, meta["condition"])
say()
say("  top split of the tree vs condition: %s"
    % ("clean" if bool((tab.max(axis=1) == tab.sum(axis=1)).all()) else "mixed"))
say("    " + tab.to_string().replace("\n", "\n    "))

stem = os.path.join(OUT, "bulk_expression_spearman_nobatch%s" % TAG)
fig.savefig(stem + ".svg")
fig.savefig(stem + ".png", dpi=200)

# ===================================================== pair scatters (matrix check)
# pairs spanning the observed range of rho: the strongest and the weakest pair plus
# evenly spaced ones in between, so the scatters cover the whole colour scale
iu = np.triu_indices(n, 1)
vals = v[iu]
# n_pairs indices evenly spaced (by rank) between the weakest and strongest pair;
# np.unique both sorts them and drops duplicate ranks that appear when n_pairs is
# large relative to the number of distinct pairs
picks = np.unique(np.linspace(0, len(vals) - 1, args.n_pairs).round().astype(int))
sel = np.argsort(vals)[::-1][picks]                    # strongest first

k = len(sel)
fig2, axes2 = plt.subplots(2, k, figsize=(3.1 * k, 6.6), squeeze=False)
say()
say("--- pair scatters " + "-" * 60)
for col, s in enumerate(sel):
    a, b = meta.index[iu[0][s]], meta.index[iu[1][s]]
    xa, yb = corrected[a].values, corrected[b].values
    rho_pair = spearmanr(xa, yb)[0]
    r_pair = np.corrcoef(xa, yb)[0, 1]
    # the scatter is drawn from the same batch-corrected columns the heatmap used,
    # so its rho must reproduce the heatmap cell exactly - assert, do not trust
    assert abs(rho_pair - rho.loc[a, b]) < 1e-9, (a, b, rho_pair, rho.loc[a, b])
    say("  %-26s vs %-26s  heatmap rho = %.4f, scatter rho = %.4f, Pearson r = %.4f"
        % (NAME[a], NAME[b], rho.loc[a, b], rho_pair, r_pair))

    same = meta.loc[a, "condition"] == meta.loc[b, "condition"]
    color = COND_COLOR[meta.loc[a, "condition"]] if same else "0.35"

    ax = axes2[0][col]
    ax.scatter(xa, yb, s=2.5, alpha=0.25, lw=0, color=color, rasterized=True)
    lo, hi = min(xa.min(), yb.min()), max(xa.max(), yb.max())
    ax.plot([lo, hi], [lo, hi], color="k", lw=0.7, ls="--", alpha=0.6)
    ax.set_xlabel(NAME[a], fontsize=8)
    ax.set_ylabel(NAME[b], fontsize=8)
    ax.set_title("rho = %.3f\nr = %.3f" % (rho_pair, r_pair), fontsize=9)
    ax.tick_params(labelsize=7)

    # rank space: this is literally what Spearman correlates, so the tightness of
    # this cloud is the heatmap cell
    ax = axes2[1][col]
    ax.scatter(rankdata(xa), rankdata(yb), s=2.5, alpha=0.25, lw=0, color=color,
               rasterized=True)
    ax.plot([1, len(xa)], [1, len(xa)], color="k", lw=0.7, ls="--", alpha=0.6)
    ax.set_xlabel("rank, " + NAME[a], fontsize=8)
    ax.set_ylabel("rank, " + NAME[b], fontsize=8)
    ax.tick_params(labelsize=7)

for ax in axes2.ravel():
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
fig2.suptitle("Gene-by-gene scatter for %d sample pairs spanning the heatmap range\n"
              "top: batch-corrected log10(normalised count + 1);  bottom: the ranks "
              "Spearman actually correlates  (%d genes)" % (k, expr.shape[0]),
              fontsize=10.5)
fig2.tight_layout(rect=(0, 0, 1, 0.93))
stem2 = os.path.join(OUT, "bulk_expression_pair_scatter%s" % TAG)
fig2.savefig(stem2 + ".svg")
fig2.savefig(stem2 + ".png", dpi=200)

# ========================================================================= report
say()
failed = [nm for nm, ok in checks if not ok]
say("%d/%d checks passed%s" % (len(checks) - len(failed), len(checks),
                               "" if not failed else "  FAILED: " + ", ".join(failed)))
say()
say("wrote %s.{svg,png}" % stem)
say("wrote %s.{svg,png}" % stem2)
say("wrote %s" % os.path.join(OUT, "bulk_expression_spearman_nobatch%s.csv" % args.suffix))

val_path = os.path.join(OUT, "bulk_expression_spearman_nobatch_validation%s.txt" % TAG)
with open(val_path, "w") as fh:
    fh.write("bulk_expression_correlation_panel.py  %s\nargs: %s\n\n"
             % (pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"), vars(args)))
    fh.write("\n".join(log) + "\n")
print("wrote", val_path)
if failed:
    raise SystemExit(1)
