"""
Sample-sample correlation of the bulk RNA-seq profiles, in two spaces.

Companion to bulk_lfc_pca.py / bulk_pca.py: the same 12-sample matrices that go
into the PCAs are here summarised as correlation matrices, which show pairwise
sample agreement directly instead of through a projection.

  --space lfc         log2 fold change of each sample vs its own exponential
                      control.  The baseline is divided out, so what is left is
                      the response to the treatment.
  --space expression  absolute expression, log10(normalised count + 1).  Here the
                      shared baseline dominates: every sample is mostly "an E.
                      coli transcriptome", so all correlations sit near 1 and the
                      structure to look for is the small residual spread.

Both correlations are reported because they answer different questions:
  Pearson  - do the samples agree on the *size* of the values?
  Spearman - do they agree on the *ranking* of genes, robust to the long tail of
             large values.  Spearman is invariant to the log, so the expression
             panel's Spearman is the correlation of the raw normalised counts.

Samples are ordered by condition, then by batch within condition, so the
condition blocks and the batch sub-blocks are both visible on the diagonal.

Input:  results/bulk_pca/bulk_log_fold_changes<SUFFIX>.csv   (--space lfc)
        results/bulk_pca/bulk_normalized_counts<SUFFIX>.csv  (--space expression)
        both genes x samples, written by
        python bulk_lfc_pca.py --bulk-control matched --shx-late
Output: results/bulk_pca/bulk_<space>_correlation_{pearson,spearman}<SUFFIX>.csv
        results/bulk_pca/bulk_<space>_correlation<SUFFIX>.{svg,png}

Usage:
    python bulk_lfc_correlation.py
    python bulk_lfc_correlation.py --space expression
    python bulk_lfc_correlation.py --space expression --remove-batch
    python bulk_lfc_correlation.py --suffix _log2_matched_late_nobatch
"""

import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "bulk_pca")

parser = argparse.ArgumentParser()
parser.add_argument("--space", choices=["lfc", "expression"], default="lfc",
                    help="lfc: log2 fold change vs own exponential control; "
                         "expression: absolute log10(normalised count + 1)")
parser.add_argument("--remove-batch", action="store_true",
                    help="subtract per-gene batch offsets before correlating, the same "
                         "limma removeBatchEffect fit the PCA scripts use: the biogroup "
                         "contrast is kept, only the fitted batch term is removed")
parser.add_argument("--suffix", default="_log2_matched_late",
                    help="suffix of the input matrix to read; the default is the "
                         "12-sample matched-control / late-SHX matrix")
args = parser.parse_args()

if args.remove_batch and "_nobatch" in args.suffix:
    raise SystemExit("--suffix already names a batch-corrected matrix; drop --remove-batch")
OUT_SUFFIX = args.suffix + ("_nobatch" if args.remove_batch else "")

STEM = {"lfc": "bulk_log_fold_changes", "expression": "bulk_normalized_counts"}
src = os.path.join(OUT, "%s%s.csv" % (STEM[args.space], args.suffix))
if not os.path.exists(src):
    raise SystemExit("%s not found - run bulk_lfc_pca.py with the matching options first"
                     % src)
lfc = pd.read_csv(src, index_col=0)          # genes x samples
if args.space == "expression":
    # +1 keeps the handful of zero counts finite; every gene here already passes the
    # mean-count filter, so the pseudocount moves nothing that carries signal
    lfc = np.log10(lfc + 1.0)
print("%d genes x %d samples from %s" % (lfc.shape + (os.path.basename(src),)))

# ------------------------------------------------- sample annotation (as in the PCA)
LABELS = {"bulk": {"Disrupted": "Dis-Arrest1", "CASP": "Reg-Arrest1"},
          "shx": "Dis-Arrest2", "casp": "Reg-Arrest2+SHX"}


def annotate(col):
    ds, s = col.split("|")
    if ds == "bulk":
        cond = "CASP" if s.startswith("CASP") else "Disrupted"
        label = LABELS["bulk"][cond]
    else:
        label = LABELS[ds]
    return {"dataset": ds, "sample": s, "label": label,
            "condition": "Reg-Arrest" if label.startswith("Reg") else "Dis-Arrest",
            "batch": "bulk" if ds == "bulk" else "timecourse"}


meta = pd.DataFrame([annotate(c) for c in lfc.columns], index=lfc.columns)

# ------------------------------------------------- order: condition > batch > sample
COND_ORDER = ["Reg-Arrest", "Dis-Arrest"]
BATCH_ORDER = ["bulk", "timecourse"]
order = sorted(lfc.columns,
               key=lambda c: (COND_ORDER.index(meta.loc[c, "condition"]),
                              BATCH_ORDER.index(meta.loc[c, "batch"]),
                              meta.loc[c, "sample"]))
lfc, meta = lfc[order], meta.loc[order]


# ------------------------------------------------------------- batch removal
if args.remove_batch:
    # limma removeBatchEffect: fit  y ~ 1 + condition + batch  per gene and subtract
    # only the fitted batch term.  `bulk` is the reference level, so the time-course
    # samples are shifted onto the bulk baseline.  Both conditions appear in both
    # batches, so the design is a balanced 2x2 and the batch offsets are identified
    # without absorbing the biological contrast.
    cond = pd.get_dummies(meta["condition"], drop_first=True).astype(float)
    batch = pd.get_dummies(meta["batch"], drop_first=False).astype(float)[["timecourse"]]
    keep_cols = np.column_stack([np.ones(len(meta)), cond.values])
    X = np.column_stack([keep_cols, batch.values])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise SystemExit("design is rank-deficient: batch and condition are confounded")
    y = lfc.T.values                                   # samples x genes
    coefs, *_ = np.linalg.lstsq(X, y, rcond=None)
    batch_fit = batch.values @ coefs[keep_cols.shape[1]:]
    lfc = pd.DataFrame((y - batch_fit).T, index=lfc.index, columns=lfc.columns)
    print("removed batch offsets (rank %d design, %d coefficients); "
          "mean |batch offset| = %.3f"
          % (np.linalg.matrix_rank(X), X.shape[1],
             np.abs(coefs[keep_cols.shape[1]:]).mean()))


def short(row):
    """Compact axis label: replicate number or time point."""
    s = re.sub(r"^(CASP|Disrupted)_biorep", "rep", row["sample"])
    s = re.sub(r"^CASP biorep1 (\d+)min", r"\1'", s)
    return "%s %s" % (row["label"], s)


names = [short(r) for _, r in meta.iterrows()]

# ------------------------------------------------------------------ correlations
cors = {"pearson": lfc.corr(method="pearson"),
        "spearman": lfc.corr(method="spearman")}
for k, c in cors.items():
    c.to_csv(os.path.join(OUT, "bulk_%s_correlation_%s%s.csv" % (args.space, k, OUT_SUFFIX)))

# within/between summaries: the numbers the heatmap is meant to make visible
off = ~np.eye(len(order), dtype=bool)
for k, c in cors.items():
    v = c.values
    same_c = (meta["condition"].values[:, None] == meta["condition"].values[None, :])
    same_b = (meta["batch"].values[:, None] == meta["batch"].values[None, :])
    print("\n%s" % k)
    print("  same condition, same batch      mean r = %.3f (n=%d)"
          % (v[off & same_c & same_b].mean(), (off & same_c & same_b).sum() // 2))
    print("  same condition, different batch mean r = %.3f (n=%d)"
          % (v[off & same_c & ~same_b].mean(), (off & same_c & ~same_b).sum() // 2))
    print("  different condition             mean r = %.3f (n=%d)"
          % (v[off & ~same_c].mean(), (off & ~same_c).sum() // 2))

# ------------------------------------------------------------------------ plot
def _boundaries(meta):
    """Row indices where a new block starts, flagged True when the condition changes."""
    out = []
    for i in range(1, len(meta)):
        prev, cur = meta.iloc[i - 1], meta.iloc[i]
        if prev["condition"] != cur["condition"]:
            out.append((i, True))
        elif prev["batch"] != cur["batch"]:
            out.append((i, False))
    return out


QUANTITY = {"lfc": "log2 fold changes",
            "expression": "log10 normalised counts"}[args.space]
SPACE_NAME = {"lfc": "log fold-change", "expression": "absolute-expression"}[args.space]

n = len(order)
fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.4))
# every pair correlates positively, so a diverging scale centred on 0 would put
# the whole matrix in one colour: stretch the scale over the observed range instead
vmin = np.floor(min(c.values[off].min() for c in cors.values()) * 20) / 20
vmax = 1.0

for ax, (k, c) in zip(axes, cors.items()):
    im = ax.imshow(c.values, cmap="viridis", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_title("%s correlation of %s" % (k.capitalize(), QUANTITY), fontsize=10)
    for i in range(n):
        for j in range(n):
            val = c.values[i, j]
            ax.text(j, i, "%.2f" % val, ha="center", va="center", fontsize=5.5,
                    color="0.15" if val > vmin + 0.75 * (vmax - vmin) else "white")
    # block separators: thick = condition, thin = batch within condition
    for edge, cond_edge in _boundaries(meta):
        ax.axhline(edge - 0.5, color="w", lw=2.0 if cond_edge else 0.9)
        ax.axvline(edge - 0.5, color="w", lw=2.0 if cond_edge else 0.9)
    ax.set_xticks(np.arange(n + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n + 1) - 0.5, minor=True)
    ax.tick_params(which="minor", length=0)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label("correlation", fontsize=8)

fig.suptitle("Bulk RNA-seq: sample-sample correlation in %s space "
             "(%d genes, grouped by condition > batch%s)"
             % (SPACE_NAME, lfc.shape[0], ", batch-corrected" if args.remove_batch else ""),
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(OUT, "bulk_%s_correlation%s.svg" % (args.space, OUT_SUFFIX)))
fig.savefig(os.path.join(OUT, "bulk_%s_correlation%s.png" % (args.space, OUT_SUFFIX)), dpi=200)
print("\nwrote", os.path.join(OUT, "bulk_%s_correlation%s.svg" % (args.space, OUT_SUFFIX)))
