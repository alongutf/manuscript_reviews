"""
PCA of bulk RNA-seq samples in *log fold-change* space.

Instead of clustering samples on their absolute expression profiles (bulk_pca.py),
every sample is expressed relative to the exponential-phase control of its own
experiment / biological replicate.  This removes the experiment-specific baseline
by construction, so no batch correction is applied here.

Datasets and their exponential controls
  1. bulk_data/bulk_count_data.csv
        technical replicates  *_biorep1a/b/c  are averaged into  *_biorep1
        CASP_biorepN, Disrupted_biorepN   ->  EXP_biorepN        (N = 1,2,3)
  2. bulk_data/time_in_shx_count_data.csv
        Disrupted_biorep1t1..t8           ->  EXP_biorep1t0
  3. bulk_data/time_in_casp_count_data.csv
        CASP biorep1 400/511/621min       ->  EXP biorep2

Pipeline: harmonise gene names -> intersect genes -> drop ERCC spike-ins ->
[optionally drop rRNA/tRNA/ncRNA] -> DESeq2 median-of-ratios size factors ->
expression transform (log2 or VST) -> subtract matched control -> PCA.

The `--transform` switch exists to answer "is the VST needed here?"; running
`--compare-transforms` prints the diagnostics that settle it.

Usage:
    python bulk_lfc_pca.py                        # log2, all shared genes
    python bulk_lfc_pca.py --transform vst
    python bulk_lfc_pca.py --compare-transforms   # log2 vs VST diagnostics
    python bulk_lfc_pca.py --drop-ncrna
"""

import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "bulk_data")
OUT = os.path.join(ROOT, "results", "bulk_pca")
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument("--transform", choices=["log2", "vst"], default="log2",
                    help="per-gene transform applied before subtracting the control")
parser.add_argument("--compare-transforms", action="store_true",
                    help="run both transforms and print the diagnostics that decide "
                         "whether the VST changes the answer")
parser.add_argument("--drop-ncrna", action="store_true",
                    help="remove rRNA, tRNA and other ncRNA features")
parser.add_argument("--min-count", type=float, default=10.0,
                    help="drop genes whose mean normalised count is below this "
                         "(low counts make log fold changes explode)")
parser.add_argument("--bulk-control", choices=["matched", "mean"], default="matched",
                    help="reference for the bulk experiment: 'matched' pairs biorep n "
                         "to EXP_biorep n; 'mean' divides every sample by the average "
                         "of EXP_biorep1/2/3 instead")
parser.add_argument("--shx-late", action="store_true",
                    help="keep only SHX time points t6-t8, the subset used by bulk_pca.py")
parser.add_argument("--remove-batch", action="store_true",
                    help="subtract per-gene batch offsets from the fold-change matrix, "
                         "protecting the CASP-like / Disrupted-like contrast")
parser.add_argument("--n-hvg", type=int, default=0,
                    help="restrict the PCA to the N most variable log fold changes "
                         "(0 = use all genes)")
args = parser.parse_args()
SUFFIX = ("_" + args.transform
          + ("_coding" if args.drop_ncrna else "")
          + "_" + args.bulk_control
          + ("_late" if args.shx_late else "")
          + ("_nobatch" if args.remove_batch else "")
          + ("_hvg%d" % args.n_hvg if args.n_hvg else ""))


def noncoding_ids(gtf=os.path.join(ROOT, "metadata", "genomic.gtf")):
    """All identifiers (symbol, gene_id, locus_tag) of non-coding features in the
    MG1655 reference annotation: rRNA, tRNA and ncRNA biotypes."""
    ids = set()
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            bt = re.search(r'gene_biotype "([^"]+)"', f[8])
            if bt is None or bt.group(1) not in ("rRNA", "tRNA", "ncRNA"):
                continue
            for pat in (r'gene_id "([^"]+)"', r'gene "([^"]+)"', r'locus_tag "([^"]+)"'):
                m = re.search(pat, f[8])
                if m:
                    ids.add(m.group(1))
    return ids


def clean_ids(df):
    """Harmonise gene-name conventions across the three count tables."""
    idx = (df.index.astype(str)
           .str.strip()
           .str.replace(r"^rna-", "", regex=True)     # casp table prefixes ncRNAs with 'rna-'
           .str.replace(r"^gene-", "", regex=True))
    df = df.copy()
    df.index = idx
    # a few symbols appear twice after stripping prefixes -> sum their counts
    return df.groupby(level=0).sum()


# ---------------------------------------------------------------- load
bulk = clean_ids(pd.read_csv(os.path.join(DATA, "bulk_count_data.csv"), index_col=0))
shx = clean_ids(pd.read_csv(os.path.join(DATA, "time_in_shx_count_data.csv"), index_col=0))
casp = clean_ids(pd.read_csv(os.path.join(DATA, "time_in_casp_count_data.csv"), index_col=0))

# --------------------------------------- collapse bulk technical replicates
# 1a/1b/1c are three sequencing runs of the same biological sample.  Averaging the
# raw counts is the same as summing them up to a factor of 3, which the size-factor
# normalisation absorbs, so this matches DESeq2's collapseReplicates.
tech = {}
for c in bulk.columns:
    m = re.match(r"^(.*_biorep\d+)[abc]$", c)
    tech.setdefault(m.group(1) if m else c, []).append(c)
bulk = pd.DataFrame({k: bulk[v].mean(axis=1) for k, v in tech.items()})
print("bulk columns after collapsing technical replicates:", list(bulk.columns))

bulk.columns = ["bulk|" + c for c in bulk.columns]
shx.columns = ["shx|" + c for c in shx.columns]
casp.columns = ["casp|" + c for c in casp.columns]

# ------------------------------------------------- intersect + drop ERCC
genes = bulk.index.intersection(shx.index).intersection(casp.index)
genes = genes[~genes.str.upper().str.startswith("ERCC")]
if args.drop_ncrna:
    nc = noncoding_ids()
    # rrsX/rrlX/rrfX (rRNA) and the tmRNA/SRP/6S RNAs, in case the annotation misses one
    hard = genes.str.match(r"^(rrs|rrl|rrf)[A-Z]$") | genes.isin(["ssrA", "ssrS", "ffs", "rnpB"])
    drop = genes.isin(nc) | hard
    print("dropping %d rRNA/tRNA/ncRNA features" % int(drop.sum()))
    genes = genes[~drop]
counts = pd.concat([bulk.loc[genes], shx.loc[genes], casp.loc[genes]], axis=1)
counts = counts.round().astype(int)
print("%d shared genes (ERCC removed) x %d samples" % counts.shape)

# ------------------------------------------------- sample -> control pairing
# each value is the list of control columns a sample is divided by; more than one
# entry means the sample is referenced to their average
bulk_exp = [c for c in counts.columns if c.startswith("bulk|EXP")]
CONTROL = {}
for c in counts.columns:
    ds, s = c.split("|")
    if ds == "bulk":
        if s.startswith("EXP"):
            continue
        # each biological replicate is referenced to the exponential control of
        # that replicate; --bulk-control mean references the average instead
        CONTROL[c] = (["bulk|EXP_biorep" + s.split("biorep")[1]]
                      if args.bulk_control == "matched" else bulk_exp)
    elif ds == "shx":
        if s.endswith("t0"):
            continue
        CONTROL[c] = ["shx|EXP_biorep1t0"]
    else:
        if s.startswith("EXP"):
            continue
        CONTROL[c] = ["casp|EXP biorep2"]
if args.shx_late:
    # bulk_pca.py used only the last three SHX time points; same subset for comparability
    for c in [c for c in CONTROL if re.search(r"t[1-5]$", c)]:
        del CONTROL[c]

missing = {x for v in CONTROL.values() for x in v} - set(counts.columns)
if missing:
    raise SystemExit("control column(s) not found: %s" % sorted(missing))
print("%d treated samples paired to their own exponential control" % len(CONTROL))

# ------------------------------------------------------------- metadata
samples = list(CONTROL)
meta = pd.DataFrame(index=samples)
meta["dataset"] = [c.split("|")[0] for c in samples]
meta["sample"] = [c.split("|")[1] for c in samples]
meta["control"] = [" + ".join(CONTROL[c]) for c in samples]


def condition(name):
    ds, s = name.split("|")
    if ds == "bulk":
        return "CASP" if s.startswith("CASP") else "Disrupted"
    if ds == "shx":
        return "SHX_" + re.sub(r".*t", "t", s)
    return "CASP_time"


meta["condition"] = [condition(c) for c in samples]

# publication labels: <biology><experiment>, experiment 1 = bulk, 2 = time course
LABELS = {"bulk": {"Disrupted": "Dis-Arrest1", "CASP": "Reg-Arrest1"},
          "shx": "Dis-Arrest2",
          "casp": "Reg-Arrest2+SHX"}
meta["label"] = [LABELS["bulk"][c] if d == "bulk" else LABELS[d]
                 for d, c in zip(meta["dataset"], meta["condition"])]
meta["biogroup"] = np.where(meta["label"].str.startswith("Reg"),
                            "CASP-like", "Disrupted-like")
meta["batch"] = np.where(meta["dataset"] == "bulk", "bulk", "timecourse")

# --------------------------------------- normalisation + expression transforms
# the design is only there to satisfy the constructor: median-of-ratios size
# factors and the blind VST below are both design-independent
norm_meta = pd.DataFrame({"dataset": [c.split("|")[0] for c in counts.columns]},
                         index=counts.columns)
dds = DeseqDataSet(counts=counts.T, metadata=norm_meta, design_factors="dataset",
                   refit_cooks=False, quiet=True)
dds.fit_size_factors()
norm = pd.DataFrame(np.asarray(dds.X) / np.asarray(dds.obsm["size_factors"])[:, None],
                    index=counts.columns, columns=counts.index)

# expressed-gene filter: log ratios of near-zero counts are pure noise
keep = norm.mean(axis=0) >= args.min_count
print("%d/%d genes with mean normalised count >= %g"
      % (int(keep.sum()), len(keep), args.min_count))

dds.vst(use_design=False)  # blind VST: dispersion trend fit on the mean only
vst = pd.DataFrame(dds.layers["vst_counts"], index=counts.columns, columns=counts.index)
log2 = np.log2(norm + 1.0)

EXPR = {"log2": log2.loc[:, keep], "vst": vst.loc[:, keep]}


def fold_changes(expr):
    """Difference of each sample from its own exponential control, on `expr`'s scale.

    Where a sample maps to several controls the reference is their mean on this
    scale -- a geometric mean of expression, which is the natural centre for a
    ratio and is far less sensitive to one noisy control than an arithmetic one.
    """
    return pd.DataFrame({s: expr.loc[s] - expr.loc[CONTROL[s]].mean(axis=0)
                         for s in samples}).T


LFC = {k: fold_changes(v) for k, v in EXPR.items()}

# ------------------------------------------- is the VST actually needed here?
if args.compare_transforms:
    print("\n--- log2 vs VST diagnostics -------------------------------------")
    mean_norm = norm.loc[:, keep].mean(axis=0)
    lo = mean_norm <= mean_norm.quantile(0.25)
    hi = mean_norm >= mean_norm.quantile(0.75)
    for k, m in LFC.items():
        slo, shi = m.loc[:, lo].std().mean(), m.loc[:, hi].std().mean()
        print("%5s: across-sample SD of fold changes -- low-expression genes %.3f, "
              "high-expression genes %.3f, ratio %.2f" % (k, slo, shi, slo / shi))
    a, b = LFC["log2"].values.ravel(), LFC["vst"].values.ravel()
    print("gene-wise agreement between the two fold-change matrices: Pearson r = %.4f"
          % np.corrcoef(a, b)[0, 1])
    scores = {}
    for k, m in LFC.items():
        p = PCA(n_components=3)
        scores[k] = p.fit_transform(m.values - m.values.mean(0))
        print("%5s: PCA explained variance (%%) = %s"
              % (k, np.round(p.explained_variance_ratio_ * 100, 1)))
    for i in range(2):
        r = abs(np.corrcoef(scores["log2"][:, i], scores["vst"][:, i])[0, 1])
        print("  PC%d sample-score correlation between transforms: |r| = %.4f" % (i + 1, r))
    print("-----------------------------------------------------------------\n")

# ------------------------------------------------------- batch removal
mat = LFC[args.transform]
if args.remove_batch:
    # limma removeBatchEffect on the fold changes: fit  lfc ~ 1 + biogroup + batch
    # per gene and subtract only the fitted batch term.  `bulk` is the reference
    # level, so the time-course samples are shifted onto the bulk baseline.
    # Any residual offset here is what survived the exponential-control division,
    # i.e. a genuine batch response difference rather than a baseline difference.
    bio = pd.get_dummies(meta["biogroup"], drop_first=True).astype(float)
    batch = pd.get_dummies(meta["batch"], drop_first=False).astype(float)[["timecourse"]]
    keep_cols = np.column_stack([np.ones(len(meta)), bio.values])
    X = np.column_stack([keep_cols, batch.values])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise SystemExit("design is rank-deficient: batch and biogroup are confounded")
    coefs, *_ = np.linalg.lstsq(X, mat.values, rcond=None)
    batch_fit = batch.values @ coefs[keep_cols.shape[1]:]
    print("batch design: %d samples, rank %d of %d coefficients; "
          "mean |batch offset| = %.3f"
          % (len(meta), np.linalg.matrix_rank(X), X.shape[1],
             np.abs(coefs[keep_cols.shape[1]:]).mean()))
    mat = pd.DataFrame(mat.values - batch_fit, index=mat.index, columns=mat.columns)

# ------------------------------------------------------------------- PCA
if args.n_hvg:
    mat = mat[mat.var(axis=0).sort_values(ascending=False).index[:args.n_hvg]]
X = mat.values - mat.values.mean(axis=0)

pca = PCA(n_components=min(5, X.shape[0]))
pcs = pca.fit_transform(X)
ev = pca.explained_variance_ratio_ * 100

pc_df = pd.DataFrame(pcs[:, :3], index=mat.index, columns=["PC1", "PC2", "PC3"])
pc_df = pd.concat([meta, pc_df], axis=1)
pc_df.to_csv(os.path.join(OUT, "bulk_lfc_pca_coordinates%s.csv" % SUFFIX))
mat.T.to_csv(os.path.join(OUT, "bulk_log_fold_changes%s.csv" % SUFFIX))
# the size-factor-normalised counts behind those fold changes, for the same samples
# and the same expressed-gene set, so absolute-expression analyses use identical input.
# these counts precede the transform, the batch removal and the HVG selection, so the
# file name drops those parts of the suffix rather than claiming a correction it lacks
EXPR_SUFFIX = re.sub(r"_hvg\d+$", "", SUFFIX).replace("_nobatch", "")
norm.loc[samples, keep].T.to_csv(
    os.path.join(OUT, "bulk_normalized_counts%s.csv" % EXPR_SUFFIX))

# ----------------------------------------------------------------- plot
colors = {"Dis-Arrest1": "#c8102e", "Dis-Arrest2": "#c8102e",
          "Reg-Arrest1": "#1f5fa9", "Reg-Arrest2+SHX": "#1f5fa9"}
markers = {"bulk": "s", "timecourse": "^"}
ORDER = ["Reg-Arrest1", "Reg-Arrest2+SHX", "Dis-Arrest1", "Dis-Arrest2"]

fig, ax = plt.subplots(figsize=(7.6, 5))
for lab in ORDER:
    sub = pc_df[pc_df["label"] == lab]
    if sub.empty:
        continue
    ax.scatter(sub["PC1"], sub["PC2"], s=80, c=colors[lab],
               marker=markers[sub["batch"].iloc[0]], edgecolor="k", linewidth=0.5,
               label="%s  (batch %d)" % (lab, 1 if sub["batch"].iloc[0] == "bulk" else 2),
               zorder=3)


def tag(s):
    s = s.replace("Disrupted_biorep", "").replace("CASP_biorep", "")
    return s.replace("CASP biorep1 ", "").replace("min", "'")


for name, row in pc_df.iterrows():
    ax.annotate(tag(row["sample"]), (row["PC1"], row["PC2"]), fontsize=6,
                xytext=(5, 4), textcoords="offset points", color="0.35")

ax.set_xlabel("PC1 (%.1f%%)" % ev[0])
ax.set_ylabel("PC2 (%.1f%%)" % ev[1])
ax.set_title("Bulk RNA-seq PCA - fold change vs own exponential control (%s)" % args.transform
             + (", coding only" if args.drop_ncrna else "")
             + (", batch-corrected" if args.remove_batch else "")
             + (", top %d variable" % args.n_hvg if args.n_hvg else ""), fontsize=10)
ax.axhline(0, color="0.85", lw=0.7, zorder=0)
ax.axvline(0, color="0.85", lw=0.7, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
fig.tight_layout()
fig.savefig(os.path.join(OUT, "bulk_lfc_pca%s.svg" % SUFFIX))
fig.savefig(os.path.join(OUT, "bulk_lfc_pca%s.png" % SUFFIX), dpi=200)
print("explained variance (%):", np.round(ev, 2))
print("wrote", os.path.join(OUT, "bulk_lfc_pca%s.svg" % SUFFIX))
