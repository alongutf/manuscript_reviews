"""
PCA of bulk RNA-seq samples across three experiments.

Datasets / columns used
  1. bulk_data/bulk_count_data.csv          -> CASP_* and Disrupted_* columns
  2. bulk_data/time_in_shx_count_data.csv   -> t6, t7, t8 columns
  3. bulk_data/time_in_casp_count_data.csv  -> CASP columns

Pipeline: harmonise gene names -> intersect genes -> drop ERCC spike-ins ->
[optionally drop rRNA/tRNA/ncRNA] -> DESeq2 VST (pydeseq2) ->
top N most variable genes -> PCA.

Usage:
    python bulk_pca.py                 # all genes
    python bulk_pca.py --drop-ncrna    # coding genes only
"""

import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet

N_HVG = 1000
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "bulk_data")
OUT = os.path.join(ROOT, "results", "bulk_pca")
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument("--drop-ncrna", action="store_true",
                    help="remove rRNA, tRNA and other ncRNA features before VST")
parser.add_argument("--remove-batch", action="store_true",
                    help="subtract per-gene dataset offsets from the VST matrix, "
                         "protecting the CASP-like / Disrupted-like contrast")
args = parser.parse_args()
SUFFIX = ("_coding" if args.drop_ncrna else "") + ("_nobatch" if args.remove_batch else "")


def noncoding_ids(gtf=os.path.join(ROOT, "metadata", "genomic.gtf")):
    """All identifiers (symbol, gene_id, locus_tag) of non-coding features in the
    MG1655 reference annotation: rRNA, tRNA and ncRNA biotypes."""
    ids = set()
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("	")
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

# keep only the CASP/Disrupted arrest samples (drop the EXP_* exponential controls);
# unlike bulk_lfc_pca.py, technical replicates (biorep1a/1b/1c) are NOT collapsed here,
# so each one enters the PCA as its own column/sample
bulk = bulk[[c for c in bulk.columns if c.startswith(("CASP", "Disrupted"))]]
shx = shx[[c for c in shx.columns if c.rstrip().endswith(("t6", "t7", "t8"))]]
casp = casp[[c for c in casp.columns if c.startswith("CASP")]]

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
    print(f"dropping {int(drop.sum())} rRNA/tRNA/ncRNA features")
    genes = genes[~drop]
counts = pd.concat([bulk.loc[genes], shx.loc[genes], casp.loc[genes]], axis=1)
counts = counts.astype(int)
print(f"{counts.shape[0]} shared genes (ERCC removed) x {counts.shape[1]} samples")

# ------------------------------------------------------------- metadata
meta = pd.DataFrame(index=counts.columns)
meta["dataset"] = [c.split("|")[0] for c in counts.columns]
meta["sample"] = [c.split("|")[1] for c in counts.columns]


def condition(name):
    """Condition label for a "<dataset>|<sample>" column name: CASP/Disrupted for the
    bulk dataset, the SHX time point suffix for the shx dataset, or a single
    "CASP_time" label for every casp time-course sample.
    """
    ds, s = name.split("|")
    if ds == "bulk":
        return "CASP" if s.startswith("CASP") else "Disrupted"
    if ds == "shx":
        return "SHX_" + s[-2:]
    return "CASP_time"


meta["condition"] = [condition(c) for c in counts.columns]

# publication labels: <biology><experiment>, experiment 1 = bulk, 2 = time course
LABELS = {"bulk": {"Disrupted": "Dis-Arrest1", "CASP": "Reg-Arrest1"},
          "shx": "Dis-Arrest2",
          "casp": "Reg-Arrest2+SHX"}
meta["label"] = [LABELS["bulk"][c] if d == "bulk" else LABELS[d]
                 for d, c in zip(meta["dataset"], meta["condition"])]

# two-level biological factor shared across batches
meta["biogroup"] = np.where(meta["label"].str.startswith("Reg"),
                            "CASP-like", "Disrupted-like")
# sequencing batch (per the experimenter): bulk_count_data.csv is one experiment,
# the two time-course tables are a second one. Both batches contain both
# biogroups, so the design is a balanced 2x2 and the batch offsets are identified
# without absorbing the biological contrast.
meta["batch"] = np.where(meta["dataset"] == "bulk", "bulk", "timecourse")

# ------------------------------------------------------------------ VST
# variance-stabilising transform on absolute counts (not fold changes -- contrast
# with bulk_lfc_pca.py, which additionally normalises by size factor and subtracts
# each sample's own exponential control before PCA); no low-count gene filter is
# applied here, unlike bulk_lfc_pca.py's --min-count
dds = DeseqDataSet(counts=counts.T, metadata=meta, design_factors="condition",
                   refit_cooks=False, quiet=True)
dds.vst(use_design=False)  # blind VST: dispersion trend fit on the mean only
vst = pd.DataFrame(dds.layers["vst_counts"], index=counts.columns, columns=counts.index)

# ------------------------------------------------------- batch removal
if args.remove_batch:
    # limma removeBatchEffect: fit  vst ~ 1 + biogroup + batch  per gene,
    # then subtract only the fitted batch term. `bulk` is the reference level,
    # so the time-course samples are shifted onto the bulk baseline.
    bio = pd.get_dummies(meta["biogroup"], drop_first=True).astype(float)
    batch = pd.get_dummies(meta["batch"], drop_first=False).astype(float)
    batch = batch[["timecourse"]]                        # bulk = reference
    keep = np.column_stack([np.ones(len(meta)), bio.values])
    X = np.column_stack([keep, batch.values])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise SystemExit("design is rank-deficient: batch and biogroup are confounded")
    # rcond=None selects numpy's current default (machine-precision-based) cutoff for
    # treating small singular values as zero, rather than the deprecated fixed cutoff
    coefs, *_ = np.linalg.lstsq(X, vst.values, rcond=None)
    batch_fit = batch.values @ coefs[keep.shape[1]:]
    vst = pd.DataFrame(vst.values - batch_fit, index=vst.index, columns=vst.columns)
    print(f"removed batch offsets (rank {np.linalg.matrix_rank(X)} design, "
          f"{X.shape[1]} coefficients)")

# ------------------------------------------------------- HVG selection + PCA
# top N_HVG genes by variance across all samples/experiments; PCA on those (rather
# than all shared genes) so PC1/PC2 are driven by the genes that actually vary
hvg = vst.var(axis=0).sort_values(ascending=False).index[:N_HVG]
X = vst[hvg].values
X = X - X.mean(axis=0)     # explicit centring (redundant with, but harmless alongside,
                           # sklearn PCA's own internal centring)

pca = PCA(n_components=min(5, X.shape[0]))
pcs = pca.fit_transform(X)
ev = pca.explained_variance_ratio_ * 100

pc_df = pd.DataFrame(pcs[:, :3], index=counts.columns, columns=["PC1", "PC2", "PC3"])
pc_df = pd.concat([meta, pc_df], axis=1)
pc_df.to_csv(os.path.join(OUT, f"bulk_pca_coordinates{SUFFIX}.csv"))

# ----------------------------------------------------------------- plot
# colour = biology (red Dis-Arrest / blue Reg-Arrest), shape = batch
colors = {"Dis-Arrest1": "#c8102e", "Dis-Arrest2": "#c8102e",
          "Reg-Arrest1": "#1f5fa9", "Reg-Arrest2+SHX": "#1f5fa9"}
markers = {"bulk": "s", "timecourse": "^"}
ORDER = ["Reg-Arrest1", "Reg-Arrest2+SHX", "Dis-Arrest1", "Dis-Arrest2"]

fig, ax = plt.subplots(figsize=(7.6, 5))
for lab in ORDER:
    sub = pc_df[pc_df["label"] == lab]
    ax.scatter(sub["PC1"], sub["PC2"], s=80, c=colors[lab],
               marker=markers[sub["batch"].iloc[0]], edgecolor="k", linewidth=0.5,
               label=f"{lab}  (batch {1 if sub['batch'].iloc[0] == 'bulk' else 2})",
               zorder=3)
# short per-point tags: replicate id (bulk) or time point (time course)
def tag(s):
    """Compact point-annotation label: strip the dataset-specific sample prefix and
    abbreviate '...min' time points to "...'" so labels fit next to scatter points.
    """
    s = s.replace("Disrupted_biorep", "").replace("CASP_biorep", "")
    return s.replace("CASP biorep1 ", "").replace("min", "'")


for name, row in pc_df.iterrows():
    ax.annotate(tag(row["sample"]), (row["PC1"], row["PC2"]), fontsize=6,
                xytext=(5, 4), textcoords="offset points", color="0.35")

ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
ax.set_title(f"Bulk RNA-seq PCA — VST, top {N_HVG} variable genes"
             + (", coding only" if args.drop_ncrna else "")
             + (", batch-corrected" if args.remove_batch else ""), fontsize=11)
ax.axhline(0, color="0.85", lw=0.7, zorder=0)
ax.axvline(0, color="0.85", lw=0.7, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
fig.tight_layout()
fig.savefig(os.path.join(OUT, f"bulk_pca{SUFFIX}.svg"))
fig.savefig(os.path.join(OUT, f"bulk_pca{SUFFIX}.png"), dpi=200)
print("explained variance (%):", np.round(ev, 2))
print("wrote", os.path.join(OUT, f"bulk_pca{SUFFIX}.svg"))
