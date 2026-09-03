"""Intersect the RegulonDB sigma-38 (RpoS) sigmulon with the DESeq2 results.

Regulon list comes from metadata/regulondb_sigma38_regulon.txt, produced by
scripts/fetch_regulondb_sigma38.py (RegulonDB GraphQL API). Nothing here is
hand-curated.

Outputs:
    results/deseq_results/rpoS_regulon_hits.csv      per-gene, per-contrast
    results/deseq_results/rpoS_regulon_summary.csv   one row per contrast

Run from the repo root or from scripts/:
    python scripts/rpos_regulon_deseq.py
"""

import glob
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DESEQ = os.path.join(ROOT, "results", "deseq_results")
REGULON = os.path.join(ROOT, "metadata", "regulondb_sigma38_regulon.txt")
ALPHA = 0.05

# Subfolders of results/deseq_results holding DESeq2 contrast tables, in report order.
FOLDERS = ["from counts", "exp0224", "aggregated_sc"]

# Significance threshold per folder. aggregated_sc is aggregated single-cell data with
# roughly double the per-gene dispersion of the bulk sets (median lfcSE 0.50 vs 0.24) and
# a visibly noisier volcano, so it is held to a stricter padj cutoff. Any comparison of
# lobe sizes ACROSS folders has to account for this -- see the doc.
ALPHA_BY_FOLDER = {"aggregated_sc": 0.01}


def alpha_for(path):
    """Significance threshold to use for a given contrast file."""
    folder = os.path.basename(os.path.dirname(path))
    return ALPHA_BY_FOLDER.get(folder, ALPHA)


def contrast_files():
    """All DESeq2 result tables, grouped by folder in FOLDERS order."""
    out = []
    for folder in FOLDERS:
        out.extend(sorted(glob.glob(os.path.join(DESEQ, folder, "*.csv"))))
    return out


def load_regulon(path):
    """Parse the regulon gene list file.

    The file is one gene name per line, with a leading block of '#'-prefixed
    comment lines; any comment line of the form '# key: value' is captured as
    provenance metadata (e.g. RegulonDB release). Returns (gene set, header dict).
    """
    header = {}
    genes = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if ":" in line:
                    k, v = line.lstrip("# ").split(":", 1)
                    header[k.strip()] = v.strip()
                continue
            genes.append(line)
    return set(genes), header


def main():
    regulon, header = load_regulon(REGULON)
    print("RegulonDB release %s, %d sigma-38 genes"
          % (header.get("RegulonDB release", "?"), len(regulon)))

    files = contrast_files()

    # ---- per-contrast: intersect regulon with tested genes, test for shift ---
    rows, summary = [], []
    for path in files:
        tag = os.path.relpath(path, DESEQ).replace(os.sep, "/")
        alpha = alpha_for(path)
        d = pd.read_csv(path, index_col=0)
        d = d[d.padj.notna()]                      # DESeq2 independent filtering
        hit = d.index.isin(regulon)
        sub = d[hit].copy()
        sub.insert(0, "contrast", tag)
        sub.insert(1, "gene", sub.index)
        rows.append(sub)

        bg = d[~hit]
        sig = sub[sub.padj < alpha]
        # Is the regulon shifted relative to every other tested gene? A two-sided
        # Mann-Whitney U on log2FoldChange asks whether regulon genes' fold
        # changes are systematically higher or lower than the background's,
        # without assuming a particular direction or a normal distribution.
        if len(sub) and len(bg):
            u_p = stats.mannwhitneyu(sub.log2FoldChange, bg.log2FoldChange,
                                     alternative="two-sided").pvalue
        else:
            u_p = np.nan

        summary.append(dict(
            contrast=tag,
            alpha=alpha,
            genes_tested=len(d),
            regulon_tested=len(sub),
            sig=len(sig),
            up=int((sig.log2FoldChange > 0).sum()),
            down=int((sig.log2FoldChange < 0).sum()),
            median_l2fc_regulon=round(sub.log2FoldChange.median(), 3),
            median_l2fc_background=round(bg.log2FoldChange.median(), 3),
            mannwhitney_p=u_p,
        ))

    # ---- write outputs ---------------------------------------------------
    allrows = pd.concat(rows)
    cols = ["contrast", "gene", "baseMean", "log2FoldChange", "lfcSE",
            "pvalue", "padj"]
    allrows[cols].to_csv(os.path.join(DESEQ, "rpoS_regulon_hits.csv"), index=False)

    S = pd.DataFrame(summary)
    S.to_csv(os.path.join(DESEQ, "rpoS_regulon_summary.csv"), index=False)
    with pd.option_context("display.width", 200):
        print(S.to_string(index=False))

    # Which regulon genes are absent from the count matrices at all? Uses just the
    # first contrast file per folder, assuming every contrast within a folder was
    # run on the same gene set (true for a shared DESeq2 dataset per folder).
    for folder in FOLDERS:
        f = sorted(glob.glob(os.path.join(DESEQ, folder, "*.csv")))[0]
        idx = set(pd.read_csv(f, index_col=0).index)
        print("\n%s: %d/%d regulon genes present in the matrix"
              % (folder, len(regulon & idx), len(regulon)))


if __name__ == "__main__":
    main()
