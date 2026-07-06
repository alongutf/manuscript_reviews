"""
GO enrichment of the leading eigenvectors at larger gene cutoffs (Reviewer #4, comment 1).

eigenvector_analysis.py runs GO on the top-50 loading genes of each mode. Because the
leading modes are strongly delocalized (the effective number of participating genes is
in the hundreds; see results/eigenvector_analysis/summary.csv), we re-run GO at larger
cutoffs -- the top 100 and top 200 genes (by |loading|) per mode -- to check whether a
coherent program emerges once more of each mode's participating genes are included.

Mirrors the preprocessing and GO machinery of eigenvector_analysis.py exactly; only the
gene cutoff changes. Eigenvalues / entropy are unchanged, so this script writes ONLY GO
results and a digest -- it does not touch summary.csv / summary.txt.

Run from the scripts/ directory:
    cd scripts
    python eigenvector_go_cutoffs.py

Outputs -> results/eigenvector_analysis/
  go_top100/<file>_mode<k>.csv   GO enrichment of each mode's top-100 loading genes
  go_top200/<file>_mode<k>.csv   GO enrichment of each mode's top-200 loading genes
  go_cutoffs_summary.txt         human-readable digest for both cutoffs
"""

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.getcwd())  # repo root when run from scripts/
sys.path.insert(0, ROOT)

import src.analysis_functions as af  # noqa: E402

# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
N_TOP = 5                 # number of leading eigenvectors to inspect per condition
GO_CUTOFFS = [100, 200]   # genes (by |loading|) fed to GO per mode

DATA_DIR = os.path.join(ROOT, "data_for_paper")
METRICS = os.path.join(ROOT, "results", "data_metrics", "test8.csv")
OUT_DIR = os.path.join(ROOT, "results", "eigenvector_analysis")
for c in GO_CUTOFFS:
    os.makedirs(os.path.join(OUT_DIR, f"go_top{c}"), exist_ok=True)


def clean_gene(col):
    """Strip the locus-tag prefix: 'LELOBEKK_araC' -> 'araC'; 'GFP' -> 'GFP'."""
    return col.split("_", 1)[1] if "_" in col else col


# condition labels + published GMP-Cor from the metrics table
metrics = pd.read_csv(METRICS, index_col=0)
cat_map = dict(zip(metrics["file_name"], metrics["category"]))
gmp_map = dict(zip(metrics["file_name"], metrics["sum_denoised_ev"]))
CAT_NAME = {"r": "regulated", "d": "dis-arrest"}

# ----------------------------------------------------------------------------
# GO machinery (reuse bulk_functions metadata + goatools)
# ----------------------------------------------------------------------------
import src.bulk_functions as bf  # noqa: E402
from goatools import obo_parser  # noqa: E402
from goatools.associations import read_gaf  # noqa: E402
from goatools.go_enrichment import GOEnrichmentStudy  # noqa: E402

_dag = obo_parser.GODag(bf.GO_OBO)
_assoc = read_gaf(bf.GAF_FILE)
_conv = bf.get_ID_conversion(bf.GTF_FILE)
print("GO machinery loaded.")


def names_to_ids(names):
    out = []
    for g in names:
        gid = _conv.get(clean_gene(g).lower())
        if gid is not None:
            out.append(gid)
    return out


def go_for_genes(study_genes, background_genes, out_path):
    """Run GO enrichment of study_genes against background_genes; write significant rows."""
    bg_ids = names_to_ids(background_genes)
    study_ids = names_to_ids(study_genes)
    if len(study_ids) < 3 or len(bg_ids) < 10:
        return None
    goea = GOEnrichmentStudy(
        bg_ids, _assoc, _dag,
        propagate_counts=False, alpha=0.05, methods=["fdr_bh"],
    )
    res = goea.run_study(study_ids, prt=None)
    sig = [r for r in res if r.enrichment == "e" and r.p_fdr_bh < 0.05]
    sig.sort(key=lambda r: r.p_fdr_bh)
    if not sig:
        return None
    df = pd.DataFrame({
        "GO_ID": [r.GO for r in sig],
        "Term": [r.name for r in sig],
        "Category": [r.NS for r in sig],
        "p_value": [r.p_uncorrected for r in sig],
        "FDR": [r.p_fdr_bh for r in sig],
        "Ratio_in_study": [r.ratio_in_study for r in sig],
        "Ratio_in_pop": [r.ratio_in_pop for r in sig],
    })
    df.to_csv(out_path, index=False)
    return df


# ----------------------------------------------------------------------------
# main loop: GO per mode at each cutoff
# ----------------------------------------------------------------------------
files = sorted(f for f in os.listdir(DATA_DIR) if f.endswith(".csv"))
text_blocks = []

for fname in files:
    df = pd.read_csv(os.path.join(DATA_DIR, fname), index_col=0)
    genes = np.array(df.columns)
    m = df.values.astype(float)

    eigvals, eigvecs, threshold, kept_cols = af.get_eig_vectors(m, n_top=N_TOP)
    kept_genes = genes[kept_cols]

    cat = cat_map.get(fname, "?")
    gmp = gmp_map.get(fname, np.nan)
    block = [
        f"\n=== {fname}  [{CAT_NAME.get(cat, cat)}]  GMP-Cor={gmp:.2f}  "
        f"lambda_max^scr={threshold:.3f}  (n_genes={len(kept_genes)}) ===",
    ]

    for k in range(len(eigvals)):
        v = eigvecs[k]
        order = np.argsort(np.abs(v))[::-1]
        above = "ABOVE" if eigvals[k] > threshold else "below"
        block.append(f"  mode {k+1}: eig={eigvals[k]:.3f} ({above} thr)")

        for cutoff in GO_CUTOFFS:
            n = min(cutoff, len(kept_genes))
            study = kept_genes[order[:n]]
            out_path = os.path.join(OUT_DIR, f"go_top{cutoff}", f"{fname[:-4]}_mode{k+1}.csv")
            go_df = go_for_genes(study, kept_genes, out_path)
            if go_df is not None and len(go_df):
                terms = "; ".join(go_df["Term"].head(4))
                block.append(f"      top{cutoff}: {len(go_df)} terms -> {terms}")
            else:
                block.append(f"      top{cutoff}: (no significant terms)")

    text_blocks.append("\n".join(block))
    print("\n".join(block))

# ----------------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------------
header = (
    "GO enrichment of leading eigenvectors at larger gene cutoffs "
    "(Reviewer #4, comment 1)\n"
    f"N_TOP={N_TOP}  GO_CUTOFFS={GO_CUTOFFS}\n"
    "GO is run on the top-N loading genes (by |loading|) of each mode, vs the gene-panel\n"
    "background. Compare with the top-50 results in results/eigenvector_analysis/go/.\n"
)
with open(os.path.join(OUT_DIR, "go_cutoffs_summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(text_blocks))

print(f"\nWrote GO outputs (top {GO_CUTOFFS}) to {OUT_DIR}")
