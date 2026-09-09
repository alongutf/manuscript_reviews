"""
GO enrichment of the top-ranked genes from the eigenvalue-weighted loading score.

Reads the per-sample tables written by eigenvector_weighted_loading.py and runs the
goatools pipeline used elsewhere in the paper on the top N genes (by
score_i = sum_{k=1..10} |v_{k,i}| * lambda_k / P) against that sample's own gene panel
as background.

This is the direct counterpart of the coordination-score enrichment in
eigenvector_leverage.py, so the two gene-level read-outs can be compared on identical
samples, background and statistics.

Run from the scripts/ directory:
    cd scripts
    python eigenvector_weighted_loading_go.py

Outputs -> results/eigenvector_analysis/weighted_loading/go/
  <file>.csv     significant GO terms (FDR < 0.05) for that sample, or absent if none
  go_summary.txt human-readable digest across all samples
"""

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.getcwd())  # repo root when run from scripts/
sys.path.insert(0, ROOT)

import src.bulk_functions as bf  # noqa: E402
from goatools import obo_parser  # noqa: E402
from goatools.associations import read_gaf  # noqa: E402
from goatools.go_enrichment import GOEnrichmentStudy  # noqa: E402

# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
N_TOP = 50         # top-ranked genes fed to GO
ALPHA = 0.05       # FDR cutoff

BASE = os.path.join(ROOT, "results", "eigenvector_analysis", "weighted_loading")
SCORE_DIR = os.path.join(BASE, "scores")
OUT_DIR = os.path.join(BASE, "go")
os.makedirs(OUT_DIR, exist_ok=True)

METRICS = os.path.join(ROOT, "results", "data_metrics", "data_metrics.csv")
CAT_NAME = {"r": "Reg-Arrest", "d": "Dis-Arrest"}
cat_map, gmp_map = {}, {}
if os.path.exists(METRICS):
    _m = pd.read_csv(METRICS, index_col=0)
    cat_map = dict(zip(_m["file_name"], _m["category"]))
    gmp_map = dict(zip(_m["file_name"], _m["sum_denoised_ev"]))

print("Loading GO machinery...")
DAG = obo_parser.GODag(bf.GO_OBO)
ASSOC = read_gaf(bf.GAF_FILE)
CONV = bf.get_ID_conversion(bf.GTF_FILE)


def names_to_ids(names):
    """Map cleaned gene names to the GAF's gene ids, dropping anything unmapped."""
    out = []
    for g in names:
        gid = CONV.get(str(g).lower())
        if gid is not None:
            out.append(gid)
    return out


def go_for_genes(study_genes, background_genes, out_path):
    bg_ids = names_to_ids(background_genes)
    study_ids = names_to_ids(study_genes)
    if len(study_ids) < 3 or len(bg_ids) < 10:
        return None, len(study_ids), len(bg_ids)
    goea = GOEnrichmentStudy(bg_ids, ASSOC, DAG, propagate_counts=False,
                             alpha=ALPHA, methods=["fdr_bh"])
    res = goea.run_study(study_ids, prt=None)
    sig = [r for r in res if r.enrichment == "e" and r.p_fdr_bh < ALPHA]
    sig.sort(key=lambda r: r.p_fdr_bh)
    if not sig:
        return None, len(study_ids), len(bg_ids)
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
    return df, len(study_ids), len(bg_ids)


# ----------------------------------------------------------------------------
# main loop
# ----------------------------------------------------------------------------
files = sorted(f for f in os.listdir(SCORE_DIR) if f.endswith(".csv"))
blocks = []
rows = []

for fname in files:
    tbl = pd.read_csv(os.path.join(SCORE_DIR, fname))
    tbl = tbl.sort_values("rank")
    study = tbl["gene_clean"].head(N_TOP).tolist()
    background = tbl["gene_clean"].tolist()

    df, n_study, n_bg = go_for_genes(study, background,
                                     os.path.join(OUT_DIR, fname))
    cat = cat_map.get(fname, "?")
    gmp = gmp_map.get(fname, np.nan)
    n_terms = 0 if df is None else len(df)
    top_terms = "" if df is None else "; ".join(df["Term"].head(5))

    rows.append({
        "file_name": fname,
        "category": cat,
        "GMP_Cor": gmp,
        "n_top_genes": len(study),
        "n_study_mapped": n_study,
        "n_background_mapped": n_bg,
        "n_significant_terms": n_terms,
        "top_terms": top_terms,
    })

    head = (f"\n=== {fname}  [{CAT_NAME.get(cat, cat)}]  GMP-Cor={gmp:.2f}  "
            f"top{N_TOP} genes, {n_study} mapped to GO ids "
            f"(background {n_bg}) ===")
    if df is None:
        body = ("  no significant terms" if n_study >= 3
                else f"  SKIPPED: only {n_study} of {len(study)} genes map to GO ids "
                     f"(gene names not in the GTF -- e.g. numeric locus ids)")
    else:
        body = "\n".join(
            f"  {r.Term}  (FDR={r.FDR:.2e}, {r.Ratio_in_study} of study, "
            f"{r.Ratio_in_pop} of panel)"
            for r in df.head(10).itertuples())
    blocks.append(head + "\n" + body)
    print(head + "\n" + body)

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT_DIR, "go_summary.csv"), index=False)

header = (
    f"GO enrichment of the top {N_TOP} genes by eigenvalue-weighted loading score\n"
    f"    score_i = sum_(k=1..10) |v_(k,i)| * lambda_k / P\n"
    f"Background: that sample's own gene panel. FDR < {ALPHA} (fdr_bh), "
    f"enrichment only.\n"
    f"Scores from results/eigenvector_analysis/weighted_loading/scores/.\n"
)
with open(os.path.join(OUT_DIR, "go_summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(blocks))

n_hit = int((summary["n_significant_terms"] > 0).sum())
print(f"\n{n_hit} of {len(summary)} samples returned significant terms")
print(f"Wrote outputs to {OUT_DIR}")
