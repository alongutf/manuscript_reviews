"""
Signal-subspace gene-level interpretation of GMP-Cor (Reviewer #4, comment 1).

Rather than inspecting individual eigenvectors (which mix within the near-degenerate
signal subspace and are not cleanly interpretable -- see eigenvector_analysis.py),
we compute a per-gene COORDINATION SCORE: the diagonal of the de-noised signal
covariance reconstructed from all modes above the scrambled threshold,

    score_i = sum_{k: eig_k > lambda_max^scr} (eig_k - lambda_max^scr) * eigvec[k,i]^2 .

This is rotation-invariant within the signal subspace, so it is stable, and it
answers "which genes carry the surviving coordination" for each condition. We then:
  1. rank genes by coordination score per condition,
  2. run GO enrichment on the top-ranked genes vs the gene-panel background,
  3. correlate the score vectors across conditions ("compendium" seed) to test
     whether regulated conditions share a coordinated program that dis-arrest loses.

Run from scripts/:
    cd scripts
    python eigenvector_leverage.py

Outputs -> results/eigenvector_analysis/leverage/
  scores/<file>.csv            every kept gene with its coordination score (ranked)
  go/<file>.csv                GO enrichment of the top-ranked genes
  score_correlation.csv        cross-condition Spearman correlation of score vectors
  summary.txt                  human-readable digest
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = os.path.dirname(os.getcwd())  # repo root when run from scripts/
sys.path.insert(0, ROOT)

import src.analysis_functions as af  # noqa: E402

# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
TOP_FRAC = 0.05    # fraction of genes (by coordination score) fed to GO
TOP_MIN = 50       # but at least this many
RUN_GO = True

DATA_DIR = os.path.join(ROOT, "data_for_paper")
METRICS = os.path.join(ROOT, "results", "data_metrics", "test8.csv")
OUT_DIR = os.path.join(ROOT, "results", "eigenvector_analysis", "leverage")
os.makedirs(os.path.join(OUT_DIR, "scores"), exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, "go"), exist_ok=True)


def clean_gene(col):
    """Strip a leading 'probe-index_' prefix from a column/gene name, e.g.
    '3_dnaA' -> 'dnaA'; names with no underscore are returned unchanged.
    """
    return col.split("_", 1)[1] if "_" in col else col


# METRICS points at test8.csv -- see the log for this file: repo convention as of the
# most recent data_metrics update is that data_metrics.csv is current and testN.csv
# files (including test8.csv) are stale scramble realisations with different
# sum_denoised_ev values; this only affects the GMP-Cor value printed in the summary
# text below, not the coordination-score computation itself
metrics = pd.read_csv(METRICS, index_col=0)
cat_map = dict(zip(metrics["file_name"], metrics["category"]))
gmp_map = dict(zip(metrics["file_name"], metrics["sum_denoised_ev"]))
CAT_NAME = {"r": "regulated", "d": "dis-arrest"}

# ----------------------------------------------------------------------------
# optional GO machinery
# ----------------------------------------------------------------------------
_go = {}
if RUN_GO:
    try:
        import src.bulk_functions as bf
        from goatools import obo_parser
        from goatools.associations import read_gaf
        from goatools.go_enrichment import GOEnrichmentStudy

        _go["dag"] = obo_parser.GODag(bf.GO_OBO)
        _go["assoc"] = read_gaf(bf.GAF_FILE)
        _go["conv"] = bf.get_ID_conversion(bf.GTF_FILE)
        _go["GOEnrichmentStudy"] = GOEnrichmentStudy
        print("GO machinery loaded.")
    except Exception as e:  # noqa: BLE001
        print(f"GO disabled ({e}); skipping enrichment.")
        RUN_GO = False


def names_to_ids(names):
    """Map cleaned, lowercased gene names to GO-database IDs via the GTF-derived
    conversion table; names with no match are silently dropped.
    """
    conv = _go["conv"]
    out = []
    for g in names:
        gid = conv.get(clean_gene(g).lower())
        if gid is not None:
            out.append(gid)
    return out


def go_for_genes(study_genes, background_genes, out_path):
    """Run GO enrichment of `study_genes` (top-scoring genes) against
    `background_genes` (all genes kept in the signal subspace for that
    condition) and write the significant, enriched terms to `out_path`.
    Returns the result DataFrame, or None if there is too little data to
    test or nothing survives FDR correction.
    """
    bg_ids = names_to_ids(background_genes)
    study_ids = names_to_ids(study_genes)
    if len(study_ids) < 3 or len(bg_ids) < 10:
        return None
    goea = _go["GOEnrichmentStudy"](
        bg_ids, _go["assoc"], _go["dag"],
        propagate_counts=False, alpha=0.05, methods=["fdr_bh"],
    )
    res = goea.run_study(study_ids, prt=None)
    # keep only over-represented ("e" = enriched, not "p" = purified) terms below the
    # BH-FDR threshold
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
# main loop: per-condition coordination scores + GO
# ----------------------------------------------------------------------------
files = sorted(f for f in os.listdir(DATA_DIR) if f.endswith(".csv"))
score_series = {}   # file -> pd.Series(score indexed by cleaned gene name) for cross-corr
text_blocks = []

for fname in files:
    df = pd.read_csv(os.path.join(DATA_DIR, fname), index_col=0)
    genes = np.array(df.columns)
    m = df.values.astype(float)

    # threshold is the scrambled-null cutoff (lambda_max^scr); kept_cols are the genes
    # that survived af's own filtering, so `genes[kept_cols]` re-aligns names to columns
    eigvals, eigvecs, threshold, kept_cols = af.get_eig_vectors(m, n_top=None)
    kept_genes = genes[kept_cols]
    # per-gene coordination score: contribution of every above-threshold ("signal")
    # eigenmode to that gene's diagonal, weighted by how far above threshold each
    # mode's eigenvalue sits -- see the module docstring for the formula
    score = af.coordination_score(eigvals, eigvecs, threshold)
    n_signal = int(np.sum(eigvals > threshold))

    order = np.argsort(score)[::-1]
    ranked = pd.DataFrame({
        "gene": kept_genes[order],
        "gene_clean": [clean_gene(g) for g in kept_genes[order]],
        "coordination_score": score[order],
    })
    ranked.to_csv(os.path.join(OUT_DIR, "scores", fname), index=False)

    # store for cross-condition correlation (clean gene name -> score).
    # Collapse duplicate cleaned gene names (e.g. multiple probes) by taking the max.
    s_cond = pd.Series(score, index=[clean_gene(g) for g in kept_genes])
    score_series[fname] = s_cond.groupby(level=0).max()

    cat = cat_map.get(fname, "?")
    gmp = gmp_map.get(fname, np.nan)
    block = [
        f"\n=== {fname}  [{CAT_NAME.get(cat, cat)}]  GMP-Cor={gmp:.2f}  "
        f"signal modes={n_signal} ===",
        f"  top genes: {', '.join(ranked['gene_clean'].head(15))}",
    ]

    if RUN_GO and n_signal > 0:
        n_top = max(TOP_MIN, int(round(TOP_FRAC * len(kept_genes))))
        study = kept_genes[order[:n_top]]
        go_path = os.path.join(OUT_DIR, "go", fname)
        go_df = go_for_genes(study, kept_genes, go_path)
        if go_df is not None and len(go_df):
            block.append(f"  GO ({len(go_df)} terms): "
                         + "; ".join(go_df["Term"].head(6)))
        else:
            block.append("  GO: (no significant terms)")
    else:
        block.append("  GO: (no signal subspace -> skipped)")

    text_blocks.append("\n".join(block))
    print("\n".join(block))

# ----------------------------------------------------------------------------
# cross-condition correlation of coordination-score vectors ("compendium" seed)
# ----------------------------------------------------------------------------
all_genes = sorted(set().union(*[s.index for s in score_series.values()]))
mat = pd.DataFrame(index=all_genes)
for f, s in score_series.items():
    mat[f] = s.reindex(all_genes).fillna(0.0)

# Spearman correlation between conditions over the shared gene axis
rank_mat = mat.apply(lambda c: rankdata(c.values), axis=0)
corr = np.corrcoef(rank_mat.values.T)
corr_df = pd.DataFrame(corr, index=mat.columns, columns=mat.columns)
corr_df.to_csv(os.path.join(OUT_DIR, "score_correlation.csv"))

# regulated-vs-regulated, dis-vs-dis, regulated-vs-dis mean correlations
cats = np.array([cat_map.get(f, "?") for f in mat.columns])
def block_mean(a, b):
    """Mean pairwise Spearman correlation between every condition in category `a`
    and every condition in category `b` (excluding self-pairs when a == b).
    """
    vals = []
    cols = list(mat.columns)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            if (cats[i] == a and cats[j] == b) or (cats[i] == b and cats[j] == a):
                vals.append(corr[i, j])
    return float(np.mean(vals)) if vals else float("nan")

rr, dd, rd = block_mean("r", "r"), block_mean("d", "d"), block_mean("r", "d")

# ----------------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------------
header = (
    "Signal-subspace coordination score -- gene-level interpretation of GMP-Cor\n"
    "(Reviewer #4, comment 1)\n"
    f"TOP_FRAC={TOP_FRAC}  TOP_MIN={TOP_MIN}  GO={RUN_GO}\n\n"
    "Cross-condition Spearman correlation of coordination-score vectors:\n"
    f"  regulated vs regulated : {rr:.3f}\n"
    f"  dis-arrest vs dis-arrest: {dd:.3f}\n"
    f"  regulated vs dis-arrest : {rd:.3f}\n"
)
with open(os.path.join(OUT_DIR, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(text_blocks))

print("\n" + header)
print(f"Wrote outputs to {OUT_DIR}")
