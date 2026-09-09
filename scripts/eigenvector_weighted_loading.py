"""
Eigenvalue-weighted absolute-loading score per gene, for every scRNA-seq sample.

An alternative gene-level read-out to the coordination score in eigenvector_leverage.py.
For each sample we take the 10 leading modes of the gene-gene correlation matrix and
score every gene by

    score_i = sum_{k=1..10} |v_{k,i}| * lambda_k / P ,

where v_k is the k-th eigenvector (gene loadings), lambda_k its eigenvalue and P the
number of genes surviving the pipeline's all-zero filter. Differences from
coordination_score:

  * |v| rather than v**2 -- linear rather than quadratic in the loading, so moderately
    loaded genes count for more;
  * a fixed 10 leading modes rather than only the modes above lambda_max^scr, so samples
    with no above-threshold mode still get a score;
  * division by P, so samples with different panel sizes are on a comparable scale.

Caveat: |v| is not preserved under rotations within a degenerate eigen-subspace, so
unlike the coordination score this quantity is NOT rotation-invariant -- the same
stability caveats that apply to individual eigenvectors apply here. The per-mode
contributions are written out so that dependence can be inspected.

Run from the scripts/ directory:
    cd scripts
    python eigenvector_weighted_loading.py

Outputs -> results/eigenvector_analysis/weighted_loading/
  scores/<file>.csv   every kept gene, ranked, with its score and per-mode contributions
  summary.csv         one row per sample: P, eigenvalues, threshold, score statistics
  summary.txt         human-readable digest with the top genes per sample
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
N_MODES = 10       # leading modes summed into the score
N_REPORT = 25      # top genes listed per sample in summary.txt

DATA_DIR = os.path.join(ROOT, "data_for_paper")
# data_metrics.csv is the current 18-dataset metrics table (test8.csv is stale)
METRICS = os.path.join(ROOT, "results", "data_metrics", "data_metrics.csv")
OUT_DIR = os.path.join(ROOT, "results", "eigenvector_analysis", "weighted_loading")
os.makedirs(os.path.join(OUT_DIR, "scores"), exist_ok=True)


def clean_gene(col):
    """Strip the locus-tag prefix: 'LELOBEKK_araC' -> 'araC'; 'GFP' -> 'GFP'."""
    return col.split("_", 1)[1] if "_" in col else col


CAT_NAME = {"r": "Reg-Arrest", "d": "Dis-Arrest"}
cat_map, gmp_map = {}, {}
if os.path.exists(METRICS):
    metrics = pd.read_csv(METRICS, index_col=0)
    cat_map = dict(zip(metrics["file_name"], metrics["category"]))
    gmp_map = dict(zip(metrics["file_name"], metrics["sum_denoised_ev"]))

# ----------------------------------------------------------------------------
# main loop
# ----------------------------------------------------------------------------
files = sorted(f for f in os.listdir(DATA_DIR) if f.endswith(".csv"))
summary_rows = []
text_blocks = []

for fname in files:
    df = pd.read_csv(os.path.join(DATA_DIR, fname), index_col=0)
    genes = np.array(df.columns)
    m = df.values.astype(float)

    eigvals, eigvecs, threshold, kept_cols = af.get_eig_vectors(m, n_top=N_MODES)
    kept_genes = genes[np.asarray(kept_cols, dtype=bool)]
    P = eigvecs.shape[1]

    score = af.weighted_loading_score(eigvals, eigvecs, n_modes=N_MODES, p=P)

    # per-mode contributions, so the dependence on any single (possibly unstable)
    # eigenvector can be checked
    contrib = {f"mode{k+1}": np.abs(eigvecs[k]) * eigvals[k] / P
               for k in range(len(eigvals))}

    order = np.argsort(score)[::-1]
    out = pd.DataFrame({
        "gene": kept_genes[order],
        "gene_clean": [clean_gene(g) for g in kept_genes[order]],
        "score": score[order],
        "rank": np.arange(1, len(order) + 1),
    })
    for k, v in contrib.items():
        out[k] = v[order]
    out.to_csv(os.path.join(OUT_DIR, "scores", fname), index=False)

    cat = cat_map.get(fname, "?")
    gmp = gmp_map.get(fname, np.nan)
    # how concentrated is the score? fraction of total score in the top 5% of genes
    top5 = int(np.ceil(0.05 * P))
    conc = float(np.sort(score)[::-1][:top5].sum() / score.sum()) if score.sum() else np.nan

    summary_rows.append({
        "file_name": fname,
        "category": cat,
        "GMP_Cor": gmp,
        "n_genes_kept": P,
        "n_modes": len(eigvals),
        "threshold_scr": threshold,
        "n_modes_above_threshold": int(np.sum(eigvals > threshold)),
        "eig_1": eigvals[0],
        "eig_10": eigvals[-1],
        "score_sum": float(score.sum()),
        "score_mean": float(score.mean()),
        "score_max": float(score.max()),
        "score_median": float(np.median(score)),
        "score_top5pct_fraction": conc,
        "top_gene": clean_gene(kept_genes[order[0]]),
    })

    block = [
        f"\n=== {fname}  [{CAT_NAME.get(cat, cat)}]  GMP-Cor={gmp:.2f}  "
        f"P={P}  lambda_1={eigvals[0]:.2f}  lambda_max^scr={threshold:.2f} "
        f"({int(np.sum(eigvals > threshold))} of {len(eigvals)} modes above) ===",
        f"  score: max={score.max():.5f}  mean={score.mean():.5f}  "
        f"median={np.median(score):.5f}  top-5% of genes carry {conc:.1%} of the total",
        f"  top {N_REPORT}: " + ", ".join(
            f"{clean_gene(g)}({s:.4f})"
            for g, s in zip(kept_genes[order[:N_REPORT]], score[order[:N_REPORT]])),
    ]
    text_blocks.append("\n".join(block))
    print("\n".join(block))

# ----------------------------------------------------------------------------
# write summaries
# ----------------------------------------------------------------------------
summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(OUT_DIR, "summary.csv"), index=False)

header = (
    "Eigenvalue-weighted absolute-loading score per gene\n"
    "    score_i = sum_{k=1..%d} |v_{k,i}| * lambda_k / P\n"
    "P = genes surviving the all-zero filter; lambda_k from get_eig_vectors (same\n"
    "preprocessing as the rest of the pipeline: drop all-zero genes/cells, row-normalize,\n"
    "z-score columns). Per-mode contributions are in scores/<file>.csv.\n"
    "NOTE: |v| is not rotation-invariant within a degenerate subspace, so this score\n"
    "inherits the stability caveats of individual eigenvectors.\n" % N_MODES
)
with open(os.path.join(OUT_DIR, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(text_blocks))

print(f"\nWrote outputs to {OUT_DIR}")
