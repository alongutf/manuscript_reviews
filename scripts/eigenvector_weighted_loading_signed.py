"""
Sign-resolved eigenvalue-weighted loading score, per sample (Reviewer #4, comment 1).

eigenvector_weighted_loading.py scores every gene by the ABSOLUTE weighted loading

    score_i = sum_{k=1..10} |v_{k,i}| * lambda_k / P ,

and eigenvector_weighted_loading_go.py runs GO on the top 50. As with the per-mode
analysis (eigenvector_go_signed.py), taking |v| pools the two poles of every mode into
one ranking, so genes that anti-correlate with each other end up adjacent in the list.

This script drops the absolute value -- same formula, same 10 modes, same normalizer:

    signed_i = sum_{k=1..10} v_{k,i} * lambda_k / P ,

and splits each sample's genes by the sign of that score: the top N most POSITIVE and the
top N most NEGATIVE are enriched separately against the sample's own gene panel. The
existing absolute-value run at the same cutoff is the control.

IMPORTANT -- the sign problem. The overall sign of each eigenvector is arbitrary (SVD may
return +v_k or -v_k), so a signed SUM ACROSS MODES is not, by itself, a well-defined
quantity: flipping any one mode changes every gene's score. That is a real weakness of
this read-out relative to the per-mode analysis, where only the pole labels were arbitrary
and the partition was not. We make it reproducible by orienting each mode the same way as
eigenvector_go_signed.py -- the pole carrying the larger squared-loading mass is the
positive one -- which removes the dependence on LAPACK's arbitrary choice but does not
make the cross-mode sum canonical. As a robustness check the script also scores every
sample under the raw (unoriented) SVD signs and reports how many genes change sign and how
the two rankings correlate; read the GO results with that number in mind.

Run from the scripts/ directory:
    cd scripts
    python eigenvector_weighted_loading_signed.py

Outputs -> results/eigenvector_analysis/weighted_loading_signed/
  scores/<file>.csv   every kept gene with its signed score, the absolute-value score,
                      per-mode signed contributions and the raw-sign variant
  go_pos/<file>.csv   GO enrichment of the most positive genes (absent if none)
  go_neg/<file>.csv   GO enrichment of the most negative genes (absent if none)
  summary.csv         one row per (sample, cutoff): pole sizes, term counts, control
  summary.txt         human-readable digest
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
N_MODES = 10          # leading modes summed into the score (as in the absolute version)
CUTOFFS = [50, 100]   # genes taken from each pole; 50 matches the existing control run
ALPHA = 0.05

DATA_DIR = os.path.join(ROOT, "data_for_paper")
METRICS = os.path.join(ROOT, "results", "data_metrics", "data_metrics.csv")
BASE = os.path.join(ROOT, "results", "eigenvector_analysis")
OUT_DIR = os.path.join(BASE, "weighted_loading_signed")
# the absolute-value run this is compared against (top 50, written by
# eigenvector_weighted_loading_go.py)
ABS_GO_DIR = os.path.join(BASE, "weighted_loading", "go")
for sub in ("scores", "go_pos", "go_neg"):
    os.makedirs(os.path.join(OUT_DIR, sub), exist_ok=True)


def clean_gene(col):
    """Strip the locus-tag prefix: 'LELOBEKK_araC' -> 'araC'; 'GFP' -> 'GFP'."""
    return col.split("_", 1)[1] if "_" in col else col


CAT_NAME = {"r": "Reg-Arrest", "d": "Dis-Arrest"}
cat_map, gmp_map = {}, {}
if os.path.exists(METRICS):
    _m = pd.read_csv(METRICS, index_col=0)
    cat_map = dict(zip(_m["file_name"], _m["category"]))
    gmp_map = dict(zip(_m["file_name"], _m["sum_denoised_ev"]))

# ----------------------------------------------------------------------------
# GO machinery (identical to eigenvector_weighted_loading_go.py)
# ----------------------------------------------------------------------------
import src.bulk_functions as bf  # noqa: E402
from goatools import obo_parser  # noqa: E402
from goatools.associations import read_gaf  # noqa: E402
from goatools.go_enrichment import GOEnrichmentStudy  # noqa: E402

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


def abs_terms(fname):
    """Terms from the absolute-value weighted-loading run (top 50), or [] if none."""
    path = os.path.join(ABS_GO_DIR, fname)
    if not os.path.exists(path):
        return []
    return list(pd.read_csv(path)["Term"])


# ----------------------------------------------------------------------------
# main loop
# ----------------------------------------------------------------------------
files = sorted(f for f in os.listdir(DATA_DIR) if f.endswith(".csv"))
rows = []
blocks = []

for fname in files:
    df = pd.read_csv(os.path.join(DATA_DIR, fname), index_col=0)
    genes = np.array(df.columns)
    m = df.values.astype(float)

    eigvals, eigvecs, threshold, kept_cols = af.get_eig_vectors(m, n_top=N_MODES)
    kept_genes = genes[np.asarray(kept_cols, dtype=bool)]
    P = eigvecs.shape[1]

    # orient each mode so its heavier squared-mass pole is positive (see module docstring)
    oriented = np.empty_like(eigvecs)
    flipped = []
    for k in range(len(eigvals)):
        v = eigvecs[k]
        flip = (v[v > 0] ** 2).sum() < (v[v < 0] ** 2).sum()
        oriented[k] = -v if flip else v
        flipped.append(bool(flip))

    # same formula as weighted_loading_score, without the absolute value
    contrib = {f"mode{k+1}": oriented[k] * eigvals[k] / P for k in range(len(eigvals))}
    signed = np.sum([c for c in contrib.values()], axis=0)
    # the published absolute-value score, for side-by-side comparison in the table
    absolute = af.weighted_loading_score(eigvals, eigvecs, n_modes=N_MODES, p=P)
    # robustness check: the same signed score under the raw, unoriented SVD signs
    signed_raw = np.sum([eigvecs[k] * eigvals[k] / P for k in range(len(eigvals))], axis=0)
    sign_agree = float(np.mean(np.sign(signed) == np.sign(signed_raw)))
    rank_corr = float(pd.Series(signed).corr(pd.Series(signed_raw), method="spearman"))

    order = np.argsort(signed)[::-1]  # most positive first, most negative last
    out = pd.DataFrame({
        "gene": kept_genes[order],
        "gene_clean": [clean_gene(g) for g in kept_genes[order]],
        "signed_score": signed[order],
        "abs_score": absolute[order],
        "signed_score_rawsign": signed_raw[order],
        "rank_positive": np.arange(1, len(order) + 1),
    })
    for k, v in contrib.items():
        out[k] = v[order]
    out.to_csv(os.path.join(OUT_DIR, "scores", fname), index=False)

    cat = cat_map.get(fname, "?")
    gmp = gmp_map.get(fname, np.nan)
    n_pos, n_neg = int((signed > 0).sum()), int((signed < 0).sum())
    a_terms = abs_terms(fname)

    block = [
        f"\n=== {fname}  [{CAT_NAME.get(cat, cat)}]  GMP-Cor={gmp:.2f}  P={P} ===",
        f"  signed score: {n_pos} positive / {n_neg} negative genes  "
        f"(range {signed.min():.5f} .. {signed.max():.5f})",
        f"  orientation: {sum(flipped)} of {len(flipped)} modes flipped; vs raw SVD signs "
        f"{sign_agree:.0%} of genes keep their sign, Spearman r={rank_corr:.2f}",
        f"  |loading| control (top50): {len(a_terms)} terms"
        + (f" -> {'; '.join(a_terms[:4])}" if a_terms else ""),
    ]

    for cutoff in CUTOFFS:
        n = min(cutoff, n_pos, n_neg)
        pos_genes = out["gene_clean"].head(n).tolist()
        neg_genes = out["gene_clean"].tail(n).tolist()[::-1]  # most negative first
        background = out["gene_clean"].tolist()

        side_terms = {}
        for side, study, sub in (("pos", pos_genes, "go_pos"), ("neg", neg_genes, "go_neg")):
            tag = fname if cutoff == CUTOFFS[0] else f"{fname[:-4]}_top{cutoff}.csv"
            go_df, n_study, n_bg = go_for_genes(
                study, background, os.path.join(OUT_DIR, sub, tag))
            side_terms[side] = list(go_df["Term"]) if go_df is not None else []
            shown = "; ".join(side_terms[side][:4]) if side_terms[side] else \
                "(no significant terms)"
            block.append(f"    top{cutoff} {side}: {len(side_terms[side])} terms -> {shown}"
                         if side_terms[side] else
                         f"    top{cutoff} {side}: (no significant terms)")

        gained = sorted((set(side_terms["pos"]) | set(side_terms["neg"])) - set(a_terms))
        if gained and cutoff == 50:
            block.append(f"      NEW vs |loading| top50: {'; '.join(gained[:6])}")

        rows.append({
            "file_name": fname,
            "category": cat,
            "GMP_Cor": gmp,
            "cutoff": cutoff,
            "n_genes_kept": P,
            "n_pos_scores": n_pos,
            "n_neg_scores": n_neg,
            "n_per_pole": n,
            "modes_flipped": sum(flipped),
            "sign_agreement_vs_raw": sign_agree,
            "spearman_vs_raw": rank_corr,
            "n_terms_pos": len(side_terms["pos"]),
            "n_terms_neg": len(side_terms["neg"]),
            "n_terms_abs50": len(a_terms),
            "terms_pos": "; ".join(side_terms["pos"]),
            "terms_neg": "; ".join(side_terms["neg"]),
            "terms_abs50": "; ".join(a_terms),
            "new_vs_abs50": "; ".join(gained),
        })

    blocks.append("\n".join(block))
    print("\n".join(block))

# ----------------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------------
summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT_DIR, "summary.csv"), index=False)

c50 = summary[summary["cutoff"] == 50]
c100 = summary[summary["cutoff"] == 100]
header = (
    "Sign-resolved eigenvalue-weighted loading score (Reviewer #4, comment 1)\n"
    f"    signed_i = sum_(k=1..{N_MODES}) v_(k,i) * lambda_k / P   "
    "(same formula as the published score, without |.|)\n"
    "Each sample's most positive and most negative genes are enriched separately vs its\n"
    "own gene panel; the absolute-value top-50 run in weighted_loading/go/ is the control.\n"
    "Modes are oriented heavier-mass-pole-positive; the cross-mode signed sum is still not\n"
    "a canonical quantity -- see the per-sample sign-agreement figures below.\n\n"
    f"samples: {len(c50)}\n"
    f"  top50  -- positive pole enriched: {int((c50.n_terms_pos > 0).sum())}, "
    f"negative pole: {int((c50.n_terms_neg > 0).sum())}, "
    f"either: {int(((c50.n_terms_pos > 0) | (c50.n_terms_neg > 0)).sum())}\n"
    f"  top100 -- positive pole enriched: {int((c100.n_terms_pos > 0).sum())}, "
    f"negative pole: {int((c100.n_terms_neg > 0).sum())}, "
    f"either: {int(((c100.n_terms_pos > 0) | (c100.n_terms_neg > 0)).sum())}\n"
    f"  |loading| top50 control        : {int((c50.n_terms_abs50 > 0).sum())}\n"
    f"  median sign agreement vs raw SVD signs: {c50.sign_agreement_vs_raw.median():.0%}\n"
)
with open(os.path.join(OUT_DIR, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(blocks))

print("\n" + header)
print(f"Wrote outputs to {OUT_DIR}")
