"""
Sign-resolved GO enrichment of the leading eigenvectors (Reviewer #4, comment 1).

eigenvector_analysis.py / eigenvector_go_cutoffs.py rank the genes of each mode by
|loading|, so the two ends of an eigenvector -- the genes that move *together* and the
genes that move *opposite* to them -- are pooled into one study set. If a mode contrasts
two coherent programs (e.g. ribosomal genes against a stress regulon), pooling them can
wash the enrichment out: each half is diluted by the other, and a term enriched on one
side is tested against a study set twice as large.

This script splits each mode instead: the top 100 genes with POSITIVE loading and the
top 100 with NEGATIVE loading are enriched separately, against the same gene-panel
background. The union of the two sets is 200 genes of both signs, i.e. the |loading|
top-200 analysis already in results/eigenvector_analysis/go_top200/ -- so that run is
the natural control, with go_top100/ as the matched-study-size control.

Sign convention: the overall sign of an eigenvector is arbitrary (SVD may return +v or
-v). Each mode is oriented so that the pole carrying the larger squared-loading mass is
the POSITIVE one. "positive"/"negative" therefore means "dominant pole" / "opposing
pole", not a direction of expression change.

Preprocessing and GO machinery mirror eigenvector_go_cutoffs.py exactly.

Run from the scripts/ directory:
    cd scripts
    python eigenvector_go_signed.py

Outputs -> results/eigenvector_analysis/
  go_signed_pos/<file>_mode<k>.csv   GO enrichment of the top-100 positive-loading genes
  go_signed_neg/<file>_mode<k>.csv   GO enrichment of the top-100 negative-loading genes
  go_signed_summary.txt              human-readable digest
  go_signed_comparison.csv           per (file, mode): term counts for pos / neg /
                                     |loading| top-100 / top-200, plus the sign-split
                                     statistics of the mode
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
N_TOP = 5        # number of leading eigenvectors to inspect per condition
N_SIDE = 100     # genes taken from each pole (positive / negative) of a mode

DATA_DIR = os.path.join(ROOT, "data_for_paper")
# same (stale) metrics table the two sibling scripts read; used only for the printed
# header labels -- see the note in eigenvector_analysis.py
METRICS = os.path.join(ROOT, "results", "data_metrics", "test8.csv")
OUT_DIR = os.path.join(ROOT, "results", "eigenvector_analysis")
POS_DIR = os.path.join(OUT_DIR, "go_signed_pos")
NEG_DIR = os.path.join(OUT_DIR, "go_signed_neg")
os.makedirs(POS_DIR, exist_ok=True)
os.makedirs(NEG_DIR, exist_ok=True)
# the |loading| runs this is compared against (written by eigenvector_go_cutoffs.py)
ABS_DIRS = {100: os.path.join(OUT_DIR, "go_top100"),
            200: os.path.join(OUT_DIR, "go_top200")}


def clean_gene(col):
    """Strip the locus-tag prefix: 'LELOBEKK_araC' -> 'araC'; 'GFP' -> 'GFP'."""
    return col.split("_", 1)[1] if "_" in col else col


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
    """Map gene-panel column names to NCBI gene IDs, dropping names absent from the
    GTF-derived conversion table (silently, e.g. reporter / spike-in probes)."""
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
    # too few mapped genes for GO enrichment to be meaningful (or valid) -- skip
    if len(study_ids) < 3 or len(bg_ids) < 10:
        return None
    goea = GOEnrichmentStudy(
        bg_ids, _assoc, _dag,
        propagate_counts=False, alpha=0.05, methods=["fdr_bh"],
    )
    res = goea.run_study(study_ids, prt=None)
    # keep only enriched ("e", vs purged "p") terms passing the 5% FDR threshold
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


def abs_terms(fname, mode, cutoff):
    """Terms found by the existing |loading| run at `cutoff`, or [] if that run found
    nothing (eigenvector_go_cutoffs.py writes a file only when something is significant)."""
    path = os.path.join(ABS_DIRS[cutoff], f"{fname[:-4]}_mode{mode}.csv")
    if not os.path.exists(path):
        return []
    return list(pd.read_csv(path)["Term"])


def write_contrast_table(comp, out_path, n_show=4):
    """Per sample and mode, what the two poles contrast: the positive-pole terms
    against the negative-pole terms. Modes enriched on BOTH poles -- a genuine
    two-program contrast rather than one program against an unstructured bulk --
    are flagged CONTRAST."""
    def fmt(terms):
        if not terms:
            return "(none)"
        parts = terms.split("; ")
        out = "; ".join(parts[:n_show])
        return out + (f"  (+{len(parts) - n_show} more)" if len(parts) > n_show else "")

    lines = [
        "Per-sample / per-mode pole contrasts (see documents/eigenvector_go_signed.md)",
        "POS = dominant (heavier squared-mass) pole, NEG = opposing pole; the labels are",
        "an orientation convention, not a direction of expression change.",
        "CONTRAST = both poles significantly enriched.",
        "",
    ]
    for fname, g in comp.groupby("file_name", sort=True):
        lines.append(f"=== {fname[:-4]}  [{CAT_NAME.get(g['category'].iloc[0], '?')}] ===")
        for _, r in g.iterrows():
            both = r["n_terms_pos100"] > 0 and r["n_terms_neg100"] > 0
            lines.append(
                "  mode %d  eig=%.2f (%s thr)  mass_pos=%.2f%s" % (
                    r["mode"], r["eigenvalue"],
                    "ABOVE" if r["above_threshold"] else "below",
                    r["mass_pos"], "   ** CONTRAST **" if both else "")
            )
            lines.append("     POS(%d): %s" % (r["n_terms_pos100"], fmt(r["terms_pos100"])))
            lines.append("     NEG(%d): %s" % (r["n_terms_neg100"], fmt(r["terms_neg100"])))
        lines.append("")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


# ----------------------------------------------------------------------------
# main loop: GO per mode, positive and negative pole separately
# ----------------------------------------------------------------------------
files = sorted(f for f in os.listdir(DATA_DIR) if f.endswith(".csv"))
text_blocks = []
rows = []

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
        # orient the mode so the heavier pole is the positive one (SVD sign is arbitrary)
        if (v[v > 0] ** 2).sum() < (v[v < 0] ** 2).sum():
            v = -v
        pos_idx = np.where(v > 0)[0]
        neg_idx = np.where(v < 0)[0]
        # within each pole, rank by |loading| and keep the strongest N_SIDE
        pos_idx = pos_idx[np.argsort(v[pos_idx])[::-1]][:N_SIDE]
        neg_idx = neg_idx[np.argsort(v[neg_idx])][:N_SIDE]

        pos_mass = float((v[v > 0] ** 2).sum())  # squared mass per pole (poles sum to 1)
        above = "ABOVE" if eigvals[k] > threshold else "below"
        block.append(
            f"  mode {k+1}: eig={eigvals[k]:.3f} ({above} thr)  "
            f"n_pos={int((v > 0).sum())} n_neg={int((v < 0).sum())}  "
            f"mass_pos={pos_mass:.2f}"
        )

        side_terms = {}
        for side, idx, out_dir in (("pos", pos_idx, POS_DIR), ("neg", neg_idx, NEG_DIR)):
            study = kept_genes[idx]
            out_path = os.path.join(out_dir, f"{fname[:-4]}_mode{k+1}.csv")
            go_df = go_for_genes(study, kept_genes, out_path)
            side_terms[side] = list(go_df["Term"]) if go_df is not None else []
            if side_terms[side]:
                shown = "; ".join(side_terms[side][:4])
                block.append(f"      {side}100 (n={len(study)}): "
                             f"{len(side_terms[side])} terms -> {shown}")
            else:
                block.append(f"      {side}100 (n={len(study)}): (no significant terms)")

        a100, a200 = abs_terms(fname, k + 1, 100), abs_terms(fname, k + 1, 200)
        block.append(f"      |loading| control: top100={len(a100)} terms, "
                     f"top200={len(a200)} terms")
        # terms the sign split recovers that the pooled |loading| top-200 did not
        gained = sorted((set(side_terms["pos"]) | set(side_terms["neg"])) - set(a200))
        if gained:
            block.append(f"      NEW vs |loading| top200: {'; '.join(gained[:6])}")

        rows.append({
            "file_name": fname,
            "category": cat,
            "mode": k + 1,
            "eigenvalue": eigvals[k],
            "threshold_scr": threshold,
            "above_threshold": eigvals[k] > threshold,
            "n_genes_kept": len(kept_genes),
            "n_pos_loadings": int((v > 0).sum()),
            "n_neg_loadings": int((v < 0).sum()),
            "mass_pos": pos_mass,
            "n_terms_pos100": len(side_terms["pos"]),
            "n_terms_neg100": len(side_terms["neg"]),
            "n_terms_abs100": len(a100),
            "n_terms_abs200": len(a200),
            "terms_pos100": "; ".join(side_terms["pos"]),
            "terms_neg100": "; ".join(side_terms["neg"]),
            "terms_abs200": "; ".join(a200),
            "new_vs_abs200": "; ".join(gained),
        })

    text_blocks.append("\n".join(block))
    print("\n".join(block))

# ----------------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------------
comp = pd.DataFrame(rows)
comp.to_csv(os.path.join(OUT_DIR, "go_signed_comparison.csv"), index=False)

n_pos = int((comp["n_terms_pos100"] > 0).sum())
n_neg = int((comp["n_terms_neg100"] > 0).sum())
n_any_signed = int(((comp["n_terms_pos100"] > 0) | (comp["n_terms_neg100"] > 0)).sum())
n_abs100 = int((comp["n_terms_abs100"] > 0).sum())
n_abs200 = int((comp["n_terms_abs200"] > 0).sum())
n_new = int((comp["new_vs_abs200"] != "").sum())

header = (
    "Sign-resolved GO enrichment of leading eigenvectors (Reviewer #4, comment 1)\n"
    f"N_TOP={N_TOP}  N_SIDE={N_SIDE} genes per pole\n"
    "Each mode's top-100 POSITIVE and top-100 NEGATIVE loading genes are enriched\n"
    "separately vs the gene-panel background. The pooled |loading| top-100 / top-200 runs\n"
    "in go_top100/ and go_top200/ are the controls (top-200 = the same 200 genes, pooled).\n"
    "Mode sign is oriented so the heavier squared-mass pole is the 'positive' one.\n\n"
    f"modes tested: {len(comp)}\n"
    f"  with >=1 term on the positive pole  : {n_pos}\n"
    f"  with >=1 term on the negative pole  : {n_neg}\n"
    f"  with >=1 term on either pole        : {n_any_signed}\n"
    f"  |loading| top-100 control           : {n_abs100}\n"
    f"  |loading| top-200 control           : {n_abs200}\n"
    f"  modes gaining a term not in top-200 : {n_new}\n"
)
with open(os.path.join(OUT_DIR, "go_signed_summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(text_blocks))

write_contrast_table(comp, os.path.join(OUT_DIR, "go_signed_contrast_table.txt"))

print("\n" + header)
print(f"Wrote sign-resolved GO outputs to {OUT_DIR}")
