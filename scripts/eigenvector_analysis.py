"""
Eigenvector / gene-level interpretation of the GMP-Cor spectrum.

Response to Reviewer #4, comment 1: can we recover gene-level information
("what is still regulated") from the eigenvectors associated with the leading
eigenvalues, rather than only the scalar GMP-Cor?

For every cell x gene matrix in data_for_paper/ we:
  1. compute the top-N eigenvectors (gene loadings) of the gene-gene correlation
     matrix (get_eig_vectors), N small (default 5) because in some Dis-Arrest
     samples NO eigenvalue exceeds the scrambled threshold -- so we deliberately
     take the top-N modes regardless of the threshold and report, per mode,
     whether it sits above or below lambda_max^scr.
  2. extract the top-loading genes of each mode (by |loading|).
  3. compute the participation ratio of each mode (localized vs delocalized).
  4. run GO enrichment on the top-loading genes of each mode (vs the gene panel
     background) to test whether the leading modes recover coherent programs.

Run from the scripts/ directory:
    cd scripts
    python eigenvector_analysis.py

Outputs -> results/eigenvector_analysis/
  top_genes/<file>.csv        top-loading genes per mode (signed loading)
  summary.csv                 one row per (file, mode): eigenvalue, threshold,
                              above_threshold, participation_ratio, GMP-Cor
  go/<file>_mode<k>.csv       GO enrichment of mode k top-loading genes (if any)
  summary.txt                 human-readable digest
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
N_TOP = 5          # number of leading eigenvectors to inspect per condition
N_GENES = 30       # top-loading genes reported per mode
RUN_GO = True      # GO enrichment on top-loading genes of each mode
GO_N_GENES = 50    # genes (by |loading|) fed to GO per mode

DATA_DIR = os.path.join(ROOT, "data_for_paper")
METRICS = os.path.join(ROOT, "results", "data_metrics", "test8.csv")
OUT_DIR = os.path.join(ROOT, "results", "eigenvector_analysis")
os.makedirs(os.path.join(OUT_DIR, "top_genes"), exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, "go"), exist_ok=True)


def clean_gene(col):
    """Strip the locus-tag prefix: 'LELOBEKK_araC' -> 'araC'; 'GFP' -> 'GFP'."""
    return col.split("_", 1)[1] if "_" in col else col


# condition labels + published GMP-Cor from the metrics table
metrics = pd.read_csv(METRICS, index_col=0)
cat_map = dict(zip(metrics["file_name"], metrics["category"]))
gmp_map = dict(zip(metrics["file_name"], metrics["sum_denoised_ev"]))
CAT_NAME = {"r": "regulated", "d": "dis-arrest"}

# ----------------------------------------------------------------------------
# optional GO machinery (reuse bulk_functions metadata + goatools)
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
    conv = _go["conv"]
    out = []
    for g in names:
        gid = conv.get(clean_gene(g).lower())
        if gid is not None:
            out.append(gid)
    return out


def go_for_genes(study_genes, background_genes, out_path):
    """Run GO enrichment of study_genes against background_genes; write significant rows."""
    bg_ids = names_to_ids(background_genes)
    study_ids = names_to_ids(study_genes)
    if len(study_ids) < 3 or len(bg_ids) < 10:
        return None
    goea = _go["GOEnrichmentStudy"](
        bg_ids, _go["assoc"], _go["dag"],
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
# main loop
# ----------------------------------------------------------------------------
files = sorted(f for f in os.listdir(DATA_DIR) if f.endswith(".csv"))
summary_rows = []
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
        f"lambda_max^scr={threshold:.3f} ===",
    ]

    top_gene_records = {}
    for k in range(len(eigvals)):
        v = eigvecs[k]
        order = np.argsort(np.abs(v))[::-1]
        pr = af.participation_ratio(v)
        above = bool(eigvals[k] > threshold)

        top_idx = order[:N_GENES]
        top_gene_records[f"mode{k+1}_gene"] = kept_genes[top_idx]
        top_gene_records[f"mode{k+1}_loading"] = v[top_idx]

        summary_rows.append({
            "file_name": fname,
            "category": cat,
            "GMP_Cor": gmp,
            "mode": k + 1,
            "eigenvalue": eigvals[k],
            "threshold_scr": threshold,
            "above_threshold": above,
            "participation_ratio": pr,
            "n_genes_kept": len(kept_genes),
        })

        flag = "ABOVE" if above else "below"
        top_named = ", ".join(clean_gene(g) for g in kept_genes[order[:12]])
        block.append(
            f"  mode {k+1}: eig={eigvals[k]:.3f} ({flag} thr)  "
            f"PR={pr:.1f}/{len(kept_genes)}  top: {top_named}"
        )

        if RUN_GO:
            go_genes = kept_genes[order[:GO_N_GENES]]
            go_path = os.path.join(OUT_DIR, "go", f"{fname[:-4]}_mode{k+1}.csv")
            go_df = go_for_genes(go_genes, kept_genes, go_path)
            if go_df is not None and len(go_df):
                terms = "; ".join(go_df["Term"].head(3))
                block.append(f"        GO: {terms}")
            else:
                block.append("        GO: (no significant terms)")

    # per-file top-gene table
    maxlen = max(len(v) for v in top_gene_records.values())
    tg = pd.DataFrame({k: pd.Series(v) for k, v in top_gene_records.items()})
    tg.to_csv(os.path.join(OUT_DIR, "top_genes", fname), index=False)

    text_blocks.append("\n".join(block))
    print("\n".join(block))

# ----------------------------------------------------------------------------
# write summaries
# ----------------------------------------------------------------------------
summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(OUT_DIR, "summary.csv"), index=False)

header = (
    "Eigenvector / gene-level interpretation of GMP-Cor (Reviewer #4, comment 1)\n"
    f"N_TOP={N_TOP}  N_GENES={N_GENES}  GO={RUN_GO}\n"
    "Note: top-N modes are taken regardless of the scrambled threshold, because in\n"
    "some Dis-Arrest samples no eigenvalue exceeds it.\n"
)
with open(os.path.join(OUT_DIR, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(text_blocks))

print(f"\nWrote outputs to {OUT_DIR}")
