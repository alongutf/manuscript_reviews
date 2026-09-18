"""
Spectral (eigenvalue-weighted) effective gene number per sample.

eigenvector_analysis.py reports an entropy PER MODE, from the probability vector
p_i = v_{k,i}^2 of a single eigenvector. This script instead builds ONE probability
vector per sample by pooling modes with eigenvalue weights:

    P_i = sum_k lambda_k * v_{i,k}^2 / P

and reports H = -sum_i P_i ln P_i, the normalized entropy H/ln(P), and the effective
gene number exp(H).

VALIDATION -- and why the sum must be truncated. Summing over ALL modes makes this
quantity degenerate. Since C = sum_k lambda_k v_k v_k^T, the inner sum is exactly the
diagonal of the correlation matrix, C_ii = 1 for every gene, so

    P_i = 1/P  for all i,   sum_i P_i = 1,   H = ln(P),   exp(H) = P

identically, for every sample. The normalization is correct (it does sum to 1) but the
vector is uniform by construction and carries no information. The script verifies this
numerically and reports it as a check.

The quantity becomes informative once the sum is restricted to the SIGNAL modes, which
is what this script reports. Truncated at K modes the numerator

    s_i^(K) = sum_{k<=K} lambda_k * v_{i,k}^2

is the diagonal of the de-noised correlation matrix -- the same quantity as
af.coordination_score -- renormalized here to a probability vector, s_i / sum_j s_j, so
its entropy is comparable across samples with different panel sizes. Three truncations
are reported:

  above  : modes with lambda_k > lambda_max^scr (the GMP-Cor signal modes); samples with
           no such mode fall back to K=1 and are flagged
  top5   : the 5 modes inspected by eigenvector_analysis.py
  top10  : the 10 modes summed by eigenvector_weighted_loading.py

Run from the scripts/ directory:
    cd scripts
    python eigenvector_spectral_entropy.py

Outputs -> results/eigenvector_analysis/spectral_entropy/
  per_gene/<file>.csv   every kept gene with its P_i at each truncation, ranked
  summary.csv           one row per (sample, truncation): H, norm_H, exp(H), checks
  summary.txt           human-readable digest
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
TOP_KS = [5, 10]   # fixed-K truncations reported alongside the above-threshold one
N_REPORT = 15      # top genes listed per sample in summary.txt

DATA_DIR = os.path.join(ROOT, "data_for_paper")
METRICS = os.path.join(ROOT, "results", "data_metrics", "data_metrics.csv")
OUT_DIR = os.path.join(ROOT, "results", "eigenvector_analysis", "spectral_entropy")
os.makedirs(os.path.join(OUT_DIR, "per_gene"), exist_ok=True)


def clean_gene(col):
    """Strip the locus-tag prefix: 'LELOBEKK_araC' -> 'araC'; 'GFP' -> 'GFP'."""
    return col.split("_", 1)[1] if "_" in col else col


CAT_NAME = {"r": "Reg-Arrest", "d": "Dis-Arrest"}
cat_map, gmp_map = {}, {}
if os.path.exists(METRICS):
    _m = pd.read_csv(METRICS, index_col=0)
    cat_map = dict(zip(_m["file_name"], _m["category"]))
    gmp_map = dict(zip(_m["file_name"], _m["sum_denoised_ev"]))


def entropy_of(p):
    """Shannon entropy (nats) of a probability vector, ignoring zero entries."""
    p = np.asarray(p, dtype=float)
    nz = p[p > 0]
    return float(-(nz * np.log(nz)).sum())


def weighted_probability(eigvals, eigvecs, k, p_total):
    """P_i from the k leading modes: sum_{j<k} lambda_j v_{i,j}^2, renormalized to 1.

    Untruncated (k = all modes) the numerator is the unit diagonal of the correlation
    matrix and p_total = P, so this returns the uniform vector -- see the module
    docstring. Truncated it is the de-noised diagonal (af.coordination_score)."""
    s = (eigvals[:k, None] * eigvecs[:k] ** 2).sum(axis=0)
    total = s.sum()
    return s / total if total > 0 else np.full(len(s), 1.0 / p_total)


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

    # all modes, so the degeneracy check below is exact
    eigvals, eigvecs, threshold, kept_cols = af.get_eig_vectors(m, n_top=None)
    kept_genes = genes[np.asarray(kept_cols, dtype=bool)]
    P = eigvecs.shape[1]

    # --- validation: the untruncated vector must be uniform and sum to 1 ---------
    p_all = (eigvals[:, None] * eigvecs ** 2).sum(axis=0) / P
    sum_all = float(p_all.sum())
    uniform_dev = float(np.max(np.abs(p_all - 1.0 / P)) * P)  # relative deviation
    H_all = entropy_of(p_all)

    n_above = int(np.sum(eigvals > threshold))
    truncations = [("above", max(n_above, 1))] + [(f"top{k}", min(k, len(eigvals)))
                                                  for k in TOP_KS]

    cat = cat_map.get(fname, "?")
    gmp = gmp_map.get(fname, np.nan)
    block = [
        f"\n=== {fname}  [{CAT_NAME.get(cat, cat)}]  GMP-Cor={gmp:.2f}  P={P} ===",
        f"  check (all {len(eigvals)} modes): sum(P_i)={sum_all:.12f}  "
        f"max|P_i - 1/P|*P={uniform_dev:.2e}  H={H_all:.4f}  ln(P)={np.log(P):.4f}  "
        f"exp(H)={np.exp(H_all):.1f}  -> uniform by construction",
        f"  lambda_max^scr={threshold:.3f}, {n_above} modes above"
        + ("  (none above -- 'above' falls back to K=1)" if n_above == 0 else ""),
    ]

    per_gene = {"gene": kept_genes,
                "gene_clean": [clean_gene(g) for g in kept_genes]}

    for label, k in truncations:
        p = weighted_probability(eigvals, eigvecs, k, P)
        H = entropy_of(p)
        per_gene[f"P_i_{label}"] = p

        rows.append({
            "file_name": fname,
            "category": cat,
            "GMP_Cor": gmp,
            "truncation": label,
            "n_modes_used": k,
            "n_genes_kept": P,
            "threshold_scr": threshold,
            "n_modes_above_threshold": n_above,
            "entropy": H,
            "norm_entropy": H / np.log(P),
            "effective_n_genes": float(np.exp(H)),
            "eff_fraction_of_P": float(np.exp(H)) / P,
            "top_gene": clean_gene(kept_genes[int(np.argmax(p))]),
            # degeneracy check, identical for every truncation of a given sample
            "check_sum_all_modes": sum_all,
            "check_uniform_dev_all_modes": uniform_dev,
            "check_eff_n_all_modes": float(np.exp(H_all)),
        })

        top = np.argsort(p)[::-1][:N_REPORT]
        block.append(
            f"  {label:>5s} (K={k:2d}): H={H:.3f}  norm={H / np.log(P):.3f}  "
            f"eff_genes={np.exp(H):.1f} ({np.exp(H) / P:.1%} of P)")
        block.append("           top: " + ", ".join(
            f"{clean_gene(kept_genes[i])}({p[i] * P:.1f}x)" for i in top))

    out = pd.DataFrame(per_gene)
    out = out.sort_values(f"P_i_{truncations[0][0]}", ascending=False)
    out.to_csv(os.path.join(OUT_DIR, "per_gene", fname), index=False)

    blocks.append("\n".join(block))
    print("\n".join(block))

# ----------------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------------
summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT_DIR, "summary.csv"), index=False)

chk = summary.drop_duplicates("file_name")
header = (
    "Spectral (eigenvalue-weighted) effective gene number\n"
    "    P_i = sum_k lambda_k * v_(i,k)^2 / P ,  H = -sum_i P_i ln P_i ,  eff = exp(H)\n\n"
    "VALIDATION over all modes: sum(P_i) = 1 exactly, but P_i = 1/P uniformly, because\n"
    "sum_k lambda_k v_(i,k)^2 is the unit diagonal of the correlation matrix. So the\n"
    "untruncated entropy is ln(P) and the effective gene number is P for every sample --\n"
    "no information. Figures below are therefore from truncated sums (signal modes only),\n"
    "where the numerator is the de-noised diagonal (= af.coordination_score).\n\n"
    f"  samples: {len(chk)}\n"
    f"  sum(P_i) over all modes: min {chk.check_sum_all_modes.min():.12f}, "
    f"max {chk.check_sum_all_modes.max():.12f}\n"
    f"  max relative deviation from uniform: "
    f"{chk.check_uniform_dev_all_modes.max():.2e}\n"
    f"  exp(H) over all modes equals P in {int((abs(chk.check_eff_n_all_modes - chk.n_genes_kept) < 1e-6).sum())}"
    f" of {len(chk)} samples\n"
)
for label in ["above"] + [f"top{k}" for k in TOP_KS]:
    sub = summary[summary["truncation"] == label]
    header += (f"  {label:>5s}: eff_genes median {sub.effective_n_genes.median():.0f} "
               f"({sub.eff_fraction_of_P.median():.1%} of P), "
               f"range {sub.effective_n_genes.min():.0f}-{sub.effective_n_genes.max():.0f}\n")

with open(os.path.join(OUT_DIR, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(blocks))

print("\n" + header)
print(f"Wrote outputs to {OUT_DIR}")
