"""Survival curves of the weighted-loading gene scores.

For each sample: y = number of genes with score > X, x = X, on log-log axes.
Colored by the sample's category in results/data_metrics/data_metrics.csv
(r = Reg-Arrest -> blue, d = Dis-Arrest -> red).

Excluded per request: adam_matrix_filtered.csv, deb_Ec_CDS_untreated.csv,
deb_KP_CDS_untreated.csv (adam_matrix_filtered2.csv is kept).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCORES = os.path.join(ROOT, "results", "eigenvector_analysis", "weighted_loading", "scores")
METRICS = os.path.join(ROOT, "results", "data_metrics", "data_metrics.csv")
OUT = os.path.join(ROOT, "results", "eigenvector_analysis", "figures")

EXCLUDE = {"adam_matrix_filtered.csv", "deb_Ec_CDS_untreated.csv", "deb_KP_CDS_untreated.csv"}
REG_COLOR = "#2166ac"   # blue  - Reg-Arrest
DIS_COLOR = "#b2182b"   # red   - Dis-Arrest

cat_map = pd.read_csv(METRICS, index_col=0).set_index("file_name")["category"].to_dict()

files = sorted(f for f in os.listdir(SCORES) if f.endswith(".csv") and f not in EXCLUDE)

fig, ax = plt.subplots(figsize=(5.2, 4.2))
rows = []
for fname in files:
    s = pd.read_csv(os.path.join(SCORES, fname))["score"].to_numpy(float)
    s = np.sort(s[s > 0])[::-1]          # descending, positive only (log x)
    counts = np.arange(1, s.size + 1)    # genes with score >= s[i]
    cat = cat_map.get(fname, "?")
    color = DIS_COLOR if cat == "d" else REG_COLOR
    ax.step(s, counts, where="post", color=color, lw=1.2, alpha=.75,
            label=fname[:-4])
    rows.append((fname, cat, s.size, float(s.max()), float(np.median(s))))

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("weighted-loading score, $X$")
ax.set_ylabel("number of genes with score $> X$")
ax.set_title("Gene-score survival curves")
ax.grid(True, which="major", lw=.4, alpha=.35)
handles = [plt.Line2D([], [], color=REG_COLOR, lw=1.6, label="Reg-Arrest (r)"),
           plt.Line2D([], [], color=DIS_COLOR, lw=1.6, label="Dis-Arrest (d)")]
ax.legend(handles=handles, frameon=False, loc="lower left")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "score_survival_loglog.svg"))
fig.savefig(os.path.join(OUT, "score_survival_loglog.png"), dpi=300)

print(pd.DataFrame(rows, columns=["file", "cat", "n_genes", "max_score", "median_score"]).to_string(index=False))
