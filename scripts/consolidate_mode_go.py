"""
Consolidate the per-mode GO enrichment results (produced by eigenvector_analysis.py,
modes 1-5) into a single table per condition, annotating in which mode(s) each GO
term was enriched and the FDR in each.

Reads:  results/eigenvector_analysis/go/<condition>_mode<k>.csv
Writes: results/eigenvector_analysis/mode_go/<condition>.csv   (one row per term)
        results/eigenvector_analysis/mode_go/all_conditions.csv (long format)
        results/eigenvector_analysis/mode_go/summary.txt

Run from scripts/:
    cd scripts
    python consolidate_mode_go.py
"""

import os
import re

import pandas as pd

ROOT = os.path.dirname(os.getcwd())
GO_DIR = os.path.join(ROOT, "results", "eigenvector_analysis", "go")
OUT_DIR = os.path.join(ROOT, "results", "eigenvector_analysis", "mode_go")
METRICS = os.path.join(ROOT, "results", "data_metrics", "test8.csv")
os.makedirs(OUT_DIR, exist_ok=True)

metrics = pd.read_csv(METRICS, index_col=0)
cat_map = dict(zip(metrics["file_name"], metrics["category"]))
gmp_map = dict(zip(metrics["file_name"], metrics["sum_denoised_ev"]))
CAT_NAME = {"r": "regulated", "d": "dis-arrest"}

# gather per-mode GO files -> {condition_csv: {mode: dataframe}}
pat = re.compile(r"^(?P<cond>.+)_mode(?P<mode>\d+)\.csv$")
by_cond = {}
for fn in sorted(os.listdir(GO_DIR)):
    mobj = pat.match(fn)
    if not mobj:
        continue
    cond = mobj.group("cond") + ".csv"          # original data file name
    mode = int(mobj.group("mode"))
    df = pd.read_csv(os.path.join(GO_DIR, fn))
    by_cond.setdefault(cond, {})[mode] = df

# every condition present in the metrics table (so empties are reported too)
all_conditions = list(metrics["file_name"])

long_rows = []
text_blocks = []
for cond in all_conditions:
    cat = cat_map.get(cond, "?")
    gmp = gmp_map.get(cond, float("nan"))
    modes = by_cond.get(cond, {})

    # term -> {Term, Category, modes={mode: fdr}}
    term_info = {}
    for mode, df in sorted(modes.items()):
        for _, r in df.iterrows():
            info = term_info.setdefault(
                r["GO_ID"], {"Term": r["Term"], "Category": r["Category"], "modes": {}})
            info["modes"][mode] = r["FDR"]
            long_rows.append({
                "file_name": cond, "category": cat, "GMP_Cor": gmp, "mode": mode,
                "GO_ID": r["GO_ID"], "Term": r["Term"], "Category": r["Category"],
                "FDR": r["FDR"], "Ratio_in_study": r.get("Ratio_in_study", ""),
            })

    rows = []
    for go_id, info in term_info.items():
        ms = sorted(info["modes"])
        rows.append({
            "GO_ID": go_id,
            "Term": info["Term"],
            "Category": info["Category"],
            "modes": ",".join(map(str, ms)),
            "n_modes": len(ms),
            "best_FDR": min(info["modes"].values()),
            "FDR_by_mode": "; ".join(f"m{x}={info['modes'][x]:.1e}" for x in ms),
        })
    cols = ["GO_ID", "Term", "Category", "modes", "n_modes", "best_FDR", "FDR_by_mode"]
    cond_df = (pd.DataFrame(rows).sort_values("best_FDR").reset_index(drop=True)
               if rows else pd.DataFrame(columns=cols))
    cond_df.to_csv(os.path.join(OUT_DIR, cond), index=False)

    block = [f"\n=== {cond}  [{CAT_NAME.get(cat, cat)}]  GMP-Cor={gmp:.2f} ==="]
    if len(cond_df):
        for _, row in cond_df.iterrows():
            block.append(f"  [{row['Category']}] {row['Term']}  "
                         f"(modes {row['modes']}; {row['FDR_by_mode']})")
    else:
        block.append("  (no significant GO terms in modes 1-5)")
    text_blocks.append("\n".join(block))
    print("\n".join(block))

pd.DataFrame(long_rows).to_csv(os.path.join(OUT_DIR, "all_conditions.csv"), index=False)

header = (
    "Per-mode GO enrichment of leading eigenvectors (modes 1-5), consolidated per "
    "condition (Reviewer #4, comment 1).\n"
    "Source: top 50 genes by |loading| per mode vs full gene-panel background, "
    "BH-FDR < 0.05.\n"
)
with open(os.path.join(OUT_DIR, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(header)
    fh.write("\n".join(text_blocks))

print("\n" + header)
print(f"Wrote outputs to {OUT_DIR}")
