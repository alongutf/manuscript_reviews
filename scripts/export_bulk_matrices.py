"""
Export the two bulk RNA-seq matrices behind the correlation figures as plain
genes x samples tables.

  bulk_matrix_expression<SUFFIX>.csv   absolute expression, log10(normalised count + 1)
  bulk_matrix_lfc<SUFFIX>.csv          log2 fold change vs each sample's own
                                       exponential control

Same 12 samples, same expressed-gene set and same row order in both files, so they
line up cell for cell.  Columns carry the publication sample labels and are ordered
condition > batch > sample, matching the correlation figures; the machine-readable
<dataset>|<sample> ids stay available in the inputs.

Input:  results/bulk_pca/bulk_log_fold_changes<SUFFIX>.csv
        results/bulk_pca/bulk_normalized_counts<SUFFIX>.csv
        written by  python bulk_lfc_pca.py --bulk-control matched --shx-late

Usage:
    python export_bulk_matrices.py
"""

import argparse
import os
import re
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "bulk_pca")

parser = argparse.ArgumentParser()
parser.add_argument("--suffix", default="_log2_matched_late",
                    help="suffix of the input matrices; the default is the 12-sample "
                         "matched-control / late-SHX set")
args = parser.parse_args()


def load(stem):
    src = os.path.join(OUT, "%s%s.csv" % (stem, args.suffix))
    if not os.path.exists(src):
        raise SystemExit("%s not found - run bulk_lfc_pca.py --bulk-control matched "
                         "--shx-late first" % src)
    return pd.read_csv(src, index_col=0)          # genes x samples


lfc = load("bulk_log_fold_changes")
# +1 keeps the handful of zero counts finite; every gene here already passes the
# mean-count filter, so the pseudocount moves nothing that carries signal
expr = np.log10(load("bulk_normalized_counts") + 1.0)
if list(lfc.columns) != list(expr.columns) or list(lfc.index) != list(expr.index):
    raise SystemExit("the two matrices do not hold the same genes and samples")

# ------------------------------------------------- sample annotation (as in the PCA)
LABELS = {"bulk": {"Disrupted": "Dis-Arrest1", "CASP": "Reg-Arrest1"},
          "shx": "Dis-Arrest2", "casp": "Reg-Arrest2+SHX"}


def annotate(col):
    ds, s = col.split("|")
    label = LABELS["bulk"]["CASP" if s.startswith("CASP") else "Disrupted"] \
        if ds == "bulk" else LABELS[ds]
    name = re.sub(r"^(CASP|Disrupted)_biorep", "rep", s)
    name = re.sub(r"^CASP biorep1 (\d+)min", r"\1min", name)
    return {"label": label, "name": "%s %s" % (label, name),
            "condition": "Reg-Arrest" if label.startswith("Reg") else "Dis-Arrest",
            "batch": "bulk" if ds == "bulk" else "timecourse"}


meta = pd.DataFrame([annotate(c) for c in lfc.columns], index=lfc.columns)

COND_ORDER = ["Reg-Arrest", "Dis-Arrest"]
BATCH_ORDER = ["bulk", "timecourse"]
order = sorted(meta.index, key=lambda c: (COND_ORDER.index(meta.loc[c, "condition"]),
                                          BATCH_ORDER.index(meta.loc[c, "batch"]),
                                          meta.loc[c, "name"]))

for stem, mat in (("expression", expr), ("lfc", lfc)):
    out = mat[order].copy()
    out.columns = [meta.loc[c, "name"] for c in order]
    out.index.name = "gene"
    dst = os.path.join(OUT, "bulk_matrix_%s%s.csv" % (stem, args.suffix))
    out.to_csv(dst, float_format="%.6g")
    print("%-10s %d genes x %d samples -> %s" % (stem, out.shape[0], out.shape[1],
                                                 os.path.basename(dst)))
print("columns:", ", ".join(meta.loc[c, "name"] for c in order))
