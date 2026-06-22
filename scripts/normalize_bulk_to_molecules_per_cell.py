"""
Normalize the bulk RNA-seq count matrix to absolute molecules per cell using
ERCC spike-in controls.

Input  : bulk_data/exp0224_count_data.csv      (raw counts, features x samples)
         metadata/ERCC_controls_analysis.txt   (ERCC Mix 1/Mix 2 concentrations)
Output : bulk_data/exp0224_molecules_per_cell.csv   (normalized matrix)
         bulk_data/exp0224_ercc_normalization_summary.csv  (per-sample QC table)

The original count matrix is NEVER modified.

See documents/bulk_ercc_normalization.md for a full written explanation of the
method, the formulas, and the assumptions behind every parameter.

Run from the scripts/ directory (or anywhere; paths are resolved relative to the
repository root inferred from this file's location).
"""

import os
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# PARAMETERS  (edit these — they are the only knobs you should ever need)
# --------------------------------------------------------------------------- #

# ERCC spike-in recipe -------------------------------------------------------
ERCC_MIX = "Mix 1"                 # which mix was spiked in ("Mix 1" or "Mix 2")

# There are two ways to define how much mix was added.  If TOTAL_MIX_VOLUME_NL is
# set (not None) it is used directly as the total volume of mix per sample, and
# the per-ng recipe below is ignored.  Set it to None to fall back to the
# per-ng recipe (SPIKE_VOLUME_NL_PER_NG x TOTAL_RNA_NG).
TOTAL_MIX_VOLUME_NL = 20.0         # total nl of ERCC mix added per sample (or None)

SPIKE_VOLUME_NL_PER_NG = 5.0       # nl of ERCC mix added per ng of total RNA
TOTAL_RNA_NG = 1000.0              # total RNA per sample (ng).  1 ug = 1000 ng.
                                   # <-- change this if the real input mass differs.

# Physical constant ----------------------------------------------------------
AVOGADRO = 6.02214076e23           # molecules per mole

# Cell-count estimates per sample group --------------------------------------
# Matched against the sample (column) name, case-insensitively, by substring.
# Order matters only if a name could match two keys; these are mutually exclusive.
CELL_COUNTS = {
    "EXP":  8e7,   # exponential-phase samples
    "VapC": 1e8,   # VapC samples
    "CASP": 3e8,   # casp samples
}

# File names -----------------------------------------------------------------
COUNTS_CSV  = "bulk_data/exp0224_count_data.csv"
ERCC_TXT    = "metadata/ERCC_controls_analysis.txt"
OUT_MATRIX  = "bulk_data/exp0224_molecules_per_cell.csv"
OUT_SUMMARY = "bulk_data/exp0224_ercc_normalization_summary.csv"

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def p(rel):
    return os.path.join(REPO_ROOT, rel)


# --------------------------------------------------------------------------- #
# 1. Load the raw count matrix (features x samples)
# --------------------------------------------------------------------------- #
counts = pd.read_csv(p(COUNTS_CSV), index_col=0)
counts.index.name = "feature"
samples = list(counts.columns)
print(f"Loaded count matrix: {counts.shape[0]} features x {counts.shape[1]} samples")

# --------------------------------------------------------------------------- #
# 2. Load ERCC concentrations and pick the spiked mix
# --------------------------------------------------------------------------- #
ercc = pd.read_csv(p(ERCC_TXT), sep="\t")
conc_col = f"concentration in {ERCC_MIX} (attomoles/ul)"
if conc_col not in ercc.columns:
    raise KeyError(f"Column '{conc_col}' not found. Available: {list(ercc.columns)}")
ercc = ercc[["ERCC ID", conc_col]].rename(columns={"ERCC ID": "feature",
                                                    conc_col: "conc_attomol_per_ul"})
ercc = ercc.set_index("feature")

# --------------------------------------------------------------------------- #
# 3. Molecules of each ERCC spiked into one sample
#
#    volume of mix per sample (ul) = (nl per ng) * (ng total RNA) / 1000
#    molecules_i = conc_i [attomol/ul] * volume [ul] * 1e-18 [mol/attomol] * N_A
#
#    Because total RNA is assumed identical across samples, the spiked molecule
#    count is the same for every sample (only the recovered ERCC *reads* differ,
#    which is exactly the per-sample capture/depth signal we calibrate on).
# --------------------------------------------------------------------------- #
if TOTAL_MIX_VOLUME_NL is not None:
    spike_volume_nl = TOTAL_MIX_VOLUME_NL
    recipe = f"{TOTAL_MIX_VOLUME_NL} nl total mix per sample (fixed volume)"
else:
    spike_volume_nl = SPIKE_VOLUME_NL_PER_NG * TOTAL_RNA_NG
    recipe = f"{SPIKE_VOLUME_NL_PER_NG} nl/ng x {TOTAL_RNA_NG} ng"
spike_volume_ul = spike_volume_nl / 1000.0
ercc["molecules_spiked"] = (
    ercc["conc_attomol_per_ul"] * spike_volume_ul * 1e-18 * AVOGADRO
)
print(f"ERCC mix volume per sample: {spike_volume_ul:g} ul  ({recipe})")
print(f"Total ERCC molecules spiked per sample: "
      f"{ercc['molecules_spiked'].sum():.3e}")

# --------------------------------------------------------------------------- #
# 4. Restrict to ERCC that are actually present in the count matrix
# --------------------------------------------------------------------------- #
ercc_in_matrix = [f for f in ercc.index if f in counts.index]
missing = sorted(set(ercc.index) - set(ercc_in_matrix))
if missing:
    print(f"WARNING: {len(missing)} ERCC in metadata not found in matrix: {missing}")
print(f"Using {len(ercc_in_matrix)} ERCC spike-ins for calibration")

ercc_counts = counts.loc[ercc_in_matrix]                 # ERCC reads, per sample
ercc_molecules = ercc.loc[ercc_in_matrix, "molecules_spiked"]   # molecules, per ERCC

# --------------------------------------------------------------------------- #
# 5. Per-sample conversion factor  k_s  =  molecules per read
#
#    k_s = (sum of spiked ERCC molecules) / (sum of ERCC reads in sample s)
#
#    The sum is dominated by the high-concentration ERCC, which are far above
#    the detection limit, so this estimator is robust to drop-out of the rare
#    low-concentration ERCC.  A log-log regression slope and R^2 are also
#    reported purely as a QC check on spike-in linearity.
# --------------------------------------------------------------------------- #
total_spiked = ercc_molecules.sum()
ercc_read_totals = ercc_counts.sum(axis=0)               # per sample
k = total_spiked / ercc_read_totals                      # molecules per read, per sample

# QC: log-log linearity of recovered reads vs input molecules, per sample
qc_r2 = {}
for s in samples:
    x = ercc_molecules.values
    y = ercc_counts[s].values
    mask = (x > 0) & (y > 0)
    if mask.sum() >= 3:
        lx, ly = np.log10(x[mask]), np.log10(y[mask])
        r = np.corrcoef(lx, ly)[0, 1]
        qc_r2[s] = r ** 2
    else:
        qc_r2[s] = np.nan

# --------------------------------------------------------------------------- #
# 6. Map each sample to a cell count
# --------------------------------------------------------------------------- #
def cells_for(sample_name):
    name = sample_name.lower()
    for key, val in CELL_COUNTS.items():
        if key.lower() in name:
            return val
    raise ValueError(f"No cell-count group matches sample '{sample_name}'. "
                     f"Known groups: {list(CELL_COUNTS)}")

n_cells = pd.Series({s: cells_for(s) for s in samples})

# --------------------------------------------------------------------------- #
# 7. Convert every feature to molecules, then to molecules per cell
#
#    molecules(g, s)          = counts(g, s) * k_s
#    molecules_per_cell(g, s) = molecules(g, s) / n_cells(s)
# --------------------------------------------------------------------------- #
molecules = counts.mul(k, axis=1)
molecules_per_cell = molecules.div(n_cells, axis=1)

molecules_per_cell.to_csv(p(OUT_MATRIX))
print(f"\nWrote normalized matrix -> {OUT_MATRIX}")

# --------------------------------------------------------------------------- #
# 8. Per-sample QC / summary table
# --------------------------------------------------------------------------- #
summary = pd.DataFrame({
    "ercc_reads": ercc_read_totals,
    "total_reads": counts.sum(axis=0),
    "ercc_read_fraction": ercc_read_totals / counts.sum(axis=0),
    "molecules_per_read_k": k,
    "n_cells": n_cells,
    "loglog_r2_qc": pd.Series(qc_r2),
})
summary.index.name = "sample"
summary.to_csv(p(OUT_SUMMARY))
print(f"Wrote per-sample summary  -> {OUT_SUMMARY}\n")
pd.set_option("display.width", 200, "display.max_columns", 20)
print(summary.round(4))