# Bulk RNA-seq normalization to molecules per cell (ERCC spike-in)

**Script:** `scripts/normalize_bulk_to_molecules_per_cell.py`
**Input matrix:** `bulk_data/exp0224_count_data.csv` (raw counts, *features × samples* — never modified)
**Outputs:**
- `bulk_data/exp0224_molecules_per_cell.csv` — normalized matrix (molecules per cell)
- `bulk_data/exp0224_ercc_normalization_summary.csv` — per-sample QC / audit table

---

## 1. Goal

Convert raw read counts into **absolute molecules per cell** for every feature in
every sample. ERCC spike-in controls are a defined number of synthetic RNA
molecules added to each sample at a fixed ratio to total RNA. Because we know how
many ERCC molecules went in and how many reads came out, the ERCC reads give a
per-sample calibration from "reads" to "molecules". Dividing by the number of
cells in each sample then gives molecules per cell.

---

## 2. Parameters (top of the script)

| Parameter | Value used | Meaning |
|---|---|---|
| `ERCC_MIX` | `"Mix 1"` | Which ERCC mix was spiked in. Concentrations read from `metadata/ERCC_controls_analysis.txt`. |
| `TOTAL_MIX_VOLUME_NL` | `20.0` | **Fixed total volume of ERCC mix added per sample, in nl.** When set, this is used directly and the per-ng recipe below is ignored. Set to `None` to use the per-ng recipe instead. |
| `SPIKE_VOLUME_NL_PER_NG` | `5.0` | (Per-ng recipe, only used when `TOTAL_MIX_VOLUME_NL is None`.) nl of ERCC mix per ng of total RNA. |
| `TOTAL_RNA_NG` | `1000.0` | (Per-ng recipe, only used when `TOTAL_MIX_VOLUME_NL is None`.) Total RNA per sample, in ng. |
| `AVOGADRO` | `6.02214076e23` | Molecules per mole. |
| `CELL_COUNTS` | EXP = 1×10⁶, VapC = 3×10⁶, CASP = 7×10⁶ | Estimated cells per sample, matched by substring in the column name. |

**Current run:** a fixed **20 nl of Mix 1 per sample** (`TOTAL_MIX_VOLUME_NL = 20`),
independent of RNA mass.

> **The two assumptions to revisit:** the spiked mix volume (`TOTAL_MIX_VOLUME_NL`,
> or `TOTAL_RNA_NG` if using the per-ng recipe) and `CELL_COUNTS`. Both are
> external to the data and set the absolute scale. See §7.

---

## 3. Step-by-step method

### Step 1 — Load the raw count matrix
`bulk_data/exp0224_count_data.csv` is read with features as rows (index) and the
24 samples as columns. 4590 features, including all 92 ERCC controls.

### Step 2 — Load ERCC concentrations
From `metadata/ERCC_controls_analysis.txt` we take, for each ERCC ID, the
**Mix 1 concentration** in attomoles/µl.

### Step 3 — Molecules of each ERCC spiked into one sample
First the volume of mix added per sample. With a fixed total volume:

```
spike_volume_ul = TOTAL_MIX_VOLUME_NL / 1000 = 20 nl / 1000 = 0.02 µl
```

(If instead the per-ng recipe is used, `spike_volume_ul = SPIKE_VOLUME_NL_PER_NG
× TOTAL_RNA_NG / 1000`.)

Then, for each ERCC *i*, the absolute number of molecules added:

```
molecules_spiked_i = conc_i [attomol/µl] × spike_volume_ul [µl]
                     × 1e-18 [mol/attomol] × N_A [molecules/mol]
```

- `1 attomole = 1e-18 mol`, and `× N_A` converts moles → molecules.
- With 0.02 µl (20 nl) of Mix 1, the **total** spiked across all 92 ERCC is
  **≈ 1.25 × 10⁹ molecules per sample**.

Because the spiked volume is the same for every sample, the spiked molecule count
is **identical across samples**. The per-sample differences in recovered ERCC
*reads* therefore reflect only capture efficiency and sequencing depth — exactly
the quantity we want to calibrate out.

### Step 4 — Restrict to ERCC present in the matrix
All 92 ERCC from the metadata are present in the count matrix, so all 92 are used
for calibration. (If any were missing the script warns and uses the intersection.)

### Step 5 — Per-sample conversion factor `k_s` (molecules per read)
For each sample *s*:

```
k_s = (Σ_i molecules_spiked_i) / (Σ_i ERCC_reads_{i,s})
```

i.e. total spiked ERCC molecules ÷ total ERCC reads in that sample.

- This **sum-based** estimator is dominated by the high-concentration ERCC, which
  sit well above the detection limit, so it is robust to drop-out of the rare,
  low-concentration ERCC (those contribute negligibly to both sums).
- As a QC check only (not used in the conversion), the script also computes the
  **log–log correlation R²** between input molecules and recovered reads per
  sample. Observed R² ≈ 0.70–0.85, confirming a roughly linear spike-in response.

### Step 6 — Map each sample to a cell count
The cell count is chosen by matching the sample name (case-insensitive substring):
`EXP → 1e6`, `VapC → 3e6`, `CASP → 7e6`. Every one of the 24 columns matches
exactly one group.

### Step 7 — Convert to molecules, then molecules per cell
For every feature *g* and sample *s*:

```
molecules(g,s)          = counts(g,s) × k_s
molecules_per_cell(g,s) = molecules(g,s) / n_cells(s)
```

The result is written to `bulk_data/exp0224_molecules_per_cell.csv` with the same
shape and labels as the input. ERCC rows are retained in the output for
transparency (their normalized values are not biologically meaningful — they are
the calibration standards — and can be dropped downstream).

---

## 4. Output files

**`exp0224_molecules_per_cell.csv`** — features × samples, in molecules per cell.

**`exp0224_ercc_normalization_summary.csv`** — one row per sample:

| Column | Meaning |
|---|---|
| `ercc_reads` | Total reads mapping to ERCC controls. |
| `total_reads` | Total reads (all features). |
| `ercc_read_fraction` | ERCC reads / total reads (0.04%–0.26% here — a normal range). |
| `molecules_per_read_k` | The per-sample conversion factor `k_s`. |
| `n_cells` | Cells assigned to the sample. |
| `loglog_r2_qc` | QC: log–log R² of ERCC reads vs input molecules. |

---

## 5. Worked sanity check (sample `T0CASP2`)

- Total ERCC molecules spiked = 1.247 × 10⁹; total ERCC reads = 20,848.
- `k = 1.247e9 / 20848 = 5.98 × 10⁴` molecules/read.
- ERCC-00002 raw count = 8983 → `8983 × 5.98e4 / 7e6 ≈ 76.7` molecules/cell
  (matches the output file).

---

## 6. Units and key conversions

- `1 attomole (attomol) = 1 × 10⁻¹⁸ mol`
- `molecules = mol × 6.02214076 × 10²³`
- `1 µl = 1000 nl`, `1 µg = 1000 ng`

---

## 7. Assumptions and how to change them

1. **A fixed 20 nl of mix is added to every sample.** Set by
   `TOTAL_MIX_VOLUME_NL`. The spiked molecule count, and therefore the absolute
   molecules-per-cell scale, is directly proportional to this volume. If the true
   volume differs, change the parameter; if it differed *between* samples, the
   current single-value assumption is not adequate and a per-sample volume would
   be required. (To switch to a per-ng recipe instead, set `TOTAL_MIX_VOLUME_NL =
   None` and use `SPIKE_VOLUME_NL_PER_NG` × `TOTAL_RNA_NG`.)
2. **Cell counts** are estimates per condition (`CELL_COUNTS`). They scale each
   group's molecules-per-cell inversely.
3. **Linear spike-in response.** A single conversion factor per sample assumes
   reads scale linearly with input molecules; the log–log R² QC supports this.
4. **Equal capture of ERCC and endogenous RNA.** Standard ERCC-normalization
   assumption — endogenous molecules are converted with the ERCC-derived factor.

To re-run after changing any parameter:

```bash
python scripts/normalize_bulk_to_molecules_per_cell.py
```

The original count matrix is never modified; outputs are overwritten in place.