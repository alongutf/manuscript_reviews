# Figure 3 — Implementation Details

## How to run

Run from `scripts/figures/` so that `os.getcwd()` resolves correctly:

```bash
cd scripts/figures
python figure3.py          # interactive window
# or save directly:
# uncomment pf.save("figure3.svg", dpi=300, transparent=True) at the bottom
```

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | Line 16 — base font size (titles = 10pt, labels/ticks/legends = 8pt) |
| `figsize` | `(7, 5.5)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(-0.04, 0.055)` | `PanelFigure(label_offset=...)` — shifts A/B/C… labels left and above panel top-left |
| `root_dir` | 2 dirs above `os.getcwd()` | Automatically resolves to repo root when running from `scripts/figures/` |
| `ev_data_dir` | `root_dir/ev_data/` | All `.npy` eigenvalue files live here |
| `RESULTS_DIR` | `root_dir/results/simulation_results/` | Calibration curve data |

---

## Global font standard

| Element | Size |
|---|---|
| Panel titles | `fsize` = 10 |
| Axis labels, tick labels, legends | `fsize - 2` = 8 |

---

## Layout map

Panel positions are `[left, bottom, width, height]` in **figure-normalized coordinates** (0–1).  
The figure is 7 × 5 inches.

```
 0.0                                                        1.0
 ┌────────────────────────────────────────────────────────────┐ 0.94
 │  A (0.08,0.57)     B (0.39,0.57)     E (0.70,0.57)         │
 │  [0.25 × 0.37]     [0.25 × 0.37]     [0.27 × 0.37]         │
 │  Reg-Arrest CCDF   Dis-Arrest CCDF   Calibration curve    │
 ├──────────────────────────────────────────────────────────── 0.47
 │  C (0.08,0.10)     D (0.39,0.10)     F (0.70,0.10)         │
 │  [0.25 × 0.37]     [0.25 × 0.37]     [0.27 × 0.37]         │
 │  Sim ρ=0.9 CCDF    Sim ρ=0.0 CCDF    Box plot              │
 └────────────────────────────────────────────────────────────┘ 0.10
```

---

## Shared helper: `_plot_ccdf(ax, npy_file, title)`

Loads a `.npy` file from `ev_data_dir`. File shape: `(2, N)` — row 0 = original eigenvalues, row 1 = scrambled.

| Parameter | Value | Notes |
|---|---|---|
| Scale | loglog | |
| Threshold `x2` | `max(scrambled)` | Dashed vertical line |
| Noise line | darkgray, α=0.7, markersize=3 | eigenvalues < x2 |
| Signal line | skyblue, markersize=3 | eigenvalues ≥ x2 |
| Scrambled line | black, α=0.5, markersize=3 | |
| CCDF formula | `1 - rank/p + 1/p` | Log-safe; computed independently for data and scrambled |
| x limits | `[0.1, max(data) × 1.5]` | |

**To swap dataset:** change the filename argument.  
**To change x-range or markersize:** edit `_plot_ccdf` directly — applies to all CCDF panels.

---

## Panel A — Reg-Arrest experimental data

**Function:** `panel_A(ax)`  
**Data file:** `ev_data/sample_15b_filtered.npy`  
**Title:** `'Reg-Arrest'`

Calls `_plot_ccdf` directly.

---

## Panel B — Dis-Arrest experimental data

**Function:** `panel_B(ax)`  
**Data file:** `ev_data/sample_15a_filtered.npy`  
**Title:** `'Dis-Arrest'`

Calls `_plot_ccdf` directly.

---

## Panel C — Simulated data, high ρ

**Function:** `panel_C(ax)`  
**Data file:** `ev_data/simulated_pcs_0.9.npy` (ρ = 0.9)  
**Title:** `r'Simulation (high $\rho$)'`

To use a different ρ value, change the filename to `simulated_pcs_0.8.npy` etc. Available files: `simulated_pcs_0.npy` through `simulated_pcs_1.npy` (step 0.1).

---

## Panel D — Simulated data, low ρ

**Function:** `panel_D(ax)`  
**Data file:** `ev_data/simulated_pcs_0.npy` (ρ = 0)  
**Title:** `r'Simulation (low $\rho$)'`

---

## Panel E — GMP-Cor calibration curve

**Function:** `panel_E(ax)`  
**Data file:** `results/simulation_results/raw/rho_sweep_summary.txt`

| Parameter | Value | Notes |
|---|---|---|
| x-axis | ρ values from summary file (0.00–1.00) | Read from index column |
| y-axis | GMP-Cor median ± SD | `median` and `std` columns |
| x ticks | `np.arange(0, 1.05, 0.2)` | Every 0.2 step, no rotation |
| x limits | `[-0.02, 1.02]` | |
| y limits | `bottom=-1, top=75` | Clips the steep rise at ρ→1 |
| Grid | dashed `'--'`, α=0.2 | |
| Fill band | ±SD, `steelblue`, α=0.15 | |
| Title | `'GMP-Cor calibration curve'` | 10pt |
| Axis labels | `fsize - 2` = 8pt | |

---

## Panel F — Box plot: Regulated vs. Dis-Arrest

**Function:** `panel_F(ax)`  
**Data file:** `results/data_metrics/test8.csv` (index_col=0)

| Parameter | Value | Notes |
|---|---|---|
| Metric | `sum_denoised_ev` column | GMP-Cor proxy |
| group1 | `category == 'r'` | Regulated samples |
| group0 | `category == 'd'` | Dis-Arrest samples |
| Box style | `patch_artist=True`, facecolor=None | Hollow boxes |
| Means | shown as solid line | `meanline=True`, `showmeans=True` |
| Significance | Mann-Whitney U test (`stats.mannwhitneyu`) | Annotated as `*`/`**`/`***`/`****`/`NS` |
| Annotation y | `max(group1, group0) + 1`, h=1 | Bracket height and spacing |
| y-label | `'GMP-Cor'` | |
| x-tick labels | `['Regulated', 'Dis-Arrest']` | |

**`format_p(p)` thresholds:**  
p < 0.0001 → `****`, p < 0.001 → `***`, p < 0.01 → `**`, p < 0.05 → `*`, else → `NS`

---

## Data files used in figure

| File | Used in | Description |
|---|---|---|
| `ev_data/sample_15b_filtered.npy` | Panel A | Reg-Arrest experimental |
| `ev_data/sample_15a_filtered.npy` | Panel B | Dis-Arrest experimental |
| `ev_data/simulated_pcs_0.9.npy` | Panel C | Simulation ρ = 0.9 (high) |
| `ev_data/simulated_pcs_0.npy` | Panel D | Simulation ρ = 0 (low) |
| `results/simulation_results/raw/rho_sweep_summary.txt` | Panel E | Calibration sweep summary |
| `results/data_metrics/test8.csv` | Panel F | Dataset metrics with categories |

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize the whole figure | `figsize=(W, H)` in `PanelFigure(...)` |
| Change all font sizes | `fsize = 10` at top of script |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Swap CCDF dataset (A/B/C/D) | change the filename string in the panel function |
| Use different simulation ρ | replace `simulated_pcs_0.9.npy` or `simulated_pcs_0.npy` with another ρ value |
| Adjust CCDF x-range | change `max(d1s)*1.5` multiplier in `_plot_ccdf` |
| Change calibration x-tick density | edit `np.arange(0, 1.05, 0.2)` in `panel_E` |
| Save as SVG | uncomment `pf.save("figure3.svg", dpi=300, transparent=True)` |