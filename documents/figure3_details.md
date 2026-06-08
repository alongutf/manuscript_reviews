# Figure 3 — Implementation Details

## How to run

Run from `scripts/figures/` so that `os.getcwd()` resolves correctly:

```bash
cd scripts/figures
python figure3.py          # interactive window
# SVG is saved automatically at the bottom of the script
```

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | Line 18 — base font size (titles = 10pt, labels/ticks/legends = 8pt) |
| `figsize` | `(7, 7)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(-0.04, 0.02)` | `PanelFigure(label_offset=...)` — shifts A/B/C… labels left and above panel top-left |
| `root_dir` | 2 dirs above `os.getcwd()` | Automatically resolves to repo root when running from `scripts/figures/` |
| `ev_data_dir` | `root_dir/ev_data/` | All `.npy` eigenvalue files live here |
| `RESULTS_DIR` | `root_dir/results/simulation_results/` | Calibration curve data |
| `REG_COLOR` | `'steelblue'` | Condition color for Regulated/high-ρ — used in panels A, B, G, H |
| `DIS_COLOR` | `'#E07B54'` | Condition color for Dis-Arrest/low-ρ — used in panels D, E, G, H |
| `_HEATMAP_PARAMS` | `dict(n=2000, shape=1.5, hub_probability=0.2, seed=31)` | Parameters for `generate_gram_hub_matrix` — match the `rho_sweep` run that generated the `.npy` files |
| `med_reg`, `med_dis` | computed at module level from `test8.csv` | Group medians shared by panels G and H |

---

## Global font standard

| Element | Size |
|---|---|
| Panel titles | `fsize` = 10 |
| Axis labels, tick labels, legends | `fsize - 2` = 8 |

---

## Layout map

Panel positions are `[left, bottom, width, height]` in **figure-normalized coordinates** (0–1).  
The figure is 7 × 7 inches.

```
 0.0                                                              1.0
 ┌──────────────────────────────────────────────────────────────────┐ ~0.96
 │  A (0.08, 0.72)          B (0.44, 0.72)        C (0.77, 0.79)  │
 │  [0.28 × 0.24]           [0.28 × 0.24]          [0.17 × 0.17]  │
 │  Reg-Arrest CCDF         Sim ρ=0.9 CCDF         Heatmap ρ=0.9  │
 ├────────────────────────────────────────────────────────────────── ~0.64
 │  D (0.08, 0.40)          E (0.44, 0.40)          F (0.77, 0.47) │
 │  [0.28 × 0.24]           [0.28 × 0.24]           [0.17 × 0.17] │
 │  Dis-Arrest CCDF         Sim ρ=0.5 CCDF          Heatmap ρ=0.5 │
 ├────────────────────────────────────────────────────────────────── ~0.30
 │  G (0.18, 0.06)                    H (0.54, 0.06)               │
 │  [0.28 × 0.24]                     [0.28 × 0.24]                │
 │  GMP-Cor calibration curve         Box plot                     │
 └──────────────────────────────────────────────────────────────────┘  0.06
```

Panels C and F are small heatmap insets (0.17 × 0.17) placed at the right of each row, vertically centered within their respective row. Panels G and H share the y-axis via `ax_H.sharey(ax_G)`.

---

## Shared helper: `_plot_ccdf(ax, npy_file, title, signal_color)`

Loads a `.npy` file from `ev_data_dir`. File shape: `(2, N)` — row 0 = original eigenvalues, row 1 = scrambled.

| Parameter | Value | Notes |
|---|---|---|
| Scale | loglog | |
| Threshold `x2` | `max(scrambled)` | Dashed vertical line |
| Threshold label | `λ_max^scr` at axes y=0.8, fontsize `fsize-2` | Uses `ax.get_xaxis_transform()` for data-x / axes-y positioning |
| Noise line | darkgray, α=0.7, markersize=3 | eigenvalues < x2 |
| Signal line | `signal_color`, markersize=3 | eigenvalues ≥ x2; color varies by panel condition |
| Scrambled line | black, α=0.5, markersize=3 | |
| CCDF formula | `1 - rank/p + 1/p` | Log-safe; computed independently for data and scrambled |
| x limits | `[0.1, 30]` | Fixed range shared across all CCDF panels |
| Legend | individual per-panel | noise / signal / scrambled |

**To swap dataset:** change the filename argument.  
**To change signal color:** pass `signal_color=` — `REG_COLOR` for Regulated/high-ρ, `DIS_COLOR` for Dis-Arrest/low-ρ.  
**To change x-range:** edit `ax.set_xlim([0.1, 30])` in `_plot_ccdf` — applies to all CCDF panels.

---

## Shared helper: `_plot_heatmap(ax, rho, title)`

Generates a correlation matrix via `generate_gram_hub_matrix` (imported from `src.simulations`) and displays the top-left 100×100 submatrix.

| Parameter | Value | Notes |
|---|---|---|
| Matrix | `generate_gram_hub_matrix(alpha=rho, **_HEATMAP_PARAMS)[:100, :100]` | First 100 genes; cluster structure is visible because genes within clusters are contiguous |
| Colormap | `RdBu_r`, vmin=-1, vmax=1 | Diverging; positive = red, negative = blue |
| Colorbar | yes, `fraction=0.046`, `pad=0.04` | Tick labels at `fsize-2`; no text label |
| Interpolation | `'nearest'` | Preserves discrete pixel boundaries |
| Axis ticks | none | `set_xticks([])`, `set_yticks([])` |
| Axis labels | none | |

**To show the full matrix:** remove `[:100, :100]` slicing.  
**To use a different ρ:** change the `rho=` argument in the panel function call.  
**To change matrix structure:** edit `_HEATMAP_PARAMS` at the top of the script.

---

## Panel A — Reg-Arrest experimental data

**Function:** `panel_A(ax)`  
**Data file:** `ev_data/sample_15b_filtered.npy`  
**Title:** `'Reg-Arrest'`  
**Signal color:** `REG_COLOR` (steelblue)

---

## Panel B — Simulated data, high ρ

**Function:** `panel_B(ax)`  
**Data file:** `ev_data/simulated_pcs_0.9.npy` (ρ = 0.9)  
**Title:** `r'Simulation ($\rho=0.9$)'`  
**Signal color:** `REG_COLOR` (steelblue)

---

## Panel C — Covariance matrix heatmap, ρ = 0.9

**Function:** `panel_C(ax)`  
**Title:** `r'$\rho=0.9$'`  
**Position:** Small inset (0.17 × 0.17) at the right of row 1, vertically centered.

---

## Panel D — Dis-Arrest experimental data

**Function:** `panel_D(ax)`  
**Data file:** `ev_data/sample_15a_filtered.npy`  
**Title:** `'Dis-Arrest'`  
**Signal color:** `DIS_COLOR` (#E07B54 coral)

---

## Panel E — Simulated data, low ρ

**Function:** `panel_E(ax)`  
**Data file:** `ev_data/simulated_pcs_0.5.npy` (ρ = 0.5)  
**Title:** `r'Simulation ($\rho=0.5$)'`  
**Signal color:** `DIS_COLOR` (#E07B54 coral)

---

## Panel F — Covariance matrix heatmap, ρ = 0.5

**Function:** `panel_F(ax)`  
**Title:** `r'$\rho=0.5$'`  
**Position:** Small inset (0.17 × 0.17) at the right of row 2, vertically centered.

---

## Panel G — GMP-Cor calibration curve

**Function:** `panel_G(ax)`  
**Data file:** `results/simulation_results/raw/rho_sweep_summary.txt`

| Parameter | Value | Notes |
|---|---|---|
| x-axis | ρ values (0.00–1.00) | Read from index column |
| y-axis | GMP-Cor median ± SD | `median` and `std` columns |
| y scale | linear | |
| y limits | `bottom=-1, top=45` | |
| x ticks | `np.arange(0, 1.05, 0.2)` | Every 0.2 step |
| x limits | `[-0.02, 1.02]` | |
| Grid | dashed `'--'`, α=0.2 | |
| Fill band | ±SD, `steelblue`, α=0.3 | |
| Reference lines | `med_reg` (REG_COLOR dashed), `med_dis` (DIS_COLOR dashed) | Connects to panel H — shows where each experimental group falls on the curve |
| Title | `'GMP-Cor calibration curve'` | 10pt |
| Axis labels | `fsize - 2` = 8pt | |

---

## Panel H — Box plot: Regulated vs. Dis-Arrest

**Function:** `panel_H(ax)`  
**Data file:** `results/data_metrics/test8.csv` (index_col=0)

| Parameter | Value | Notes |
|---|---|---|
| Metric | `sum_denoised_ev` column | GMP-Cor proxy |
| group1 | `category == 'r'` | Regulated samples |
| group0 | `category == 'd'` | Dis-Arrest samples |
| Box edge color | `REG_COLOR` / `DIS_COLOR` | Condition-coded; facecolor α=0.12 tint |
| Strip plot | `s=10`, α=0.5, jitter ±0.12, seed=42 | Individual samples overlaid on boxes |
| Means | shown as solid line | `meanline=True`, `showmeans=True` |
| Reference lines | `med_reg` (REG_COLOR dashed), `med_dis` (DIS_COLOR dashed) | Same y-values as panel G lines |
| Significance | Mann-Whitney U test (`stats.mannwhitneyu`) | Annotated as `*`/`**`/`***`/`****`/`NS` |
| Annotation y | `max(group1, group0) + h*0.5`, h=1 | |
| y limits | `bottom=-1, top=45` | Fixed to match panel G |
| y-label | `''` | Shared with panel G via `sharey` |
| x-tick labels | `['Regulated', 'Dis-Arrest']` | |
| Shared y-axis | `ax_H.sharey(ax_G)` | Called after both panels are drawn |

**`format_p(p)` thresholds:**  
p < 0.0001 → `****`, p < 0.001 → `***`, p < 0.01 → `**`, p < 0.05 → `*`, else → `NS`

---

## Data files used in figure

| File | Used in | Description |
|---|---|---|
| `ev_data/sample_15b_filtered.npy` | Panel A | Reg-Arrest experimental |
| `ev_data/simulated_pcs_0.9.npy` | Panel B | Simulation ρ = 0.9 |
| `ev_data/sample_15a_filtered.npy` | Panel D | Dis-Arrest experimental |
| `ev_data/simulated_pcs_0.5.npy` | Panel E | Simulation ρ = 0.5 |
| `results/simulation_results/raw/rho_sweep_summary.txt` | Panel G | Calibration sweep summary |
| `results/data_metrics/test8.csv` | Panel H | Dataset metrics with categories |

Panels C and F generate their matrices at runtime via `generate_gram_hub_matrix` — no pre-saved file.

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize the whole figure | `figsize=(W, H)` in `PanelFigure(...)` |
| Change all font sizes | `fsize = 10` at top of script |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Move panel label | `label_offset=(-0.04, 0.02)` in `PanelFigure(...)` |
| Swap CCDF dataset | change the filename string in the panel function |
| Change condition colors | `REG_COLOR` / `DIS_COLOR` at top of script — propagates to all panels |
| Use different simulation ρ (CCDF) | replace `simulated_pcs_0.9.npy` or `simulated_pcs_0.5.npy` with another value |
| Use different simulation ρ (heatmap) | change `rho=` argument in `panel_C` or `panel_F` |
| Show full heatmap (2000×2000) | remove `[:100, :100]` in `_plot_heatmap` |
| Change heatmap structure | edit `_HEATMAP_PARAMS` at top of script |
| Fix x-range across all CCDF panels | edit `ax.set_xlim([0.1, 30])` in `_plot_ccdf` |
| Move threshold label vertically | change `0.8` in `ax.text(x2 * 1.1, 0.8, ...)` — axes coords (0=bottom, 1=top) |
| Change calibration x-tick density | edit `np.arange(0, 1.05, 0.2)` in `panel_G` |
| Re-enable separate y-axis on H | uncomment the two `ax_H` lines after `ax_H.sharey(ax_G)` |
| Save as SVG | `pf.save("figure3.svg", dpi=300, transparent=True)` at bottom of script |
