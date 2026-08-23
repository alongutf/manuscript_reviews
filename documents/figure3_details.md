# Figure 3 — Implementation Details

## How to run

Run from `scripts/figures/` so that `os.getcwd()` resolves correctly:

```bash
cd scripts/figures
python figure3.py          # interactive window
# figure3.pdf + figure3_preview.png are saved automatically at the bottom
```

---

## Layout at a glance

Figure 3 is a **2×2 grid of CCDF panels (A–D)** on the left with the **GMP-Cor box
plot (E) to the right**, spanning the full height of the A–D block.

The GMP-Cor **calibration curve** and the **box plot** used to be the bottom row
(old panels E/F). They now live side-by-side in **Supplementary Figure S10**
(`scripts/supplementary_figures/figure_s10.py`). The box plot is reproduced in
both figures (panel E here, panel B in S10).

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | base font size (titles = 10pt, labels/ticks/legends = 8pt) |
| `figsize` | `(9.5, 4.5)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(-0.04, 0.02)` | `PanelFigure(label_offset=...)` — shifts A/B/C… labels left and above panel top-left |
| `root_dir` | 2 dirs above `os.getcwd()` | Automatically resolves to repo root when running from `scripts/figures/` |
| `ev_data_dir` | `root_dir/ev_data/` | All `.npy` eigenvalue files live here |
| `REG_COLOR` | `'steelblue'` | Condition color for Regulated / high-χ — panels A, C, E |
| `DIS_COLOR` | `'#E07B54'` | Condition color for Dis-Arrest / low-χ — panels B, D, E |
| `REF_LW`, `REF_ALPHA` | `1`, `0.85` | Style of the median reference lines in the box plot |
| `med_reg`, `med_dis` | computed at module level from `test8.csv` | Group medians drawn as reference lines in panel E |

---

## Global font standard

| Element | Size |
|---|---|
| Panel titles | `fsize` = 10 |
| Axis labels, tick labels, legends | `fsize - 2` = 8 |

---

## Layout map

Panel positions are `[left, bottom, width, height]` in **figure-normalized coordinates** (0–1).
The figure is 9.5 × 4.5 inches.

```
 0.0                                                          1.0
 ┌──────────────────────────────────────────────────────────────┐ ~0.93
 │  A (0.06, 0.57)      B (0.37, 0.57)                          │
 │  [0.28 × 0.36]       [0.28 × 0.36]        E (0.78, 0.11)     │
 │  Reg-Arrest CCDF     Dis-Arrest CCDF      [0.17 × 0.82]      │
 │                                           Box plot           │
 │  C (0.06, 0.11)      D (0.37, 0.11)       (Regulated vs      │
 │  [0.28 × 0.36]       [0.28 × 0.36]         Dis-Arrest)       │
 │  Sim χ=0.9 CCDF      Sim χ=0.5 CCDF                          │
 └──────────────────────────────────────────────────────────────┘  0.11
```

Panels A–D are added first (auto-labeled A–D); the box plot is added last and is
auto-labeled **E**.

---

## Shared helper: `_plot_ccdf(ax, npy_file, title, signal_color)`

Loads a `.npy` file from `ev_data_dir`. File shape: `(2, N)` — row 0 = original eigenvalues, row 1 = scrambled.

| Parameter | Value | Notes |
|---|---|---|
| Scale | loglog | |
| Threshold `x2` | `max(scrambled)` | Dashed vertical line |
| Threshold label | `λ_max^scr` at axes y=0.8, fontsize `fsize-2` | Uses `ax.get_xaxis_transform()` for data-x / axes-y positioning |
| Spurious line | darkgray, α=0.7, markersize=3 | eigenvalues < x2 |
| Signal line | `signal_color`, markersize=3 | eigenvalues ≥ x2; color varies by panel condition |
| Scrambled line | black, α=0.5, markersize=3 | |
| CCDF formula | `1 - rank/p + 1/p` | Log-safe; computed independently for data and scrambled |
| x limits | `[0.1, 30]` | Fixed range shared across all CCDF panels |
| Legend | individual per-panel | spurious / signal / scrambled |

**To swap dataset:** change the filename argument.
**To change signal color:** pass `signal_color=` — `REG_COLOR` for Regulated/high-χ, `DIS_COLOR` for Dis-Arrest/low-χ.
**To change x-range:** edit `ax.set_xlim([0.1, 30])` in `_plot_ccdf` — applies to all CCDF panels.

---

## Panel A — Reg-Arrest experimental data

**Function:** `panel_A(ax)`
**Data file:** `ev_data/sample_15b_filtered.npy`
**Title:** `'Reg-Arrest'` · **Signal color:** `REG_COLOR` (steelblue)

## Panel B — Dis-Arrest experimental data

**Function:** `panel_B(ax)`
**Data file:** `ev_data/sample_15a_filtered.npy`
**Title:** `'Dis-Arrest'` · **Signal color:** `DIS_COLOR` (#E07B54 coral)

## Panel C — Simulated data, high χ

**Function:** `panel_C(ax)`
**Data file:** `ev_data/simulated_pcs_0.9.npy` (χ = 0.9)
**Title:** `r'Simulation ($\chi=0.9$)'` · **Signal color:** `REG_COLOR`

## Panel D — Simulated data, low χ

**Function:** `panel_D(ax)`
**Data file:** `ev_data/simulated_pcs_0.5.npy` (χ = 0.5)
**Title:** `r'Simulation ($\chi=0.5$)'` · **Signal color:** `DIS_COLOR`

---

## Panel E — Box plot: Regulated vs. Dis-Arrest

**Function:** `panel_E(ax)` (formerly `panel_F`)
**Data file:** `results/data_metrics/data_metrics.csv` (index_col=0)

| Parameter | Value | Notes |
|---|---|---|
| Metric | `sum_denoised_ev` column | GMP-Cor proxy |
| Excluded rows | `adam_matrix_filtered.csv`, `deb_Ec_CDS_untreated.csv`, `deb_KP_CDS_untreated.csv` | Dropped before plotting |
| group1 | `category == 'r'` | Regulated samples |
| group0 | `category == 'd'` | Dis-Arrest samples |
| Box edge color | `REG_COLOR` / `DIS_COLOR` | Condition-coded; facecolor α=0.12 tint |
| Strip plot | `s=10`, α=0.5, jitter ±0.12, seed=42 | Individual samples overlaid on boxes |
| Means | solid line | `meanline=True`, `showmeans=True` |
| Reference lines | `med_reg` (REG_COLOR dashed), `med_dis` (DIS_COLOR dashed) | Group medians from `test8.csv` |
| Significance | Mann-Whitney U (`stats.mannwhitneyu`) | Annotated `*`/`**`/`***`/`****`/`NS` via `format_p` |
| y limits | `bottom=-1, top=45` | |
| y-label | `'GMP-Cor'` | Standalone (no longer shares a y-axis) |
| x-tick labels | `['Regulated', 'Dis-Arrest']` | |

**`format_p(p)` thresholds:** p < 0.0001 → `****`, p < 0.001 → `***`, p < 0.01 → `**`, p < 0.05 → `*`, else → `NS`

---

## Data files used in figure

| File | Used in | Description |
|---|---|---|
| `ev_data/sample_15b_filtered.npy` | Panel A | Reg-Arrest experimental |
| `ev_data/sample_15a_filtered.npy` | Panel B | Dis-Arrest experimental |
| `ev_data/simulated_pcs_0.9.npy` | Panel C | Simulation χ = 0.9 |
| `ev_data/simulated_pcs_0.5.npy` | Panel D | Simulation χ = 0.5 |
| `results/data_metrics/data_metrics.csv` | Panel E | Dataset metrics with categories |
| `results/data_metrics/test8.csv` | Panel E (via `_load_group_medians`) | Group medians for reference lines |

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize the whole figure | `figsize=(W, H)` in `PanelFigure(...)` |
| Change all font sizes | `fsize = 10` at top of script |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Move the box plot relative to A–D | edit `panel_pos[4]` (the last entry) |
| Move panel label | `label_offset=(-0.04, 0.02)` in `PanelFigure(...)` |
| Swap CCDF dataset | change the filename string in the panel function |
| Change condition colors | `REG_COLOR` / `DIS_COLOR` at top of script — propagates to all panels |
| Use different simulation χ (CCDF) | replace `simulated_pcs_0.9.npy` or `simulated_pcs_0.5.npy` with another value |
| Fix x-range across all CCDF panels | edit `ax.set_xlim([0.1, 30])` in `_plot_ccdf` |
| Edit the calibration curve | it moved to `scripts/supplementary_figures/figure_s10.py` (panel A) |
| Save as SVG | `pf.save("figure3.svg", dpi=300, transparent=True)` at bottom of script |

---

## Permutation p-value indicator

Every CCDF panel carries the empirical permutation p-value that the observed signal
exceeds the scrambled null, as the last entry inside the legend box:

    p = (1 + #{lambda_1^perm >= lambda_1^obs}) / (B + 1),   B = 2000

Rendered `p < 5x10^-4` when no permutation reached the observed lambda_1 (the value
is censored at the 1/(B+1) resolution floor and is not resolved further), otherwise
`p = 0.006` style.

**Source:** `results/data_metrics/data_metrics.csv`, column `permutation_p`, written by
`scripts/add_permutation_metrics.py` from the B=2000 run
(`scripts/eigenvalue_permutation_full_B2000.py`).

**Implementation:** `scripts/figures/permutation_pvalues.py`, shared by figure2/3/5 and
figure_s5. `pv.legend_with_p(ax, npy_file, fontsize=...)` replaces the plain
`ax.legend(...)` call. The entry uses a zero-width legend handler (`_ZeroWidthHandle`) so
the text sits flush with the left edge of the box rather than indented into the label
column, and its font is set two points smaller than the series labels.

Datasets with no permutation test (the simulated spectra, `simulated_pcs_*.npy`) get an
ordinary three-entry legend -- the helper is a no-op for them.

| Goal | What to change |
|---|---|
| Reword or reformat the p-value | `p_label()` in `scripts/figures/permutation_pvalues.py` |
| Change its font size | `p_fontsize=` argument of `pv.legend_with_p` (default: legend size - 2) |
| Refresh the values | rerun `scripts/add_permutation_metrics.py` |

---
