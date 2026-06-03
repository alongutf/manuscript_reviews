# Figure 2 — Implementation Details

## How to run

Run from `scripts/figures/` so that `os.getcwd()` resolves correctly:

```bash
cd scripts/figures
python figure2.py          # interactive window
# or save directly:
# uncomment pf.save("figure2.svg", dpi=300) at the bottom
```

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | Line 16 — base font size |
| `figsize` | `(7, 7.5)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(-0.04, 0.04)` | `PanelFigure(label_offset=...)` — shifts A/B/C… labels relative to panel top-left corner |
| `root_dir` | 2 dirs above `os.getcwd()` | Automatically resolves to repo root when running from `scripts/figures/` |
| `ev_data_dir` | `root_dir/ev_data/` | All `.npy` eigenvalue files live here |

---

## Global font standard

| Element | Size |
|---|---|
| Panel titles | `fsize` = 10 |
| Axis labels, tick labels, legends | `fsize - 2` = 8 |

---

## Layout map

Panel positions are `[left, bottom, width, height]` in **figure-normalized coordinates** (0–1).  
The figure is 7 × 7.5 inches; 1 unit ≈ 7 in wide, 7.5 in tall.

```
 0.0                                                        1.0
 ┌────────────────────────────────────────────────────────────┐ 0.97
 │  A (0.08,0.80)     C (0.32,0.54)   D (0.56,0.54)          │
 │  [0.19 × 0.15]     [0.15 × 0.40]   [0.40 × 0.40]          │
 │  single histogram  2×1 grid        2-row grid              │
 │                    ─ orig sparse   ─ sim eigenvalues       │
 │  B (0.08,0.50)     ─ scrambled     ─ sim scrambled         │
 │  [0.19 × 0.20]                                             │
 │  single GMP curves                                         │
 ├──────────────────────────────────────────────────────────── ~0.47
 │                                                            │
 │  E (0.08,0.07)                     F (0.66,0.30)           │
 │  [0.47 × 0.35]                     [0.30 × 0.15]           │
 │  single (wide)                     CCDF only               │
 │  PDF + inset CCDF                                          │
 │                                    G (0.66,0.06)           │
 │                                    [0.30 × 0.15]           │
 │                                    CCDF only               │
 └────────────────────────────────────────────────────────────┘ 0.06
```

To shift a panel: edit `panel_pos[i]` in the assembly block.  
To resize: change `width` or `height` (4th value).  
Keep ≥ 0.03 gap between adjacent panels to avoid overlap.

---

## Panel A — Random matrix intro

**Function:** `panel_A(ax)` — single panel  
**Panel call:** `pf.add_panel(panel_pos[0], draw_func=panel_A)`

| Parameter | Value | Notes |
|---|---|---|
| Matrix size | N=500 rows, P=1000 cols | γ = P/N = 2 |
| Histogram bins | `np.linspace(0, 6, 40)` | Adjust range/count to change resolution |
| MP curve | `af.mp_distribution(val, P/N)` | x from 0 to 6 in 100 steps; plotted with `linestyle='--'` |
| Colors | histogram `'r'` α=0.5, MP line `'r'` | |
| Title | `'Random matrix (RM)\nMP distribution'` | |

**To change:** increase N and P together (keep γ=P/N=2) for a smoother histogram.

---

## Panel B — Generalized MP

**Function:** `panel_B(ax)` — single panel  
**Panel call:** `pf.add_panel(panel_pos[1], draw_func=panel_B)`

| Parameter | Value | Notes |
|---|---|---|
| Data files | `model_alpha2_sigma0/07/08/09.txt` | In `root_dir/model fit/`, col 0 = λ, col 1 = ρ(λ) |
| χ values plotted | 0, 0.7, 0.8, 0.9 | Change list `files` + `labels` |
| Colors | `plt.cm.RdBu` at [0.1, 0.7, 0.85, 0.95] | Adjust the index list for different shades |
| x limits | 0 – 8.5 | `ax.set_xlim` |
| y limits | 0 – 0.35 | `ax.set_ylim` |
| x ticks | [0, 2, 4, 6, 8] | |
| y ticks | [0, 0.1, 0.2, 0.3] | |
| Title | `'RM with correlations\nGeneralized MP'` | |

---

## Panel C — Simulation data structure

**Function:** `panel_C(axes)` — called with a 2×1 grid  
**Grid call:** `pf.add_grid_panel(panel_pos[2], 2, 1, hspace=0.30)`  
`axes[0, 0]` = original matrix (top); `axes[1, 0]` = scrambled matrix (bottom)

### Synthetic matrix
| Parameter | Value | Notes |
|---|---|---|
| Random seed | `np.random.seed(42)` | Change seed to get a different illustrative pattern |
| Matrix size | 22 cells × 14 genes | Adjust `n_cells`, `n_genes` |
| Sparsity | 15–38% non-zero per column (`np.random.uniform(0.15, 0.38)`) | Increase range for more/less sparsity |
| Non-zero values | `np.random.exponential(2.5)` | Adjust scale for brighter/dimmer non-zero entries |
| Colormap | `'Greys'`, vmin=0, vmax=matrix.max() | Both panels share the same vmax |
| Top title | `'Synthetic data\nwith correlations'` | |
| Scrambled title | none (commented out) | Bottom panel has ylabel `'cells'` but no title |

No arrows are drawn. Change `hspace` in the grid call to widen the gap between the two heatmaps.

---

## Panel D — Simulation eigenvalue distributions

**Function:** `panel_D(axes)` — called with a 2×1 grid  
**Grid call:** `pf.add_grid_panel(panel_pos[3], 2, 1, hspace=0.3)`  
**Data file:** `ev_data/simulated_pcs.npy` — shape (2, N): row 0 = pcs, row 1 = pcs1 (scrambled)

### Thresholds
| Threshold | Formula | Meaning |
|---|---|---|
| `x1` | `(1 + √2)²` ≈ 5.83 | Analytical MP upper edge (γ = P/N = 2) |
| `x2` | `max(pcs1)` | Maximum scrambled eigenvalue = GMP-Cor threshold |

### Color coding (top subplot bars)
| Color | Condition | Legend label |
|---|---|---|
| `'darkgray'` | bin_x < x1 | `'MP noise'` |
| `'salmon'` | x1 ≤ bin_x < x2 | `'sparsity induced correlations'` |
| `'skyblue'` | bin_x ≥ x2 | `'true correlations'` |

### Other settings
| Parameter | Value | Notes |
|---|---|---|
| bin_width | 0.15 | Shared by both subplots |
| bar width | `bin_width * 0.8` | Slightly narrower than bin for visual separation |
| align | `'right'` | Bars anchor to right edge of each bin |
| x-axis | fixed `[2, 12]` on both subplots | |
| y-axis | fixed `[0, 0.2]` on both subplots | |
| Threshold lines | top: `x1` black `'--'`; both: `x2` dimgray `':'` | Bottom vline at `x2 + bin_width` |
| yticks | [0, 0.1, 0.2, 0.3, 0.4] | Both subplots (clipped by ylim) |
| Top title | `'Correlation eigenvalue density'` | |
| Bottom title | none (empty string) | |

**To change:** swap `simulated_pcs.npy` for `simulated_pcs_0.9.npy` etc. to show a specific correlation strength. Adjust `bin_width` for finer/coarser histogram.

---

## Panel E — Regulated dataset: PDF + inset CCDF

**Function:** `panel_E(ax)` — single wide panel  
**Panel call:** `pf.add_panel(panel_pos[4], draw_func=panel_E)`  
**Data file:** `ev_data/Expira_biorep_t0A_filtered.npy` — shape (2, N): row 0 = data, row 1 = scrambled  
**Title:** `r'Exponential $\it{E. coli}$'`

### Main histogram
| Parameter | Value | Notes |
|---|---|---|
| bin_width | 0.4 | |
| bin_edges upper cap | 12 | `np.arange(..., 12 + bin_width, ...)` — change to extend x-axis |
| Data bars | width=`bin_width*0.5`, `align='left'` | Left half of each bin |
| Scrambled bars | shifted by `+bin_width*0.5`, `align='right'` | Right half of each bin, solid black |
| Color threshold | `x2 = max(scrambled)` | Bars < x2: darkgray; bars ≥ x2: skyblue |
| yticks | [0, 0.1, 0.2, 0.3, 0.4] | |

### Inset CCDF
| Parameter | Value | Notes |
|---|---|---|
| Inset position | `[0.40, 0.38, 0.56, 0.56]` | `[left, bottom, width, height]` in axes-fraction coords |
| Scale | loglog | |
| Noise line | darkgray, α=0.7, markersize=3 | eigenvalues < x2 |
| Signal line | skyblue, markersize=3 | eigenvalues ≥ x2 |
| Scrambled line | black, α=0.5, markersize=3 | |
| CCDF formula | `1 - rank/p + 1/p` | Avoids zero (log-safe). Computed independently for data and scrambled |
| x limits | `[0.1, max(data) * 1.5]` | |
| Threshold vline | x2, black dashed | |

**To change dataset:** replace `'Expira_biorep_t0A_filtered.npy'` with any `.npy` in `ev_data/`. Update the title string accordingly.

---

## Panels F & G — CCDF-only panels

Both call the shared helper `_plot_ccdf(ax, npy_file, title)`.

| Panel | Data file | Title |
|---|---|---|
| F | `deb_Ec_CDS_untreated.npy` | `"Exponential $\it{E. coli}$, Ma et. al"` |
| G | `deb_KP_CDS_untreated.npy` | `"Exponential $\it{K. pneumoniae}$, Ma et. al"` |

`_plot_ccdf` uses loglog scale, same color scheme as panel E's inset (darkgray noise / skyblue signal / black scrambled), same CCDF formula (`1 - rank/p + 1/p`), markersize=3, x limits `[0.1, max(data)*1.5]`. Note: title fontsize in `_plot_ccdf` is `fsize-2` (8pt) rather than `fsize`.

**To swap dataset:** change the filename string in `panel_F` or `panel_G`.  
**To change x-limits or markersize:** edit `_plot_ccdf` directly — it applies to both F and G.

---

## Data files used in figure (`ev_data/`)

All files: shape `(2, N)` — row 0 = original eigenvalues, row 1 = scrambled eigenvalues.

| File | Used in | Condition |
|---|---|---|
| `simulated_pcs.npy` | Panel D | Synthetic simulation |
| `simulated_pcs_0.npy` – `simulated_pcs_1.npy` | (swap for D) | Simulation at χ = 0 … 1 |
| `Expira_biorep_t0A_filtered.npy` | Panel E | Exponential *E. coli* (our data) |
| `deb_Ec_CDS_untreated.npy` | Panel F | Exponential *E. coli* (Ma et al.) |
| `deb_KP_CDS_untreated.npy` | Panel G | Exponential *K. pneumoniae* (Ma et al.) |

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize the whole figure | `figsize=(W, H)` in `PanelFigure(...)` in assembly block |
| Change all font sizes at once | `fsize = 10` at line 16 |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Change panel E dataset | replace `'Expira_biorep_t0A_filtered.npy'` and update the title string |
| Change panel F or G dataset | replace filename in `panel_F` / `panel_G` and update title |
| Show more/fewer GMP curves in B | edit the `files` and `labels` lists in `panel_B` |
| Change simulation correlation strength in D | replace `'simulated_pcs.npy'` with `'simulated_pcs_0.9.npy'` etc. |
| Adjust inset position in E | `ax.inset_axes([left, bottom, width, height])` in `panel_E` |
| Adjust CCDF x-range | change `np.max(d1s)*1.5` multiplier in `_plot_ccdf` or `panel_E` inset |
| Save as SVG | uncomment `pf.save("figure2.svg", dpi=300)` and comment out `plt.show()` |
