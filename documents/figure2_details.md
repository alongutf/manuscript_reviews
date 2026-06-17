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
 │  A (0.08,0.80)     C (0.33,0.50)   D (0.57,0.50)          │
 │  [0.19 × 0.15]     [0.15 × 0.45]   [0.39 × 0.45]          │
 │  single histogram  2×1 grid        2-row grid              │
 │                    ─ orig sparse   ─ sim eigenvalues+CCDF  │
 │  B (0.08,0.50)     ─ scrambled       inset (top)          │
 │  [0.19 × 0.20]                     ─ sim scrambled         │
 │  single GMP curves                                         │
 │  C & D span 0.50→0.95 = A+B combined height               │
 ├──────────────────────────────────────────────────────────── ~0.45
 │                                                            │
 │  E (0.08,0.09)     F (0.40,0.09)   G (0.72,0.09)          │
 │  [0.24 × 0.32]     [0.24 × 0.32]   [0.24 × 0.32]          │
 │  CCDF only         CCDF only       CCDF only              │
 │  ── single row of three CCDF panels ──                    │
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

Panel C now spans `0.50 → 0.95` in height (`[0.33, 0.50, 0.15, 0.45]`) so its total
height matches the combined height of panels A (`0.80→0.95`) and B (`0.50→0.70`).

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

Panel D now spans `0.50 → 0.95` in height (`[0.57, 0.50, 0.39, 0.45]`) so its total
height matches the combined height of panels A + B.

### Thresholds
| Threshold | Formula | Meaning |
|---|---|---|
| `x1` | `(1 + √2)²` ≈ 5.83 | Analytical MP upper edge (γ = P/N = 2) |
| `x2` | `max(pcs1)` | Maximum scrambled eigenvalue = GMP-Cor threshold |

### Color coding — **top subplot** (original pcs)
| Color | Condition | Category |
|---|---|---|
| `'darkgray'` | bin_x < x1 | MP spurious correlations |
| `'salmon'` | x1 ≤ bin_x < x2 | Sparsity induced spurious correlations |
| `'skyblue'` | bin_x ≥ x2 | true correlations |

### Color coding — **bottom subplot** (scrambled pcs1)
Same scheme **without** the `skyblue` class (scrambled data has no eigenvalue ≥ x2):
| Color | Condition |
|---|---|
| `'darkgray'` | bin_x < x1 |
| `'salmon'` | bin_x ≥ x1 (intermediate / sparsity induced) |

Coloring the intermediate (salmon) bars in the bottom plot shows that sparsity-induced
correlations appear even in scrambled data.

### Shared legend
A single legend (drawn on the **bottom** subplot, `loc='upper right'`, fontsize `fsize-3`)
serves both subplots, with three `Patch` handles:
`'MP spurious correlations'` (darkgray), `'Sparsity induced spurious correlations'`
(salmon), `'true correlations'` (skyblue).

### CCDF inset (top subplot)
`ax_top.inset_axes([0.46, 0.42, 0.52, 0.55])`, drawn via the shared `_draw_ccdf` helper
(`show_legend=False`, `markersize=2`). Uses the shared CCDF color scheme
(grey spurious / blue true / black scrambled). Inset font is `fsize-3`.

### Other settings
| Parameter | Value | Notes |
|---|---|---|
| bin_width | 0.15 | Shared by both subplots |
| bar width | `bin_width * 0.8` | Slightly narrower than bin for visual separation |
| align | `'right'` | Bars anchor to right edge of each bin |
| x-axis | fixed `[2, 12]` on both subplots | |
| y-axis | fixed `[0, 0.2]` on both subplots | |
| Threshold lines | **both subplots**: `x1` black `'--'`; `x2` dimgray `':'` | Extended to both subplots |
| yticks | [0, 0.1, 0.2, 0.3, 0.4] | Both subplots (clipped by ylim) |
| Top title | `'Correlation eigenvalue density'` | |
| Bottom title | none | |

**To change:** swap `simulated_pcs.npy` for `simulated_pcs_0.9.npy` etc. to show a specific correlation strength. Adjust `bin_width` for finer/coarser histogram.

---

## Panels E, F & G — CCDF-only panels (single row)

E, F, and G are now **three CCDF-only panels in a single bottom row**
(`[0.08/0.40/0.72, 0.09, 0.24, 0.32]`). Panel E used to be a wide PDF histogram
with a CCDF inset — that is gone; it is now just a CCDF plot like F and G.

All three call the shared helper `_plot_ccdf(ax, npy_file, title)`, which in turn
calls `_draw_ccdf(...)`.

| Panel | Data file | Title |
|---|---|---|
| E | `Expira_biorep_t0A_filtered.npy` | `r'Exponential $\it{E. coli}$'` (our data) |
| F | `deb_Ec_CDS_untreated.npy` | `"Exponential $\it{E. coli}$, Ma et. al"` |
| G | `deb_KP_CDS_untreated.npy` | `"Exponential $\it{K. pneumoniae}$, Ma et. al"` |

### `_draw_ccdf(ax, data1, data2, show_legend=True, markersize=3)`
Shared CCDF drawing routine (used by E/F/G **and** the panel D inset).

| Element | Color | Legend label |
|---|---|---|
| eigenvalues < x2 (`x2 = max(scrambled)`) | `'darkgray'`, α=0.7 | `'spurious correlations'` |
| eigenvalues ≥ x2 | `'skyblue'` | `'true correlations'` |
| scrambled data | `'black'`, α=0.5 | `'scrambled data'` |

- Scale: loglog. CCDF formula `1 - rank/p + 1/p` (log-safe, computed independently for data and scrambled).
- x limits `[0.1, max(data)*1.5]`; threshold `x2` drawn as a black dashed vline.
- `show_legend=False` suppresses the legend (used for the compact panel D inset).

`_plot_ccdf` adds the `λ` / `CCDF` axis labels and the title (fontsize `fsize-2`).

**To swap a dataset:** change the filename string in `panel_E` / `panel_F` / `panel_G`.  
**To change x-limits, markersize, colors, or legend labels:** edit `_draw_ccdf` — it applies to E, F, G and the panel D inset at once.

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
| Change panel E/F/G dataset | replace filename in `panel_E` / `panel_F` / `panel_G` and update title |
| Show more/fewer GMP curves in B | edit the `files` and `labels` lists in `panel_B` |
| Change simulation correlation strength in D | replace `'simulated_pcs.npy'` with `'simulated_pcs_0.9.npy'` etc. |
| Adjust panel D CCDF inset position | `ax_top.inset_axes([left, bottom, width, height])` in `panel_D` |
| Adjust CCDF x-range / colors / labels | edit `_draw_ccdf` — applies to E, F, G and the D inset |
| Save as SVG | uncomment `pf.save("figure2.svg", dpi=300)` and comment out `plt.show()` |
