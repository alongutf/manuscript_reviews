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
| `fsize` | `10` | Line 16 — base font size (all labels/ticks are `fsize`, `fsize-2`, `fsize-3`, etc.) |
| `figsize` | `(7, 9.5)` inches | `PanelFigure(figsize=...)` at line 292 |
| `label_offset` | `(-0.04, 0.04)` | `PanelFigure(label_offset=...)` — shifts A/B/C… labels relative to panel top-left corner |
| `root_dir` | 2 dirs above `os.getcwd()` | Automatically resolves to repo root when running from `scripts/figures/` |
| `ev_data_dir` | `root_dir/ev_data/` | Line 19 — all `.npy` eigenvalue files live here |

---

## Layout map

Panel positions are `[left, bottom, width, height]` in **figure-normalized coordinates** (0–1).  
The figure is 7 × 9.5 inches; 1 unit ≈ 7 in wide, 9.5 in tall.

```
 0.0                                                        1.0
 ┌────────────────────────────────────────────────────────────┐ 0.97
 │  A (0.04,0.61)     B (0.27,0.71)   C (0.52,0.61)          │
 │  [0.18 × 0.35]     [0.21 × 0.23]   [0.44 × 0.35]          │
 │  2-row grid        single          1×2 grid                │
 │  ─ heatmap         GMP curves      ─ orig sparse           │
 │  ─ MP hist                         ─ scrambled             │
 ├──────────────────────────────────────────────────────────── 0.61
 │                                                            │
 │  D (0.04,0.05)     E (0.36,0.30)                           │
 │  [0.26 × 0.50]     [0.59 × 0.24]                           │
 │  2-row grid        single (wide)                           │
 │  ─ sim pcs         PDF + inset CCDF                        │
 │  ─ sim scrambled                                           │
 │                    F (0.36,0.05)   G (0.69,0.05)           │
 │                    [0.28 × 0.21]   [0.27 × 0.21]           │
 │                    CCDF only       CCDF only               │
 └────────────────────────────────────────────────────────────┘ 0.05
```

To shift a panel: edit `panel_pos[i]` in the assembly block (lines 294–302).  
To resize: change `width` or `height` (4th value).  
Keep ≥ 0.03 gap between adjacent panels to avoid overlap.

---

## Panel A — Random matrix intro

**Function:** `panel_A(axes)` — called with a 2×1 grid  
**Grid call:** `pf.add_grid_panel(panel_pos[0], 2, 1, hspace=0.4)`

| Parameter | Value | Notes |
|---|---|---|
| Matrix size | N=500 rows, P=1000 cols | γ = P/N = 2 |
| Heatmap subset | first 20×20 of correlation matrix | Change `:20, :20` for more/fewer entries |
| Heatmap colormap | `'bwr'`, center=0, vmin/vmax ±0.5 | |
| Colorbar ticks | `[-0.5, 0, 0.5]` | |
| Histogram bins | `np.linspace(0, 6, 40)` | Adjust range/count to change resolution |
| MP curve | `af.mp_distribution(val, P/N)` | x from 0 to 6 in 100 steps |
| Colors | histogram `'red'` α=0.5, MP line `'red'` | |
| hspace | 0.4 | vertical gap between the two sub-rows |

**To change:** increase N and P together (keep γ=P/N=2) for a smoother histogram. Change `hspace` to push the two sub-panels apart.

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
| Title | `'Generalized MP\n(with correlations)'` | |

---

## Panel C — Simulation data structure

**Function:** `panel_C(axes)` — called with a 1×2 grid  
**Grid call:** `pf.add_grid_panel(panel_pos[2], 1, 2, wspace=0.12)`

### Synthetic matrix
| Parameter | Value | Notes |
|---|---|---|
| Random seed | `np.random.seed(42)` | Change seed to get a different illustrative pattern |
| Matrix size | 22 cells × 14 genes | Adjust `n_cells`, `n_genes` |
| Sparsity | 15–38% non-zero per column (`np.random.uniform(0.15, 0.38)`) | Increase range for more/less sparsity |
| Non-zero values | `np.random.exponential(2.5)` | Adjust scale for brighter/dimmer non-zero entries |
| Colormap | `'YlOrRd'`, vmin=0, vmax=matrix.max() | Both panels share the same vmax |

### Arrows
| Feature | Where | Key parameters |
|---|---|---|
| **Within-panel arc arrows** (in scrambled subplot) | `axes[0,1].annotate(...)` | Columns targeted: `[2, 6, 11]`; threshold: value > 0.5; min movement: 4 rows; arc radius: `rad=0.4`; color: `#2171b5` |
| **Cross-panel schematic arrow** | `fig.add_artist(FancyArrowPatch(...))` | Uses `axes.get_position()` to compute x-coordinates; `mutation_scale=10`; color `#636363` |
| "shuffle" label | `fig.text(...)` | y offset +0.007 above arrow midpoint; `fsize-3` |

**To change arrows:** edit columns list `[2, 6, 11]` or lower the `> 0.5` threshold to show more arrows. Increase `rad` in `arc3,rad=X` for wider arcs. Change `lw` for thicker arrows.

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
| Color | Condition |
|---|---|
| `'darkgray'` | bin_x < x1 (within MP noise) |
| `'salmon'` | x1 ≤ bin_x < x2 (above MP but below scrambled max) |
| `'skyblue'` | bin_x ≥ x2 (signal: above scrambled threshold) |

### Other settings
| Parameter | Value | Notes |
|---|---|---|
| bin_width | 0.2 | Shared by both subplots |
| bar width | `bin_width * 0.8` | Slightly narrower than bin for visual separation |
| align | `'right'` | Bars anchor to right edge of each bin |
| x-axis shared | Manually via `set_xlim(xlim)` on both subplots | xlim derived from combined data range |
| Threshold lines | `x1`: black `'--'`; `x2`: dimgray `':'` | |
| yticks | [0, 0.1, 0.2, 0.3, 0.4] | Both subplots |
| Top xticks | hidden (`set_xticks([])`) | Only bottom subplot shows λ axis |

**To change:** swap `simulated_pcs.npy` for `simulated_pcs_0.9.npy` etc. to show a specific correlation strength. Change `x1` formula if γ ≠ 2. Adjust `bin_width` for finer/coarser histogram.

---

## Panel E — Regulated dataset: PDF + inset CCDF

**Function:** `panel_E(ax)` — single wide panel  
**Panel call:** `pf.add_panel(panel_pos[4], draw_func=panel_E)`  
**Data file:** `ev_data/sample_13b_filtered.npy` — shape (2, N): row 0 = data, row 1 = scrambled

### Main histogram
| Parameter | Value | Notes |
|---|---|---|
| bin_width | 0.2 | |
| bin_edges upper cap | 10 | `np.arange(..., 10 + bin_width, ...)` — change to extend x-axis |
| Data bars | width=`bin_width*0.5`, `align='left'` | Left half of each bin |
| Scrambled bars | shifted by `+bin_width*0.5`, `align='right'` | Right half of each bin, solid black |
| Color threshold | `x2 = max(scrambled)` | Bars < x2: darkgray; bars ≥ x2: skyblue |
| yticks | [0, 0.1, 0.2, 0.3, 0.4] | |

### Inset CCDF
| Parameter | Value | Notes |
|---|---|---|
| Inset position | `[0.40, 0.38, 0.56, 0.56]` | `[left, bottom, width, height]` in axes-fraction coords |
| Scale | loglog | |
| Noise line | darkgray, α=0.7 | eigenvalues < x2 |
| Signal line | skyblue | eigenvalues ≥ x2 |
| Scrambled line | black, α=0.5 | |
| CCDF formula | `1 - rank/p + 1/p` | Avoids zero (log-safe). Computed independently for data and scrambled |
| x limits | `[0.1, max(data)]` | Change 0.1 to extend left tail |
| Threshold vline | x2, black dashed | |

**To change dataset:** replace `'sample_13b_filtered.npy'` with any `.npy` in `ev_data/`. Available regulated files: `sample_13b`, `sample_13a`, `sample_15a`, `sample_15b`. Update the title string accordingly.

---

## Panels F & G — CCDF-only panels

Both call the shared helper `_plot_ccdf(ax, npy_file, title)`.

| Panel | Data file | Title |
|---|---|---|
| F | `sample_2b_filtered.npy` | `'Exponential'` |
| G | `sample_15b_filtered.npy` | `'Reg-Arrest (rep. 2)'` |

`_plot_ccdf` is identical to panel E's inset CCDF, but full-panel. Same color scheme (darkgray noise / skyblue signal / black scrambled), same CCDF formula.

**To swap dataset:** change the filename string in `panel_F` or `panel_G`.  
**To add a legend title or x-limits:** edit `_plot_ccdf` directly — it applies to both F and G.

---

## Available eigenvalue data files (`ev_data/`)

| File | Condition | GMP-Cor (approx.) |
|---|---|---|
| `simulated_pcs.npy` | Synthetic simulation | — |
| `simulated_pcs_0.npy` – `simulated_pcs_1.npy` | Simulation at χ = 0 … 1 | — |
| `sample_2b_filtered.npy` | Exponential *E. coli* (our data) | ~0 |
| `sample_13a_filtered.npy` | Dis-Arrest rep. A | ~1.8 |
| `sample_13b_filtered.npy` | Reg-Arrest rep. B | ~7.9 |
| `sample_15a_filtered.npy` | Dis-Arrest rep. A | ~0 |
| `sample_15b_filtered.npy` | Reg-Arrest rep. B | ~11.2 |

All files: shape `(2, N)` — row 0 = original eigenvalues, row 1 = scrambled eigenvalues.

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize the whole figure | `figsize=(W, H)` in `PanelFigure(...)` |
| Change all font sizes at once | `fsize = 10` at line 16 |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Change panel E to a different sample | replace `'sample_13b_filtered.npy'` and update the title string |
| Show more/fewer GMP curves in B | edit the `files` and `labels` lists in `panel_B` |
| Change simulation correlation strength in D | replace `'simulated_pcs.npy'` with `'simulated_pcs_0.9.npy'` etc. |
| Adjust inset position in E | `ax.inset_axes([left, bottom, width, height])` in `panel_E` |
| Make C arrows more visible | lower threshold `> 0.5`, increase `lw`, or change `#2171b5` to a brighter color |
| Save as SVG | uncomment `pf.save("figure2.svg", dpi=300)` and comment out `plt.show()` |
