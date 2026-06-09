# Figure 5 — Implementation Details

## How to run

Run from `scripts/figures/` so that `os.getcwd()` resolves correctly:

```bash
cd scripts/figures
python figure5.py          # interactive window
# uncomment pf.save(...) at the bottom to write SVG
```

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | Line ~308 — base font size (titles = 10pt, labels/ticks/legends = 8pt) |
| `figsize` | `(7, 6)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(0, 0.03)` | `PanelFigure(label_offset=...)` — shifts A/B/C… labels above panel top-left |
| `root_dir` | 2 dirs above `os.getcwd()` | Automatically resolves to repo root when running from `scripts/figures/` |
| `ev_data_dir` | `root_dir/ev_data/` | All `.npy` eigenvalue files live here |
| `REG_COLOR` | `'steelblue'` | Condition color for Regulated/high-ρ — used in panel D |
| `DIS_COLOR` | `'#E07B54'` | Condition color for Dis-Arrest/low-ρ — used in panel E |

---

## Global font standard

| Element | Size |
|---|---|
| Panel titles | `fsize` = 10 |
| Axis labels, tick labels, legends | `fsize - 2` = 8 |

---

## Layout map

Panel positions are `[left, bottom, width, height]` in **figure-normalized coordinates** (0–1).  
The figure is 7 × 6 inches.

```
 0.0                                                              1.0
 ┌──────────────────────────────────────────────────────────────────┐ ~0.94
 │  A (0.075, 0.60)  [0.50 × 0.34]                                 │
 │  BioRender schematic (VapC experiment timeline)                  │
 ├────────────────────────────────────────────────────────────────── ~0.65
 │  B (0.075, 0.45)            C (0.40, 0.45)                       │
 │  [0.24 × 0.20]              [0.24 × 0.20]                        │
 │  UMAP by time point         UMAP by cluster                      │
 ├────────────────────────────────────────────────────────────────── ~0.36
 │  D (0.075, 0.08)            E (0.40, 0.08)                       │
 │  [0.24 × 0.28]              [0.24 × 0.28]                        │
 │  Early VapC CCDF (2h)       Late VapC CCDF (24h)                 │
 └──────────────────────────────────────────────────────────────────┘ 0.08
                   Right column (x=0.70, width=0.275):
 ┌──────────────────────────────────────────────────────────────────┐ ~0.94
 │  F (0.70, 0.72)  [0.275 × 0.22]   GMP-Cor bar plot              │
 ├────────────────────────────────────────────────────────────────── ~0.63
 │  G (0.70, 0.41)  [0.275 × 0.22]   SDS sensitivity bar plot      │
 ├────────────────────────────────────────────────────────────────── ~0.30
 │  H (0.70, 0.08)  [0.275 × 0.22]   Lag time histograms           │
 └──────────────────────────────────────────────────────────────────┘ 0.08
```

---

## Shared helper: `_load_svg_image(svg_path)`

Parses the BioRender SVG file and extracts the embedded PNG (stored as a base64 data URI in the `xlink:href` attribute of the `<image>` element). Returns a PIL `Image` object.

| Detail | Value |
|---|---|
| Input file | `scripts/figures/figure5/biorender2.svg` |
| Encoding | PNG base64 data URI (`data:image/png;base64,...`) |
| Return type | `PIL.Image` |
| Output format | Raster (PNG content); final figure container stays SVG |

**Note**: The BioRender SVG already wraps a raster PNG — it is not a true vector image internally. Displaying via `ax.imshow` is equivalent to displaying the original SVG in terms of fidelity.

---

## Shared helper: `_plot_ccdf(ax, npy_file, title, signal_color)`

Loads a `.npy` file from `ev_data_dir`. File shape: `(2, N)` — row 0 = original eigenvalues, row 1 = scrambled.

| Parameter | Value | Notes |
|---|---|---|
| Scale | loglog | |
| Threshold `x2` | `max(scrambled)` | Dashed vertical line |
| Threshold label | `λ_max^scr` at axes y=0.8, fontsize `fsize-2` | Uses `ax.get_xaxis_transform()` for data-x / axes-y positioning |
| Noise line | darkgray, α=0.7, markersize=3 | eigenvalues < x2 |
| Signal line | `signal_color`, markersize=3 | eigenvalues ≥ x2; color varies by condition |
| Scrambled line | black, α=0.5, markersize=3 | |
| CCDF formula | `1 - rank/p + 1/p` | Log-safe |
| x limits | `[0.1, 30]` | Fixed range |
| Legend | individual per-panel | noise / signal / scrambled |

**To change x-range:** edit `ax.set_xlim([0.1, 30])` in `_plot_ccdf`.

---

## Panel A — BioRender schematic

**Function:** `panel_A(ax)`  
**Source file:** `scripts/figures/figure5/biorender2.svg`  
**Display:** `ax.imshow(img)` + `ax.axis('off')`

Shows the experimental timeline schematic: Exponential → VapC induction at T=0 → scRNA-seq sampling at 2h, 5h, 24h.

---

## Panel B — UMAP by time point

**Function:** `panel_B(ax)`  
**Data file:** `scanpy/umap_coordinates_vapc.csv`

| Batch label | Condition | Color |
|---|---|---|
| `exp` | Exponential | `#4393c3` |
| `T2` | VapC 2h | `#f4a582` |
| `T5A` / `T5B` | VapC 5h | `#d6604d` |
| `TON` | VapC 24h | `#b2182b` |

Legend placed at `lower right`, `bbox_to_anchor=(1.2, -0.2)`, `markerscale=4`, no frame.

---

## Panel C — UMAP by cluster

**Function:** `panel_C(ax)`  
**Data file:** `scanpy/umap_coordinates_vapc.csv`

Colors 5 clusters by `cluster` column (0–4). Cluster index annotated as text at centroid. Colors: `['#8073ac', '#b2182b', '#4393c3', '#d6604d', '#92c5de']`.

---

## Panel D — Early VapC CCDF (2h)

**Function:** `panel_D(ax)`  
**Data file:** `ev_data/VapC_biorep_t2A_filtered.npy`  
**Title:** `'Early VapC (2h)'`  
**Signal color:** `REG_COLOR` (steelblue) — this sample is regulated/high-ρ

---

## Panel E — Late VapC CCDF (24h)

**Function:** `panel_E(ax)`  
**Data file:** `ev_data/VapC_biorep_tONA_filtered.npy`  
**Title:** `'Late VapC (24h)'`  
**Signal color:** `DIS_COLOR` (#E07B54 coral) — this sample is dis-arrested/low-ρ

---

## Panel F — GMP-Cor bar plot (VapC time course)

**Function:** `panel_F(ax)`  
**Data file:** `results/data_metrics/test8.csv` (index_col=0)  
**Metric:** `sum_denoised_ev` column

| `file_name` in test8.csv | Label | Color | Value |
|---|---|---|---|
| `Expira_biorep_t0A_filtered.csv` | Exponential | `#4393c3` | ~31.9 |
| `VapC_biorep_t2A_filtered.csv` | VapC 2h | `#f4a582` | ~27.5 |
| `VapC_biorep_t5A_filtered.csv` | VapC 5h | `#d6604d` | ~13.2 |
| `VapC_biorep_tONA_filtered.csv` | VapC 24h | `#b2182b` | ~4.9 |

| Parameter | Value |
|---|---|
| Bar width | `0.25` |
| Gap between bars | `0.4` |
| y limits | `[0, 45]` |
| y label | `'GMP-Cor'` |
| Error bars | none (single samples) |

Samples are matched from test8.csv by the `file_name` column. The order (Exp → 2h → 5h → 24h) reflects increasing dysregulation over VapC induction time.

---

## Panel G — SDS sensitivity bar plot

**Function:** `panel_G(ax)`  
**Data file:** `scripts/figures/figure5/normalizedOD_at_20h.csv`

5 conditions: Reg 2h, Reg 24h, VapC 2h, VapC 5h, VapC 24h. Y-axis: Normalized OD (log scale). Error bars: SE (biological vs technical replicate rule). Significance annotations via Welch t-test between adjacent bars.

---

## Panel H — Lag time distribution

**Function:** `panel_H(ax)`  
**Data files:** `scripts/figures/figure5/CTRLt0.csv`, `VAPCt240.csv`, `VAPCt1400.csv`

Overlapping histograms of single-cell lag times for 3 conditions (Reg-Arrest, Early VapC, Late VapC). X: 300–1100 min. Y: frequency density (×10⁻² scaling label). Colors: `['#9ecae1', '#fb6a4a', '#a50f15']`.

---

## Data files used in figure

| File | Used in | Description |
|---|---|---|
| `scripts/figures/figure5/biorender2.svg` | Panel A | BioRender experiment schematic |
| `scanpy/umap_coordinates_vapc.csv` | Panels B, C | UMAP coordinates + batch/cluster labels |
| `ev_data/VapC_biorep_t2A_filtered.npy` | Panel D | Eigenvalues: Early VapC (2h) |
| `ev_data/VapC_biorep_tONA_filtered.npy` | Panel E | Eigenvalues: Late VapC (24h) |
| `results/data_metrics/test8.csv` | Panel F | GMP-Cor (`sum_denoised_ev`) per sample |
| `scripts/figures/figure5/normalizedOD_at_20h.csv` | Panel G | SDS growth assay OD data |
| `scripts/figures/figure5/CTRLt0.csv` | Panel H | Control lag times |
| `scripts/figures/figure5/VAPCt240.csv` | Panel H | Early VapC lag times |
| `scripts/figures/figure5/VAPCt1400.csv` | Panel H | Late VapC lag times |

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize the whole figure | `figsize=(W, H)` in `PanelFigure(...)` |
| Change all font sizes | `fsize = 10` in the assembly block |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Move panel label | `label_offset=(0, 0.03)` in `PanelFigure(...)` |
| Swap CCDF dataset | change the filename string in `panel_D` or `panel_E` |
| Change condition colors | `REG_COLOR` / `DIS_COLOR` in the assembly block |
| Change CCDF x-range | edit `ax.set_xlim([0.1, 30])` in `_plot_ccdf` |
| Move threshold label vertically | change `0.8` in `ax.text(x2 * 1.1, 0.8, ...)` — axes coords |
| Change GMP-Cor y-range | edit `ax.set_ylim([0, 45])` in `panel_F` |
| Add error bars to panel F | add `yerr=` to `ax.bar(...)` and supply per-sample SD/SE |
| Save as SVG | uncomment `pf.save("figure5.svg", dpi=300)` at bottom |
