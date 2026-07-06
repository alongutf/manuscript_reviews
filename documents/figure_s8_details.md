# Figure S8 — Implementation Details

Supplementary figure applying the **correlation-spectrum / GMP-Cor** framework (developed for the
bacterial datasets) to the **PC9 cancer cell line** (day-14 scRNA-seq). It demonstrates that the
same global-dysregulation readout transfers to a mammalian system.

**Important label convention:** `low`/`med`/`high` refer to the **fluorescent marker** level; the
proliferation state is the **opposite** — low marker = **proliferating** (cycling), high marker =
**arrested** (non-cycling). Only the proliferation names (`Proliferating`/`Intermediate`/`Arrested`)
are shown to readers; the marker level is not printed anywhere in the figure. The **proliferating**
fraction carries the strongest coordinated correlation structure (highest GMP-Cor); the **arrested**
fraction the weakest.

**Layout (2×2):** A Proliferating CCDF · B Arrested CCDF · C UMAP · D GMP-Cor bar plot.

## How to run

The script lives in `scripts/supplementary_figures/`.

```bash
cd scripts/supplementary_figures
python figure_s8.py        # writes figure_s8.pdf + figure_s8_preview.png next to the script
```

Import bootstrap (path-independent), identical to the other S-figures:

```python
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                # repo root
sys.path.insert(0, _REPO)                                      # -> import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))  # -> import figure_functions
```

---

## Inputs

| Input | Path | Used for |
|---|---|---|
| Cell metadata + UMAP | `pc9_data/obs_metadata.csv` | **Panel C.** Columns `UMAP_X`, `UMAP_Y`, `sample` (`high1/high2/med1/med2/low1/low2`), `leiden` (clusters `0,1,2,3`). 26,990 cells. |
| Eigenvalue spectra | `results/cancer_pc9/ev_data/<sample>.npy` | **Panels A & B.** One array per day-14 sample, shape `(2, P)`: row 0 = empirical eigenvalues, row 1 = scrambled. |
| GMP-Cor metrics | `results/cancer_pc9/pc9_day14_gmp_cor.csv` | **Panels A, B & D.** Column `GMP_Cor` per `sample` (`14_rep{1,2}_{high,med,low}`). |

All three ev_data / metrics inputs are produced by **`scripts/cancer_pc9_correlation_analysis.py`**
(the upstream pipeline: gene panel = union of top-400-scoring genes in marker groups 0 & 2 →
800 genes; cells = top n_genes/2 = 400 by total UMI per sample; `af.get_eig_dist` with
`norm='sum'`, `norm_sum=1`; GMP-Cor = `sum(max(λ − λ_max^scr, 0))`). Re-run that script first if the
underlying data changes.

---

## Global settings

| Variable | Value | Meaning |
|---|---|---|
| `fsize` | `10` | base font size (panel letters 10pt; titles 9pt; axis labels/ticks/annotations 8pt; CCDF legend 7pt) |
| `figsize` | `(7, 7.5)` in | 2×2 portrait layout |
| `COND_COLOR` | `{'low':'steelblue','med':'#fec44f','high':'salmon'}` | marker-level colors, shared across all panels |
| `COND_ORDER` | `['low','med','high']` | plotting / bar order (marker level) |
| `COND_DISPLAY` | `{'low':'Proliferating','med':'Intermediate','high':'Arrested'}` | marker level → proliferation display name (the only label shown to readers) |
| `CLUSTER_LABEL` | `{2:'Cycling', 0:'Non-cycling'}` | defined but **not currently drawn** (UMAP cluster annotations removed) |
| `REP_HIGH` / `REP_LOW` | `14_rep1_high` / `14_rep1_low` | the two CCDF samples (Arrested / Proliferating) |

The fluorescent-marker level (`low/med/high`) is **not** shown to readers anywhere — only the
proliferation display names (`Proliferating`/`Intermediate`/`Arrested`) appear in legends, titles,
and tick labels.

`_condition(tag)` maps any sample tag (`high1`, `14_rep1_high`, …) to `low`/`med`/`high` by substring;
`COND_DISPLAY` then converts that to the proliferation label shown in every panel.

---

## Layout

Four `add_panel` calls (normalized figure coords `[left, bottom, width, height]`), a 2×2 grid:

```python
pf = PanelFigure(figsize=(7, 7.5), label_offset=(-0.05, 0.015))
# Row 1 — representative CCDFs (Proliferating | Arrested)
pf.add_panel([0.10, 0.60, 0.36, 0.30], label='A', draw_func=panel_proliferating_ccdf)
pf.add_panel([0.57, 0.60, 0.36, 0.30], label='B', draw_func=panel_arrested_ccdf)
# Row 2 — UMAP (smaller) | GMP-Cor bar plot
pf.add_panel([0.06, 0.08, 0.46, 0.36], label='C', draw_func=panel_A)   # UMAP
pf.add_panel([0.64, 0.12, 0.28, 0.26], label='D', draw_func=panel_C)   # GMP-Cor bar plot
```

- **A** Proliferating CCDF (`REP_LOW`, keeps the y-label), **B** Arrested CCDF (`REP_HIGH`, y-label dropped).
- **C** UMAP (smaller than before), **D** GMP-Cor bar plot.
- `panel_A` (draw function) still renders the UMAP and `panel_C` the bar plot — only their figure
  positions / letters changed.

---

## Panel content

### Panels A & B — representative CCDFs — `_plot_ccdf(ax, sample, title, gmp_cor, signal_color, show_ylabel)`
Same construction/colors as `figure2.py` / `figure3.py`:
1. `arr = np.load(ev_data/<sample>.npy)`; `data1 = arr[0]` (empirical), `data2 = arr[1]` (scrambled), both `> 0`.
2. `x2 = max(data2)` = scrambled maximum = GMP-Cor threshold (`λ_max^scr`).
3. CCDF each: sort ascending, `ccdf = 1 - arange(1,P+1)/P + 1/P`.
4. Three loglog series: `data1 < x2` → **grey** (spurious); `data1 ≥ x2` → **condition color** (signal);
   `data2` → **black** α0.5 (scrambled).
5. Dashed vertical line at `x2`, labeled `λ_max^scr`.
6. **GMP-Cor annotation**: `ax.text(0.04, 0.05, f'GMP-Cor: {gmp_cor:.1f}', weight='bold')` (lower-left);
   legend at **upper-right**.
7. `set_xlim([0.1, 200])`, `ylim(top=1.5)`.

- **A** = `panel_proliferating_ccdf` → `REP_LOW`, title "Proliferating", steelblue signal, keeps the
  y-label, GMP-Cor ≈ 181.7.
- **B** = `panel_arrested_ccdf` → `REP_HIGH`, title "Arrested", salmon signal, y-label dropped,
  GMP-Cor ≈ 125.6.

### Panel C — UMAP (`panel_A` draw function)
- Scatter of precomputed `UMAP_X/Y`, colored by subpopulation; drawn `high → med → low` so the
  sparser groups stay visible on top. `s=2`, `alpha=0.5`, `rasterized=True` (keeps the vector PDF small).
- Legend entries are the plain `COND_DISPLAY[cond]` names (no marker level), `markerscale=4`, lower-left.
- Axes ticks removed. (Cluster `Cycling`/`Non-cycling` annotations were removed — `CLUSTER_LABEL` is
  no longer drawn.)

### Panel D — GMP-Cor bar plot (`panel_C` draw function)
- For each marker level in `COND_ORDER`: mean of the two replicate `GMP_Cor` values, error bar = **SEM**
  (`std(ddof=1)/sqrt(n)`, n = 2). Bars colored by `COND_COLOR`, black edge, `width=0.3`, `ylim(bottom=50)`.
  x-tick labels are the plain `COND_DISPLAY` names (no marker level); no x-axis label.
- Individual replicate points overlaid as black dots (`scatter`, zorder 3).
- Trend: **Proliferating > Intermediate > Arrested** — the proliferating fraction is the most
  globally coordinated.

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Change representative CCDF samples | `REP_HIGH` (Arrested, panel B) / `REP_LOW` (Proliferating, panel A) |
| Recolor subpopulations | `COND_COLOR` |
| Rename proliferation display labels | `COND_DISPLAY` |
| Reorder bars | `COND_ORDER` |
| Re-add UMAP cluster annotations | draw from `CLUSTER_LABEL` inside `panel_A` (currently not drawn) |
| Rearrange panels / resize UMAP | the four `pf.add_panel([...])` rects in the assembly block |
| Move GMP-Cor annotation | the `ax.text(0.04, 0.05, ...)` in `_plot_ccdf` |
| Change error bars (e.g. std instead of SEM) | the `sems.append(...)` line in `panel_C` |
| Adjust UMAP point size / opacity | `s` / `alpha` in `panel_A` scatter |
| Regenerate the underlying data | re-run `scripts/cancer_pc9_correlation_analysis.py` |
| Save as SVG | change `pf.save("figure_s8.pdf", ...)` to `.svg` |

---

## Notes
- CCDF construction and color scheme are kept identical to `figure2.py` / `figure3.py` / `figure_s5.py`.
- UMAP coordinates are **precomputed** in `obs_metadata.csv`; the script does not run UMAP.
- Preview is saved alongside the script at `scripts/supplementary_figures/figure_s8_preview.png`.
