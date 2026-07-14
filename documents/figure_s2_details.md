# Figure S2 — Implementation Details

Supplementary figure explaining the **synthetic scRNA-seq data generator** end to end:
a schematic of the correlation-structure design, representative covariance heatmaps, the
count-generation pipeline, and representative rendered output (cell/gene rank plots)
compared against a real dataset.

## How to run

The script lives in `scripts/supplementary_figures/` (not `scripts/figures/`).

```bash
cd scripts/supplementary_figures
python figure_s2.py        # interactive window; also writes figure_s2.pdf + figure_s2_preview.png
```

An import bootstrap at the top of the script makes it path-independent:

```python
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                # repo root
sys.path.insert(0, _REPO)                                      # -> import src.simulations
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))  # -> import figure_functions
```

Inserting `_REPO` first also shadows any installed `src` package copy that lacks `simulations.py`.
The figure needs the real dataset `data_for_paper/sample_2b_filtered.csv` (for the row-3
comparison rank plots); every other panel is generated at runtime.

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | base font size (titles = 10pt, labels/ticks = 8pt, inset/annotation = 7pt) |
| `figsize` | `(7, 6.5)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(-0.02, 0.02)` | `PanelFigure(label_offset=...)` — panel A/B/C… letters |
| `GENE_C` | `#9ecae1` | ordinary gene node color (schematic) |
| `HUB_C` | `#E07B54` | cluster-hub color |
| `GHUB_C` | `#c0392b` | global-hub color |
| `NB_C` | `steelblue` | count-distribution color (flow-chart insets only) |
| `RANK_C` | `#404040` | dark grey — **all rank plots** (row 3) |
| `_HEATMAP_PARAMS` | `dict(n=2000, shape=1.5, hub_probability=0.2, seed=31)` | params for the row-1 heatmaps |
| `_N_SMALL` | `400` | size of the correlation matrix + simulation behind the examples |
| `_REAL_DATA_PATH` | `data_for_paper/sample_2b_filtered.csv` | real dataset for the comparison rank plots |

---

## Pre-computed data (module level, generated once)

| Object | How | Used by |
|---|---|---|
| `R_high` | `generate_gram_hub_matrix(alpha=0.9, **_HEATMAP_PARAMS)` | Panel B |
| `R_low` | `generate_gram_hub_matrix(alpha=0.5, **_HEATMAP_PARAMS)` | Panel C |
| `R_small` | `generate_gram_hub_matrix(n=400, alpha=0.9, shape=1.5, hub_probability=0.2, seed=31)` | drives the sim |
| `true_counts`, `obs_counts` | `simulate_scRNA_data(n_cells=400, n_genes=400, sigma=R_small, dropout_rate=1, ...)` | `obs_counts` → Panel D (observed inset) |
| `rank_obs_counts` | `simulate_scRNA_data(n_cells=1000, n_genes=2000, sigma=R_high, dropout_rate=1, seed=31)` — defaults `inv_gamma_shape=1.5, inv_gamma_scale=0.01`, matching the `rho_sweep` GMP-Cor calibration | Panels E, F (simulated) — shaped and parameterized to match the data |
| `sim_cell_totals`, `sim_gene_totals` | row/column sums of `rank_obs_counts` | Panels E, F (simulated) |
| `real_cell_totals`, `real_gene_totals` | row/column sums of `sample_2b_filtered.csv` (1041×2071) | Panels E, F (data) |
| `_CELL_YLIM`, `_GENE_YLIM` | `_pair_ylim(...)` — shared y-limits per sim↔data pair | Panels E, F |

`generate_gram_hub_matrix` and `simulate_scRNA_data` are imported from `src/simulations.py`.

---

## Layout map

Positions are `[left, bottom, width, height]` in figure-normalized coords. 3-row narrative;
bold step headers are placed at the left margin (`pf.fig.text`, x=0.02) above each row.

```
 1 · Design correlation structure                         (header y≈0.925)
   A [0.02,0.64,0.30,0.24] network schematic
   B [0.41,0.66,0.23,0.20] heatmap χ=0.9   C [0.72,0.66,0.23,0.20] heatmap χ=0.5
 2 · Generate counts                                      (header y≈0.575)
   D [0.05,0.36,0.90,0.17] count-generation flow (3 inset mini-plots + arrows)
 3 · Representative output vs. data                       (header y≈0.255)
   E [0.09,0.07,0.16,0.15] sim cell rank + [0.28,...] data cell rank
   F [0.57,0.07,0.16,0.15] sim gene rank + [0.76,...] data gene rank
```

Schematic panels (A, D) are added with `hide_axis=True` and draw matplotlib primitives /
`inset_axes`. Panel letters are passed explicitly (`label='A'` …). The two **data** rank
panels are paired next to their simulated counterparts and are added with `label=' '`
(a single space) to suppress the auto-label — otherwise `PanelFigure` auto-assigns stray
letters to any panel called with `label=None`.

---

## Panel A — Cluster + hub factor-model network

`hide_axis=True`. Draws 4 gene clusters (7 dots each, star edges to a cluster hub) and a central
**global hub** (`*` marker) linked to a subset of cluster hubs (`connect = [T, F, T, T]`, i.e. the
`hub_probability` idea). Annotations: "gene", "cluster hub", "global hub". RNG seed `1` (fixed layout).
Title `'Cluster + hub factor model'`.

## Panels B, C — Representative correlation matrices (covariance heatmaps)

Shared helper `_heatmap(ax, R, title)`: `imshow(R[:100,:100], cmap='RdBu_r', vmin=-1, vmax=1)`,
colorbar (`fraction=0.046, pad=0.04`, ticks ±1, 0). B = χ=0.9 (clear block structure),
C = χ=0.5 (weaker). Same color scale for honest comparison. (These are the heatmaps that used to
live in Figure 3.)

## Panel D — Count generation

`hide_axis=True`. Three `inset_axes` (`xs=[0.04,0.39,0.74]`, `w=0.22, h=0.60, y0=0.08`)
joined by arrows with operator labels:
1. correlated latent **MVN** scatter (2 genes, corr 0.8) — `np.random.default_rng(7)`
2. → `Φ, NB⁻¹` → **NB counts** histogram (`nbinom.ppf`, μ=3, r=0.5)
3. → `scale, drop` → **observed counts** heatmap (`obs_counts[:40,:40]`, `Greys`,
   vmax=95th pct, annotated zero-fraction) — the actual generator output.

The intermediate uniform-copula scatter was removed; the "Gaussian copula" label was dropped
from both the panel and the step header. No panel-level title (the step header serves as the title).

## Panels E, F — Representative output vs. data (rank plots)

Shared helper `_rank_plot(ax, totals, title, xlabel, show_ylabel=True, drop_top=False, ylim=None)`:
sorts totals descending (`_rank`), plots on a **log-scale y-axis** (`semilogy`) vs. rank, dark
grey (`RANK_C`). Each panel is a **pair** — simulated (left) next to real data (right):

| Panel | Left (Simulated) | Right (Data) |
|---|---|---|
| E | per-cell total (`sim_cell_totals`) | `real_cell_totals` |
| F | per-gene total (`sim_gene_totals`) | `real_gene_totals` |

Three details make the comparison honest:
- **Shape- and parameter-matched simulation.** The simulated side uses `rank_obs_counts` — a
  dedicated 1000-cell × 2000-gene simulation (`obs_counts`, i.e. post library-scaling and
  dropout) whose count-model parameters mirror the `rho_sweep` GMP-Cor calibration
  (`dropout_rate=1`, `seed=31`, default `inv_gamma_shape=1.5`/`inv_gamma_scale=0.01`). This
  makes both the axis extent and the total-expression marginals track the real data.
- **Matched y-axis.** Each pair shares y-limits (`_CELL_YLIM`, `_GENE_YLIM` from `_pair_ylim`);
  the right (data) panel passes `show_ylabel=False`, which also hides its (now duplicate) y-tick
  labels via `tick_params(labelleft=False)`.
- **Outlier trimming (genes).** The gene-rank panels pass `drop_top=True`, so `_rank` omits the
  single highest-expression gene (which otherwise stretches the log scale); ranks keep their
  original index and therefore start at 2.

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize whole figure | `figsize=(W,H)` in `PanelFigure(...)` |
| Change all font sizes | `fsize` at top |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Move/relabel a step header | edit the `for y, txt in [...]` list near the bottom |
| Change heatmap χ values | edit `R_high` / `R_low` `alpha=` |
| Change heatmap structure | edit `_HEATMAP_PARAMS` |
| Change the representative simulation | edit `_N_SMALL`, `R_small`, or the `simulate_scRNA_data(...)` call |
| Change the real comparison dataset | edit `_REAL_DATA_PATH` |
| Change copula example correlation | edit `cov` in `panel_E` |
| Recolor rank plots | `RANK_C` |
| Recolor schematic nodes | `GENE_C` / `HUB_C` / `GHUB_C` |
| Save as SVG | change `pf.save("figure_s2.pdf", ...)` to `.svg` |

---

## Notes

- This figure absorbs the two covariance heatmaps removed from Figure 3.
- Panel labels are contiguous **A–F**; the old "A→R construction" panel and the standalone
  observed-count matrix panel were removed, and the observed-count render now lives inside the
  Panel D flow.
- Preview is saved alongside the script at `scripts/supplementary_figures/figure_s2_preview.png`
  (written **before** the PDF, so an open/locked `figure_s2.pdf` doesn't block the preview).
