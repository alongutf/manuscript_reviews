# Figure S2 — Implementation Details

Supplementary figure explaining the **synthetic scRNA-seq data generator** end to end:
schematics of the correlation-structure design + the Gaussian-copula count pipeline, plus
representative rendered examples (covariance heatmaps, count marginal, heavy tail, sparse matrix).

## How to run

The script lives in `scripts/supplementary_figures/` (not `scripts/figures/`).

```bash
cd scripts/supplementary_figures
python figure_s2.py        # interactive window; also writes figure_s2.pdf next to the script
```

An import bootstrap at the top of the script makes it path-independent:

```python
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                # repo root
sys.path.insert(0, _REPO)                                      # -> import src.simulations
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))  # -> import figure_functions
```

Inserting `_REPO` first also shadows any installed `src` package copy that lacks `simulations.py`.
The figure needs **no** files from `ev_data/` or `results/` — every panel is generated at runtime.

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | base font size (titles = 10pt, labels/ticks = 8pt, inset/annotation = 7pt) |
| `figsize` | `(7.5, 9)` inches | `PanelFigure(figsize=...)` in assembly block |
| `label_offset` | `(-0.025, 0.008)` | `PanelFigure(label_offset=...)` — panel A/B/C… letters |
| `GENE_C` | `#9ecae1` | ordinary gene node color (schematic) |
| `HUB_C` | `#E07B54` | cluster-hub color |
| `GHUB_C` | `#c0392b` | global-hub color |
| `NB_C` | `steelblue` | count-distribution color |
| `_HEATMAP_PARAMS` | `dict(n=2000, shape=1.5, hub_probability=0.2, seed=31)` | params for the row-2 heatmaps |
| `_N_SMALL` | `400` | size of the correlation matrix + simulation behind row-4 examples |

---

## Pre-computed data (module level, generated once)

| Object | How | Used by |
|---|---|---|
| `R_high` | `generate_gram_hub_matrix(alpha=0.9, **_HEATMAP_PARAMS)` | Panel C |
| `R_low` | `generate_gram_hub_matrix(alpha=0.5, **_HEATMAP_PARAMS)` | Panel D |
| `R_small` | `generate_gram_hub_matrix(n=400, alpha=0.9, shape=1.5, hub_probability=0.2, seed=31)` | drives the sim |
| `true_counts`, `obs_counts` | `simulate_scRNA_data(n_cells=400, n_genes=400, sigma=R_small, dropout_rate=0.5, seed=0)` | Panels E (iv), F, G, H |

`generate_gram_hub_matrix` and `simulate_scRNA_data` are imported from `src/simulations.py`.

---

## Layout map

Positions are `[left, bottom, width, height]` in figure-normalized coords. 4-row narrative;
bold step headers are placed at the left margin (`pf.fig.text`, x=0.02) above each row.

```
 1 · Design correlation structure                         (header y≈0.965)
   A [0.05,0.78,0.42,0.14] network schematic   B [0.55,0.78,0.42,0.14] A→R construction
 2 · Representative correlation matrices                  (header y≈0.715)
   C [0.10,0.54,0.34,0.14] heatmap ρ=0.9       D [0.58,0.54,0.34,0.14] heatmap ρ=0.5
 3 · Generate counts (Gaussian copula)                    (header y≈0.465)
   E [0.04,0.30,0.93,0.13] copula flow (4 inset mini-plots + arrows)
 4 · Representative simulated output                      (header y≈0.225)
   F [0.09,0.05,0.21,0.15]  G [0.41,0.05,0.21,0.15]  H [0.71,0.05,0.25,0.15]
```

Schematic panels (A, B, E) are added with `hide_axis=True` and draw matplotlib primitives /
`inset_axes`. Panel letters are passed explicitly (`label='A'` …).

---

## Panel A — Cluster + hub factor-model network

`hide_axis=True`. Draws 4 gene clusters (7 dots each, star edges to a cluster hub) and a central
**global hub** (`*` marker) linked to a subset of cluster hubs (`connect = [T, F, T, T]`, i.e. the
`hub_probability` idea). Annotations: "gene", "cluster hub", "global hub". RNG seed `1` (fixed layout).
Title `'Cluster + hub factor model'`.

## Panel B — Build the correlation matrix

`hide_axis=True`. Two `inset_axes`: left = illustrative loading matrix `A` (`magma`), right =
correlation `R` (`RdBu_r`, vmin/vmax ±1, colorbar). Arrow between them labeled `AAᵀ / normalize`.
`A` and `R` come from `_illustrative_loading_matrix(n=40, alpha=0.8, seed=3)`, which rebuilds the
same `A = √(1-α)I | cluster cols | global-hub col` structure at small scale for legibility.

## Panels C, D — Representative correlation matrices (the covariance heatmaps)

Shared helper `_heatmap(ax, R, title)`: `imshow(R[:100,:100], cmap='RdBu_r', vmin=-1, vmax=1)`,
colorbar (`fraction=0.046, pad=0.04`, ticks ±1, 0). C = ρ=0.9 (clear block structure),
D = ρ=0.5 (weaker). Same color scale for honest comparison. (These are the heatmaps that used to
live in Figure 3.)

## Panel E — Gaussian-copula count generation

`hide_axis=True`. Four `inset_axes` (`xs=[0.025,0.285,0.545,0.805]`, `w=0.165, h=0.58, y0=0.10`)
joined by arrows with operator labels:
1. correlated latent **MVN** scatter (2 genes, ρ=0.8) — `np.random.default_rng(7)`
2. → `Φ` → **uniform** copula scatter (`norm.cdf`)
3. → `NB⁻¹` → **NB counts** histogram (`nbinom.ppf`, μ=3, r=0.5)
4. → `scale, drop` → **sparse** observed matrix (`obs_counts[:40,:40] > 0`)

No panel-level title (the step header serves as the title).

## Panels F, G, H — Representative simulated output

| Panel | Content | Notes |
|---|---|---|
| F | Per-gene **NB count** histogram | gene at the 80th percentile of `true_counts` totals; 20 bins, `NB_C` |
| G | Gene-total **CCDF** (log-log) | heavy tail from inverse-Gamma gene means; helper `_ccdf` (`1 - rank/p + 1/p`) |
| H | Observed **sparse** count matrix | `obs_counts[:120,:120]`, `Greys`, vmax=95th pct; annotated zero-fraction (dropout) |

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Resize whole figure | `figsize=(W,H)` in `PanelFigure(...)` |
| Change all font sizes | `fsize` at top |
| Move a panel | edit `panel_pos[i]` — `[left, bottom, width, height]` |
| Move/relabel a step header | edit the `for y, txt in [...]` list near the bottom |
| Change heatmap ρ values | edit `R_high` / `R_low` `alpha=` |
| Change heatmap structure | edit `_HEATMAP_PARAMS` |
| Change the representative simulation | edit `_N_SMALL`, `R_small`, or the `simulate_scRNA_data(...)` call |
| Change copula example correlation | edit `cov` in `panel_E` |
| Recolor schematic nodes | `GENE_C` / `HUB_C` / `GHUB_C` |
| Save as SVG | change `pf.save("figure_s2.pdf", ...)` to `.svg` |

---

## Notes

- This figure absorbs the two covariance heatmaps removed from Figure 3; the standalone heatmap
  stub `scripts/figures/figure s5.py` was deleted (its code is reproduced here).
- Preview is saved alongside the script at `scripts/supplementary_figures/figure_s2_preview.png`.
