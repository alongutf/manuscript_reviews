# Figure S5 — Implementation Details

Supplementary figure showing the **correlation-spectrum CCDF** for every dataset in the paper.
It is the modernized replacement for the `# additional data` cell in
`scripts/supplementary figures.ipynb`, which plotted an eigenvalue-density (PDF) panel with a
small CCDF inset. Figure S5 keeps **only the CCDF** (complementary CDF, `1 - CDF`) of the
eigenvalue spectrum, drawn in the same loglog style as `figure2.py` / `figure3.py`, and reports
the **current GMP-Cor metric** for each dataset.

## How to run

The script lives in `scripts/supplementary_figures/` (not `scripts/figures/`).

```bash
cd scripts/supplementary_figures
python figure_s5.py        # interactive window; also writes figure_s5.pdf + figure_s5_preview.png
```

An import bootstrap at the top makes the script path-independent:

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
| Sample list / titles | `ev_data/titles.xlsx` | **the driver.** Columns `file_name`, `title`, `category`. Defines *which* samples appear, their display **title**, the `category` (`r` = regulated, `d` = dis-arrest) used to color the signal, and the **plotting order** (row order in the file). |
| Eigenvalue spectra | `ev_data/<dataset>.npy` | one array per sample, shape `(2, P)`: row 0 = empirical eigenvalues, row 1 = scrambled eigenvalues. Drives every CCDF curve. |
| Metrics table | `results/data_metrics/data_metrics.csv` | source of the **GMP-Cor** value per sample (column `sum_denoised_ev`), looked up by filename stem. |

The figure is **driven by `titles.xlsx`**: it iterates over that file's rows in order. For each row it
strips the extension from `file_name` (handles both `.npy` and `.csv`) to get the stem, loads
`ev_data/<stem>.npy`, and looks up `sum_denoised_ev` for that stem from `data_metrics.csv`. A row is
skipped with a printed `WARNING` if either the `.npy` spectrum or the GMP-Cor value is missing. The
current `titles.xlsx` lists 17 samples, all of which resolve.

---

## Global settings

| Variable | Value | Where to change |
|---|---|---|
| `fsize` | `10` | base font size (panel letters = 10pt, axis labels = 8pt, titles/ticks/annotations = 7pt, legend = 8pt) |
| `figsize` | `(7, 1.55*nrows + 0.4)` inches | auto-scales to the number of rows in the assembly block |
| `label_offset` | `(-0.06, 0.03)` | `PanelFigure(label_offset=...)` (grid uses per-axes letters instead, see below) |
| `REG_COLOR` | `steelblue` | signal color for `category == 'r'` (regulated) |
| `DIS_COLOR` | `#E07B54` | signal color for `category == 'd'` (dis-arrest) |
| `CATEGORY_COLOR` | `{'r': REG_COLOR, 'd': DIS_COLOR}` | maps category → signal color |

---

## Layout

A single `add_grid_panel` call builds a `nrows × ncols` grid of small-multiple CCDF panels:

```python
ncols = 4
nrows = ceil((n + 1) / ncols)          # +1 cell reserved for the shared legend
pf.add_grid_panel([0.07, 0.05, 0.90, 0.91], nrows, ncols,
                  wspace=0.3, hspace=0.45, label=" ")
```

With the current 17 datasets this is a **5×4** grid (20 cells): 17 CCDFs + 1 legend cell + 2 blank.
Panels are filled row-major in `titles.xlsx` order. Each gets a bold letter (`A`, `B`, …) drawn
manually at `ax.text(-0.18, 1.12, ...)` (the grid helper only places one label, so per-panel letters
are added in the loop). The `add_grid_panel` call passes `label=" "` (a non-empty blank string) so
the helper does **not** auto-draw an "A" for the whole grid — which would otherwise collide with the
per-panel "A". Trailing empty cells are turned off with `set_axis_off()`; the first empty cell hosts
the shared legend.

To reduce clutter, axis labels are shown selectively: `CCDF` (y-label) only on the **left column**,
`λ` (x-label) only on the **bottom panel of each column**.

---

## Per-panel content — `_plot_ccdf(ax, npy_path, title, gmp_cor, signal_color, ...)`

Same construction as the CCDF helpers in `figure2.py` / `figure3.py`:

1. Load `arr = np.load(npy_path)`; `data1 = arr[0]` (empirical), `data2 = arr[1]` (scrambled), both filtered to `> 0`.
2. `x2 = max(data2)` — the scrambled maximum, i.e. the GMP-Cor threshold (`λ_max^scr`).
3. CCDF of each: sort ascending, `ccdf = 1 - arange(1, P+1)/P + 1/P`.
4. Three loglog series:
   - `data1 < x2` → **grey** (`darkgray`, α 0.7) = spurious correlations
   - `data1 ≥ x2` → **signal color by category** = true correlation signal
   - `data2` → **black** (α 0.5) = scrambled
5. Dashed vertical line at `x2` (`λ_max^scr`).
6. `set_xlim([0.1, 30])`, title = the display `title` from `titles.xlsx`.
7. **GMP-Cor annotation**: `ax.text(0.04, 0.05, f'GMP-Cor: {gmp_cor:.2f}', weight='bold')` — the
   value is `sum_denoised_ev` from `data_metrics.csv`.

## Shared legend

Built from `matplotlib.lines.Line2D` proxies placed in the first empty grid cell: spurious,
signal (regulated), signal (dis-arrest), scrambled, and the `λ_max^scr` dashed line. Title
`'Correlation spectrum'`, `frameon=False`.

---

## Common tweaks quick-reference

| Goal | What to change |
|---|---|
| Add/remove a sample, change order, or retitle | edit `ev_data/titles.xlsx` (`file_name`, `title`, `category` columns) — and ensure a matching `ev_data/<stem>.npy` and a `sum_denoised_ev` row exist |
| Change the reported metric | edit the `gmp_by_stem` lookup (currently `sum_denoised_ev` from `data_metrics.csv`) |
| Change grid width | `ncols` in the assembly block (`nrows` auto-recomputes) |
| Tune panel spacing | `wspace` / `hspace` in `add_grid_panel` |
| Recolor regulated / dis-arrest signal | `REG_COLOR` / `DIS_COLOR` |
| Change all font sizes | `fsize` at top |
| Move the GMP-Cor annotation | the `ax.text(0.04, 0.05, ...)` call in `_plot_ccdf` |
| Save as SVG | change `pf.save("figure_s5.pdf", ...)` to `.svg` |

---

## Notes

- CCDF construction and color scheme are kept identical to `figure2.py` / `figure3.py` so the
  supplementary small-multiples match the main-text panels.
- Signal color encodes `category` (regulated vs dis-arrest), mirroring `figure3.py`.
- Preview is saved alongside the script at `scripts/supplementary_figures/figure_s5_preview.png`.
