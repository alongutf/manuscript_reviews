"""
Supplementary Figure S5 — Correlation-spectrum CCDFs across all datasets.

Recreates the "# additional data" panel from `scripts/supplementary figures.ipynb`,
but updated to the current figure conventions:

  * Each dataset is shown as a CCDF (complementary CDF, 1 - CDF) of the
    correlation eigenvalue spectrum ONLY -- no eigenvalue-density (PDF) panel --
    in the same loglog style as figure2.py / figure3.py (spurious / signal /
    scrambled, with the scrambled-maximum threshold marked).
  * Each panel is annotated with the GMP-Cor metric read from
    `results/data_metrics/data_metrics.csv`, column `sum_denoised_ev`.
    The signal portion of the CCDF is colored by the dataset `category`
    (regulated vs. dis-arrest), matching figure3.py.

Eigenvalue spectra are loaded from `ev_data/<dataset>.npy`, each an array of
shape (2, P): row 0 = empirical eigenvalues, row 1 = scrambled eigenvalues.

Run from this directory:
    cd scripts/supplementary_figures
    python figure_s5.py
The figure is written next to this script as figure_s5.pdf.
"""
import os
import sys

# --- import bootstrap: this script lives in scripts/supplementary_figures/ ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                       # repo root
sys.path.insert(0, _REPO)                                            # import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))        # figure_functions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Global style
# ------------------------------------------------------------------
fsize = 10
plt.close("all")

EV_DATA_DIR = os.path.join(_REPO, 'ev_data')
METRICS_CSV = os.path.join(_REPO, 'results', 'data_metrics', 'data_metrics.csv')
TITLES_XLSX = os.path.join(_REPO, 'ev_data', 'titles.xlsx')

# Signal color by dataset category (matches figure3.py: regulated vs dis-arrest)
REG_COLOR = 'steelblue'
DIS_COLOR = '#E07B54'
CATEGORY_COLOR = {'r': REG_COLOR, 'd': DIS_COLOR}


def _plot_ccdf(ax, npy_path, title, gmp_cor, signal_color='skyblue',
               show_xlabel=True, show_ylabel=True, show_legend=False):
    """Draw the loglog CCDF of empirical (row 0) vs scrambled (row 1) eigenvalues.

    Color scheme (shared with figure2.py / figure3.py):
      grey  -> spurious correlations (eigenvalues below the scrambled maximum)
      color -> true correlation signal (eigenvalues at/above the scrambled max)
      black -> scrambled data
    """
    arr = np.load(npy_path)
    data1 = arr[0, :]; data1 = data1[data1 > 0]      # empirical eigenvalues
    data2 = arr[1, :]; data2 = data2[data2 > 0]      # scrambled eigenvalues
    x2 = float(np.max(data2))                         # scrambled maximum (GMP-Cor threshold)

    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2

    noise = d1s < x2
    ax.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='spurious', markersize=2)
    ax.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
              color=signal_color, label='signal', markersize=2)
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled', markersize=2)

    ax.set_xlim([0.1, 30])
    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)

    ax.set_title(title, fontsize=fsize - 3, pad=2)
    ax.tick_params(labelsize=fsize - 3, pad=1)
    if show_xlabel:
        ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    if show_ylabel:
        ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)

    # GMP-Cor annotation (the new metric: sum_denoised_ev), on a rounded,
    # semi-transparent box tinted by category (salmon = dis-arrest,
    # steelblue = regulated) via the panel's signal_color.
    ax.text(0.04, 0.05, f'GMP-Cor: {gmp_cor:.2f}', transform=ax.transAxes,
            fontsize=fsize - 3, weight='bold', va='bottom', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=signal_color,
                      alpha=0.5, edgecolor='none'))

    if show_legend:
        ax.legend(fontsize=fsize - 4, loc='upper right', framealpha=0.9)


# ------------------------------------------------------------------
# Collect the datasets to plot.
#
# `ev_data/titles.xlsx` is the driver: it lists which samples to show, the
# display `title` for each, the `category` (signal color), and the order.
# The GMP-Cor value for each sample is looked up from data_metrics.csv
# (column `sum_denoised_ev`) by matching the file-name stem.
# ------------------------------------------------------------------
def _stem(name):
    """Strip a trailing .npy / .csv extension."""
    for ext in ('.npy', '.csv'):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


metrics = pd.read_csv(METRICS_CSV, index_col=0)
gmp_by_stem = {_stem(fn): float(v)
               for fn, v in zip(metrics['file_name'], metrics['sum_denoised_ev'])}

titles = pd.read_excel(TITLES_XLSX)

datasets = []
for _, row in titles.iterrows():
    stem = _stem(str(row['file_name']))
    npy_path = os.path.join(EV_DATA_DIR, stem + '.npy')
    if not os.path.exists(npy_path):
        print(f"WARNING: no eigenvalue data for '{stem}' (expected {npy_path}); skipping")
        continue
    if stem not in gmp_by_stem:
        print(f"WARNING: no GMP-Cor (sum_denoised_ev) for '{stem}' in data_metrics.csv; skipping")
        continue
    datasets.append({
        'stem': stem,
        'title': str(row['title']),
        'npy_path': npy_path,
        'gmp_cor': gmp_by_stem[stem],
        'category': str(row['category']),
    })

n = len(datasets)
print(f"Plotting CCDFs for {n} datasets")

# ------------------------------------------------------------------
# Assemble — portrait grid of small-multiple CCDF panels
# ------------------------------------------------------------------
ncols = 3
nrows = int(np.ceil((n + 1) / ncols))   # +1 cell reserved for a shared legend

pf = PanelFigure(figsize=(7, 1.7 * nrows + 0.3), label_offset=(-0.06, 0.03))
# label=" " (a non-empty, blank string) suppresses the grid's auto "A" label;
# per-panel letters are drawn individually in the loop below.
axes = pf.add_grid_panel([0.07, 0.05, 0.90, 0.91], nrows, ncols,
                         wspace=0.3, hspace=0.45, label=" ")

axes_flat = axes.flatten()
for i, ds in enumerate(datasets):
    r, c = divmod(i, ncols)
    ax = axes_flat[i]
    _plot_ccdf(
        ax, ds['npy_path'], ds['title'], ds['gmp_cor'],
        signal_color=CATEGORY_COLOR.get(ds['category'], 'skyblue'),
        show_xlabel=(r == nrows - 1) or (i + ncols >= n),  # bottom of its column
        show_ylabel=(c == 0),
    )
    # Per-panel letter (A, B, C, ...)
    ax.text(-0.18, 1.12, chr(ord('A') + i), transform=ax.transAxes,
            fontsize=fsize, fontweight='bold', va='top', ha='left')

# Use the trailing empty cell(s) for a single shared legend
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0], [0], marker='.', linestyle='-', color='darkgray',
           label='spurious', markersize=6),
    Line2D([0], [0], marker='.', linestyle='-', color=REG_COLOR,
           label='signal (regulated)', markersize=6),
    Line2D([0], [0], marker='.', linestyle='-', color=DIS_COLOR,
           label='signal (dis-arrest)', markersize=6),
    Line2D([0], [0], marker='.', linestyle='-', color='black', alpha=0.5,
           label='scrambled', markersize=6),
    Line2D([0], [0], linestyle='--', color='k', alpha=0.6,
           label=r'$\lambda_\mathrm{max}^\mathrm{scr}$'),
]
for j in range(n, nrows * ncols):
    axes_flat[j].set_axis_off()
legend_ax = axes_flat[n] if n < nrows * ncols else None
if legend_ax is not None:
    legend_ax.legend(handles=legend_handles, loc='center', fontsize=fsize - 2,
                     frameon=False, title='Correlation spectrum',
                     title_fontsize=fsize - 2)

pf.save("figure_s5.pdf", dpi=300, transparent=True)
pf.fig.savefig("figure_s5_preview.png", dpi=200)
plt.show()
