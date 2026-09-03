"""
Split figure2 panels C/D into three standalone presentation figures.

Reuses figure2.py's exact data generation and the `_draw_ccdf` helper by
exec'ing the script up to its `pf = PanelFigure` line (no full-figure build).
figure2.py itself is left untouched.

Outputs (transparent SVG, GRC folder):
  1. figure2C_heatmap_tall.svg  - the 2C top matrix only, tall & narrow.
  2. figure2D_main_hist.svg     - the 2D top histogram only, no lambda_max^MP
                                  line; grey 'spurious' up to lambda_max^scr,
                                  blue 'signal' above it.
  3. figure2D_ccdf_inset.svg    - the CCDF inset of 2D top as its own figure.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
os.chdir(HERE)  # figure2.py's data paths are relative to scripts/figures/

import numpy as np
import matplotlib
matplotlib.use('Agg')  # headless backend; this script only ever writes files
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# external, machine-specific export target (a synced cloud drive), not part of
# this repo -- update or redirect this path to run the export elsewhere
OUT_DIR = r'G:\Other computers\My MacBook Air\Alon\PhD\documents\GRC conference\figures'
FS = 12   # base font size for these presentation figures


def load_figure2_namespace():
    """Execute figure2.py's module-level code up to (but not including) the
    line that builds its `pf = PanelFigure(...)` full-figure layout, and return
    the resulting namespace.

    This reuses figure2.py's data loading and its `_draw_ccdf` helper without
    running or reproducing its panel-layout logic, and without ever importing
    figure2.py as a module (which would execute the whole file, panels
    included). It relies on the literal text 'pf = PanelFigure' appearing at
    the start of a line in figure2.py -- if that line is reworded or the
    variable renamed, `next()` raises StopIteration and this script breaks.
    """
    path = os.path.join(HERE, 'figure2.py')
    with open(path, 'r', encoding='utf-8') as fh:
        lines = fh.readlines()
    cut = next(i for i, ln in enumerate(lines) if ln.lstrip().startswith('pf = PanelFigure'))
    ns = {'__name__': '_grc_export', '__file__': path}
    exec(compile(''.join(lines[:cut]), path, 'exec'), ns)
    ns['fsize'] = FS
    return ns


ns = load_figure2_namespace()
ev_data_dir = ns['ev_data_dir']  # directory of pre-computed eigenvalue-spectrum .npy files, from figure2.py
_draw_ccdf = ns['_draw_ccdf']    # figure2.py's shared CCDF-plotting helper


def save(fig, name, tight_bbox=True):
    # writes the real (transparent SVG) output to OUT_DIR, plus a PNG preview
    # alongside this script for a quick local look without opening OUT_DIR
    bbox = 'tight' if tight_bbox else None
    svg = os.path.join(OUT_DIR, name + '.svg')
    fig.savefig(svg, transparent=True, bbox_inches=bbox)
    fig.savefig(os.path.join(HERE, '_grc_preview_' + name + '.png'),
                dpi=200, bbox_inches=bbox)
    plt.close(fig)
    print('wrote', svg)


# ------------------------------------------------------------------
# 1. Heatmap (2C top) — tall & narrow
# ------------------------------------------------------------------
def make_heatmap():
    # this panel is a schematic, not real data: a synthetic sparse cells x genes
    # matrix (per-gene detection rate 15-38%, exponential nonzero values) just to
    # illustrate what a raw scRNA-seq count matrix looks like
    np.random.seed(42)
    n_cells, n_genes = 30, 12
    matrix = np.zeros((n_cells, n_genes))
    for j in range(n_genes):
        n_nz = max(1, int(n_cells * np.random.uniform(0.15, 0.38)))
        rows = np.random.choice(n_cells, size=n_nz, replace=False)
        matrix[rows, j] = np.random.exponential(2.5, size=n_nz)
    vmax = matrix.max()

    fig, ax = plt.subplots(figsize=(1.2, 3.0))
    ax.imshow(matrix, aspect='auto', cmap='Greys', vmin=0, vmax=vmax,
              interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('genes', fontsize=FS, labelpad=3)
    ax.set_ylabel('cells', fontsize=FS, labelpad=3)
    ax.set_title('scRNA-seq data', fontsize=FS, pad=4)
    save(fig, 'figure2C_heatmap_tall')


# ------------------------------------------------------------------
# Shared 2D-top data
# ------------------------------------------------------------------
pcs_data = np.load(os.path.join(ev_data_dir, 'simulated_pcs.npy'))
pcs, pcs1 = pcs_data[0], pcs_data[1]  # row 0 = real-data eigenvalues, row 1 = scrambled-null eigenvalues
data1 = pcs[pcs > 0]
data2 = pcs1[pcs1 > 0]
bin_width = 0.15
# scrambled maximum = the GMP-Cor threshold; offset by one bin width so the
# dividing line sits at a histogram bin edge rather than mid-bin
x2 = float(np.max(pcs1)) + bin_width

SPURIOUS_COLOR = 'darkgray'
SIGNAL_COLOR = 'skyblue'


# ------------------------------------------------------------------
# 2. Main histogram (2D top) — grey spurious up to x2, blue signal above
# ------------------------------------------------------------------
def make_hist():
    # bin edges span both the real and scrambled data so the two histograms
    # would align if plotted together, even though only data1 (real) is drawn
    all_data = np.concatenate([data1, data2])
    bin_edges = np.arange(min(all_data), max(all_data) + bin_width, bin_width)

    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    _, _, patches = ax.hist(
        data1, bins=bin_edges, width=bin_width * 0.8, align='right',
        edgecolor='black', color='#d9d9d9', alpha=0.7, density=True)
    # recolor bars by whether they lie below (spurious/noise) or above
    # (signal) the scrambled-null threshold x2
    for patch in patches:
        patch.set_facecolor(SPURIOUS_COLOR if patch.get_x() < x2 else SIGNAL_COLOR)

    # keep only the scrambled-max threshold line (MP line removed)
    ax.axvline(x2, color='dimgray', linestyle=':', alpha=0.8)

    ax.set_xlabel(r'$\lambda$', fontsize=FS, labelpad=0)
    ax.set_ylabel(r'$\rho(\lambda)$', fontsize=FS)
    ax.set_title('Correlation eigenvalue density', fontsize=FS, pad=4)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xlim([2, 12])
    ax.set_ylim(0, 0.2)
    ax.tick_params(axis='both', which='major', labelsize=FS - 2)
    ax.grid(False)

    legend_handles = [
        Patch(facecolor=SPURIOUS_COLOR, edgecolor='black', label='spurious'),
        Patch(facecolor=SIGNAL_COLOR, edgecolor='black', label='signal'),
    ]
    ax.legend(handles=legend_handles, fontsize=FS - 2, loc='upper right',
              framealpha=1)
    save(fig, 'figure2D_main_hist')


# ------------------------------------------------------------------
# 3. CCDF inset of 2D top, as its own figure
# ------------------------------------------------------------------
def make_ccdf():
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    _draw_ccdf(ax, data1, data2, show_legend=True, markersize=3)
    ax.set_xlabel(r'$\lambda$', fontsize=FS, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=FS, labelpad=0)
    ax.tick_params(labelsize=FS - 2)
    ax.set_xlim([0.2, 14])
    ax.set_title('Transform to CCDF', fontsize=FS, pad=4)
    save(fig, 'figure2D_ccdf_inset')


make_heatmap()
make_hist()
make_ccdf()
print('done')
