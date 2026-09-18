"""Clean two-panel version of the dis-arrest cluster figure.

Reads the per-cell table and the JSON log written by
scripts/dis_arrest_cluster_test.py (results/dis_arrest_clusters/) and redraws only
the two UMAP panels of that run's six-panel summary:

  A  UMAP colored by leiden cluster  (cluster 1 = leiden 0, cluster 2 = leiden 1)
  B  UMAP colored by sample of origin (sample 1 = dis1, sample 2 = dis2)

The two clusters of panel A carry a GMP-Cor text box in the style of supplementary
figure S12 (scripts/supplementary_figures/figure_s12.py): the matched-n GMP-Cor of
that cluster, in a rounded white box sitting on the cluster's own centroid. The value
is the MEDIAN over the run's independent cell draws, not the mean -- the draw
distribution is right-skewed, so the mean is pulled up by a few outlier-cell draws.
Panel B is unannotated -- a per-sample GMP-Cor is not what that panel is about. The
leiden resolution is not shown.

Outputs (results/dis_arrest_clusters/figures/):
  <stem>_panelsAB_<timestamp>.svg / .png
"""
import os
import sys
import json
import argparse
import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'figures'))
from figure_functions import PanelFigure

OUTDIR = os.path.join(ROOT, 'results', 'dis_arrest_clusters')
DEFAULT_RUN = 'dis_arrest_cluster_test_20260830_150710'

fsize = 10

# palette: the two clusters keep the colors of the original run, the two samples
# get the matplotlib default pair the original panel B used
C_C1 = '#2980b9'        # cluster 1 = leiden 0 (the larger, low-depth cluster)
C_C2 = '#c0392b'        # cluster 2 = leiden 1
C_S1 = '#1f77b4'        # sample 1 = dis1
C_S2 = '#ff7f0e'        # sample 2 = dis2
BOX = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.72, linewidth=0.8)

# the run labels its groups by what they turned out to be; the figure labels them
# neutrally, so keep the translation in one place
CLUSTER_LABEL = {0: 'cluster 1', 1: 'cluster 2'}
CLUSTER_GROUP = {0: 'low-depth cluster', 1: 'high-depth cluster'}
SAMPLE_LABEL = {'dis1': 'sample 1', 'dis2': 'sample 2'}


def _load(path):
    with open(path, encoding='utf8') as fh:
        return json.load(fh)


def _gmp(log, group):
    """Matched-n GMP-Cor of one group: the median over the run's cell draws.

    Falls back to the draws themselves for runs written before the runner logged a
    median, so older runs still redraw.
    """
    d = log['q2_gmp_cor']['per_group'][group]
    if 'gmp_cor_median' in d:
        return float(d['gmp_cor_median'])
    return float(np.median(d['gmp_cor_draws']))


def _bare(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])


def _make_room(ax, xy, xpad=0.06, ypad=0.10):
    """A little headroom around the cloud; the boxes sit on the centroids."""
    x, y = np.asarray(xy[:, 0]), np.asarray(xy[:, 1])
    x0, x1, y0, y1 = x.min(), x.max(), y.min(), y.max()
    xr, yr = x1 - x0, y1 - y0
    ax.set_xlim(x0 - xpad * xr, x1 + xpad * xr)
    ax.set_ylim(y0 - 0.08 * yr, y1 + ypad * yr)
    return x0, x1, y0, y1, xr, yr


def _annotate(ax, box_pt, label, value, color):
    """One S12-style GMP-Cor box, placed directly on the cluster centroid."""
    ax.text(*box_pt, f'{label}\nGMP-Cor = {value:.1f}', fontsize=fsize - 3,
            ha='center', va='center', zorder=3, bbox=dict(**BOX, edgecolor=color))


def _umap_panel(ax, cells, log, key, groups, title, annotate=True):
    """Scatter the UMAP split by `key`, optionally boxing each group's GMP-Cor."""
    for val, label, group, color in groups:
        sub = cells[cells[key] == val]
        ax.scatter(sub.UMAP_1, sub.UMAP_2, color=color, alpha=.8, s=1.2,
                   linewidths=0, zorder=1, label=label)

    _bare(ax)
    xy = cells[['UMAP_1', 'UMAP_2']].to_numpy(dtype=float)
    _make_room(ax, xy)

    if annotate:
        for val, label, group, color in groups:
            cen = cells.loc[cells[key] == val, ['UMAP_1', 'UMAP_2']].mean().to_numpy()
            _annotate(ax, tuple(cen), label, _gmp(log, group), color)

    ax.set_title(title, fontsize=fsize - 1, pad=2)
    leg = ax.legend(fontsize=fsize - 3, frameon=False, loc='lower center',
                    ncol=2, handletextpad=0.2, borderpad=0.1, columnspacing=1.0,
                    markerscale=6)
    for h in leg.legend_handles:
        h.set_alpha(1)


def panel_A(ax, cells, log):
    _umap_panel(ax, cells, log, 'leiden',
                [(0, CLUSTER_LABEL[0], CLUSTER_GROUP[0], C_C1),
                 (1, CLUSTER_LABEL[1], CLUSTER_GROUP[1], C_C2)],
                'Leiden clusters')


def panel_B(ax, cells, log):
    _umap_panel(ax, cells, log, 'batch',
                [('dis1', SAMPLE_LABEL['dis1'], 'dis1', C_S1),
                 ('dis2', SAMPLE_LABEL['dis2'], 'dis2', C_S2)],
                'Sample of origin', annotate=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', default=DEFAULT_RUN,
                    help='run stem in results/dis_arrest_clusters/logs (no extension)')
    args = ap.parse_args()

    log = _load(os.path.join(OUTDIR, 'logs', args.run + '.json'))
    cells = pd.read_csv(os.path.join(OUTDIR, 'logs', args.run + '_cells.csv'),
                        index_col=0)

    pf = PanelFigure(figsize=(9, 3.9), label_offset=(-0.035, 0.03))
    pf.add_panel([0.065, 0.09, 0.39, 0.82], label='A',
                 draw_func=lambda ax: panel_A(ax, cells, log))
    pf.add_panel([0.565, 0.09, 0.39, 0.82], label='B',
                 draw_func=lambda ax: panel_B(ax, cells, log))

    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    stem = os.path.join(OUTDIR, 'figures',
                        'dis_arrest_cluster_panelsAB_' + stamp)
    pf.save(stem + '.svg')
    pf.save(stem + '.png', dpi=300)
    print('wrote:\n  %s.svg\n  %s.png' % (stem, stem))
    print('source run: %s  (matched n = %d cells per group)'
          % (args.run, log['q2_gmp_cor']['matched_n']))
    for label, group in ((CLUSTER_LABEL[0], CLUSTER_GROUP[0]),
                         (CLUSTER_LABEL[1], CLUSTER_GROUP[1])):
        print('  %-10s GMP-Cor = %6.2f' % (label, _gmp(log, group)))


if __name__ == '__main__':
    main()
