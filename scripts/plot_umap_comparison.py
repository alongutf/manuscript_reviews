"""
Compare the SHX/casp UMAP across three runs, all drawn in the same format and
colours:

  1. Published         data_for_umap barcodes, as used by figure1.py panels G/H
                       (scanpy/umap_coordinates.csv, i.e. data/scanpy_shx.h5ad)
  2. Original barcodes  data_for_umap barcodes, uniform target_sum, genes
                       removed before normalisation
                       (umap_paper_barcodes.py --barcodes umap)
  3. Paper barcodes    data_for_paper barcodes, same pipeline as row 2
                       (umap_paper_barcodes.py --barcodes paper)

Row 1 vs 2 shows what remains of the published embedding once the pipeline is
re-run on the same cells; rows 2 vs 3 isolate the effect of the cell set.

Similarity of each row to row 1 is reported as a Procrustes disparity (0 =
identical embedding up to translation/rotation/scale) and as the adjusted Rand
index between cluster labels.  Note that `scanpy/umap_coordinates_shx_scaled.csv`
(with the per-condition 3802/3507/46 scalings) is a DIFFERENT run and is not the
published one.

Every panel is drawn by the same function from saved coordinate CSVs, so any
visual difference is a difference in the data, not in the plotting.

Colours: scanpy's default categorical palette (sc.pl.palettes.default_20),
assigned in the same category order used by the pipeline
(batches: exp, dis1, dis2, reg1, reg2; clusters: 0, 1, 2, ...).

Outputs (timestamped):
  scanpy/umap_shx_published_<ts>.svg / .png         (row 1 alone)
  scanpy/umap_shx_three_way_comparison_<ts>.svg / .png
"""

import glob
import os
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import procrustes
from sklearn.metrics import adjusted_rand_score

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(REPO, 'scanpy')

PUBLISHED_CSV = os.path.join(OUT_DIR, 'umap_coordinates.csv')
BATCH_ORDER = ['exp', 'dis1', 'dis2', 'reg1', 'reg2']
PALETTE = list(sc.pl.palettes.default_20)

TS = datetime.now().strftime('%Y%m%d_%H%M%S')


def latest_csv(src):
    """Most recent coordinates CSV for a barcode source ('paper' or 'umap')."""
    hits = sorted(glob.glob(os.path.join(
        OUT_DIR, f'umap_coordinates_shx_{src}_barcodes_*.csv')))
    if not hits:
        raise FileNotFoundError(
            f'no {src}-barcode coordinates found; run '
            f'`python umap_paper_barcodes.py --barcodes {src}` first')
    return hits[-1]


POINT_SIZE = 4
POINT_ALPHA = 0.7


def panel(ax, df, color_by, title, categories=None):
    """One UMAP panel: small semi-transparent points, drawn in shuffled order.

    Plotting each category as its own scatter puts the last category entirely on
    top of the earlier ones, which hides whatever it overlaps.  Instead every
    point goes into a single scatter in a fixed random order, so overlapping
    groups are interleaved and the transparency shows what is underneath.
    """
    cats = categories if categories is not None else sorted(df[color_by].unique())
    colors = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(cats)}

    order = np.random.default_rng(0).permutation(len(df))
    d = df.iloc[order]
    ax.scatter(d.UMAP_1, d.UMAP_2, s=POINT_SIZE, linewidths=0, alpha=POINT_ALPHA,
               c=[colors[v] for v in d[color_by]])
    # legend proxies, since the single scatter carries no per-category label
    for c in cats:
        ax.scatter([], [], s=40, linewidths=0, color=colors[c], label=str(c))
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),
              frameon=False, handletextpad=0.1)


def two_panel(df, tag, path_stem):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    panel(axes[0], df, 'cluster', 'Leiden clusters')
    panel(axes[1], df, 'batch', 'Sample', categories=BATCH_ORDER)
    fig.suptitle(tag)
    fig.tight_layout()
    fig.savefig(path_stem + '.svg')
    fig.savefig(path_stem + '.png', dpi=200)
    plt.close(fig)
    print('wrote', path_stem + '.svg')


def similarity(ref, df):
    """Procrustes disparity and cluster ARI against the published embedding."""
    shared = ref.index.intersection(df.index)
    _, _, disp = procrustes(ref.loc[shared, ['UMAP_1', 'UMAP_2']].values,
                            df.loc[shared, ['UMAP_1', 'UMAP_2']].values)
    ari = adjusted_rand_score(ref.loc[shared, 'cluster'], df.loc[shared, 'cluster'])
    return disp, ari, len(shared)


def main():
    runs = [
        ('Published\n(data_for_umap barcodes,\nfigure 1G/H)', PUBLISHED_CSV),
        ('Original barcodes\n(data_for_umap,\nuniform, genes removed first)',
         latest_csv('umap')),
        ('Paper barcodes\n(data_for_paper,\nuniform, genes removed first)',
         latest_csv('paper')),
    ]
    frames = []
    for label, path in runs:
        df = pd.read_csv(path, index_col=0)
        frames.append((label, df))
        print(f'{label.splitlines()[0]:<20} {len(df)} cells, '
              f'{df.cluster.nunique()} clusters ({os.path.basename(path)})')

    ref = frames[0][1]
    print('\nsimilarity to the published embedding:')
    for label, df in frames[1:]:
        disp, ari, n = similarity(ref, df)
        print(f'  {label.splitlines()[0]:<20} procrustes={disp:.4f}  '
              f'ARI={ari:.3f}  (n={n} shared cells)')

    two_panel(ref, 'Published UMAP (scanpy/umap_coordinates.csv)',
              os.path.join(OUT_DIR, f'umap_shx_published_{TS}'))

    # 3x2 comparison: rows = run, columns = colouring
    fig, axes = plt.subplots(3, 2, figsize=(13, 16.5))
    for row, (label, df) in enumerate(frames):
        panel(axes[row, 0], df, 'cluster', 'Leiden clusters')
        panel(axes[row, 1], df, 'batch', 'Sample', categories=BATCH_ORDER)
        axes[row, 0].set_ylabel(f'{label}\nn={len(df)}', fontsize=11)
        axes[row, 0].yaxis.set_visible(True)
        axes[row, 0].set_yticks([])
    fig.tight_layout()
    stem = os.path.join(OUT_DIR, f'umap_shx_three_way_comparison_{TS}')
    fig.savefig(stem + '.svg')
    fig.savefig(stem + '.png', dpi=200)
    plt.close(fig)
    print('wrote', stem + '.svg')


if __name__ == '__main__':
    main()
