"""
UMAP of the concatenated SHX/casp datasets.

Pipeline follows `scripts/scanpy_analysis.ipynb` (first cell, the SHX/casp UMAP)
with three deliberate differences.

(1) Order of operations: the genes in GENES_TO_REMOVE are dropped BEFORE
normalisation.  The notebook normalises each sample first and drops them after,
so `normalize_total` divides by a total that is ~92% 16s rRNA -- the surviving
genes then carry a scale factor driven almost entirely by rRNA content.
Matching the published UMAP (scanpy/umap_coordinates.csv) requires this order:
removal-then-normalisation reaches a Procrustes disparity of 0.054 against it,
versus 0.30-0.45 for normalisation-then-removal.

(2) Normalisation (`--target-sum`, default 1e4): every sample is normalised to
the SAME total.  The notebook used hard-coded per-condition scalings
(exp = 3802, SHX = 3507, casp = 46), a ~76x scale gap imposed before log1p.
Those belong to a different run (scanpy_shx_scaled_by_total_rna.h5ad /
umap_coordinates_shx_scaled.csv), not to the published figure.

(3) Cell set (`--barcodes`):
  paper  the unfiltered matrix subset to the barcodes of
         data_for_paper/<sample>_filtered.csv
  umap   the barcodes used by the original notebook, i.e.
         data_for_umap/<sample>_filtered.csv

`data_for_umap/*_filtered.csv` is *cell*-filtered only: its gene columns are
identical to the unfiltered matrix.  Both modes therefore see the same gene
space, and differ only in which cells are included.

Outputs (timestamped, <src> = paper_barcodes | umap_barcodes):
  scanpy/adata_shx_<src>_<ts>.h5ad
  scanpy/umap_coordinates_shx_<src>_<ts>.csv
  scanpy/umap_shx_<src>_<ts>.svg / .png
  scanpy/marker_genes_shx_<src>_<ts>.xlsx
  scripts/logs/umap_<src>_<ts>.json
"""

import argparse
import json
import os
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# --------------------------------------------------------------------------
# Parameters (kept identical to scanpy_analysis.ipynb)
# --------------------------------------------------------------------------
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# NOTE: absolute path to a sibling repo on this machine that holds the raw
# per-sample unfiltered probe-count CSVs; only used when --barcodes=paper.
# Running this script on another machine requires this path to be updated.
UNFILTERED_DIR = r'C:\Users\owner\Documents\Projects\rnaseq_correlations\data'
PAPER_DIR = os.path.join(REPO, 'data_for_paper')
UMAP_DIR = os.path.join(REPO, 'data_for_umap')
OUT_DIR = os.path.join(REPO, 'scanpy')
LOG_DIR = os.path.join(REPO, 'scripts', 'logs')

# sample -> batch key.  Every sample is normalised to the SAME total (TARGET_SUM);
# the original notebook instead used hard-coded per-condition "total RNA" scalings
# (exp = 3802, SHX = 3507, casp = 46), which imposed a ~76x scale gap between
# conditions before log1p.  That is deliberately not done here.
SAMPLES = [
    ('sample_2b',  'exp'),
    ('sample_13a', 'dis1'),
    ('sample_15a', 'dis2'),
    ('sample_13b', 'reg1'),
    ('sample_15b', 'reg2'),
]

TARGET_SUM = 1e4  # uniform across all samples (override with --target-sum)

# reporter/rRNA/tRNA-like genes dropped from every sample before normalisation
# (see docstring point (1) for why removal must happen before, not after)
GENES_TO_REMOVE = ['16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR', 'mCherry']
MIN_CELLS = 3         # scanpy filter_genes: drop genes detected in fewer cells than this
N_TOP_GENES = 2000     # highly-variable-gene subset size passed to --n-top-genes
N_COMPS = 50           # number of PCA components computed by sc.tl.pca
N_NEIGHBORS = 40       # neighbourhood size for the kNN graph feeding UMAP/Leiden
N_PCS = 40             # number of PCs used to build that neighbour graph
MIN_DIST = 0.3         # UMAP min_dist (point spread) parameter
RANDOM_STATE = 0       # UMAP RNG seed, for a reproducible embedding
LEIDEN_RESOLUTION = 0.4  # Leiden clustering resolution (higher = more clusters)
# all of the above are kept identical to scanpy_analysis.ipynb (see module docstring)

POINT_SIZE = 4       # figure only; scanpy's default would be ~24 at this n
POINT_ALPHA = 0.35

TS = datetime.now().strftime('%Y%m%d_%H%M%S')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)


def parse_args():
    # CLI overrides for the three deliberate differences from the notebook
    # (target-sum, barcode set, HVG count) plus one convenience switch (mCherry).
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target-sum', type=float, default=TARGET_SUM,
                        help='uniform normalize_total target for every sample '
                             f'(default: {TARGET_SUM:g})')
    parser.add_argument('--barcodes', choices=['paper', 'umap'], default='paper',
                        help='cell set: data_for_paper barcodes (default) or the '
                             'original data_for_umap barcodes')
    parser.add_argument('--n-top-genes', type=int, default=N_TOP_GENES,
                        help=f'HVGs to keep; 0 disables subsetting '
                             f'(default: {N_TOP_GENES}). The published run used no '
                             f'HVG subsetting.')
    parser.add_argument('--keep-mcherry', action='store_true',
                        help='keep the mCherry reporter, as the published run did '
                             '(default: remove it with the other reporters)')
    return parser.parse_args()


def load_sample(sample, barcodes_src, log):
    """Raw counts for one sample under the chosen cell set.

    Deliberately NOT normalised here: normalisation happens on the concatenated
    matrix after GENES_TO_REMOVE is dropped (see module docstring).
    """
    if barcodes_src == 'umap':
        mat = pd.read_csv(os.path.join(UMAP_DIR, f'{sample}_filtered.csv'), index_col=0)
    else:
        barcodes = pd.read_csv(
            os.path.join(PAPER_DIR, f'{sample}_filtered.csv'), index_col=0, usecols=[0]
        ).index
        mat = pd.read_csv(os.path.join(UNFILTERED_DIR, f'{sample}_unfiltered.csv'),
                          index_col=0)
        missing = barcodes.difference(mat.index)
        if len(missing):
            raise ValueError(f'{sample}: {len(missing)} paper barcodes absent from the '
                             f'unfiltered matrix (first few: {list(missing[:5])})')
        mat = mat.loc[barcodes]

    log[sample] = {
        'n_cells': int(mat.shape[0]),
        'n_genes_in': int(mat.shape[1]),
        'median_total_counts': float(mat.sum(1).median()),
    }
    print(f'{sample}: {mat.shape[0]} cells x {mat.shape[1]} genes')

    adata = sc.AnnData(mat.astype('float32'))
    adata.var_names = mat.columns.astype(str)
    adata.obs_names = mat.index.astype(str)
    return adata


def main():
    args = parse_args()
    target_sum = args.target_sum
    n_top_genes = args.n_top_genes
    src = f'{args.barcodes}_barcodes'

    # everything below is recorded so a given output set can be traced back to
    # exactly which code path and parameters produced it
    log = {'timestamp': TS,
           'script': os.path.basename(__file__),
           'barcode_set': args.barcodes,
           'barcode_source': PAPER_DIR if args.barcodes == 'paper' else UMAP_DIR,
           'matrix_source': UNFILTERED_DIR if args.barcodes == 'paper' else UMAP_DIR,
           'scanpy_version': sc.__version__,
           'params': {
               'order': 'remove genes -> filter_genes -> normalize_total -> log1p',
               'normalisation': 'uniform target_sum across all samples',
               'target_sum': target_sum,
               'genes_to_remove': GENES_TO_REMOVE,
               'keep_mcherry': bool(args.keep_mcherry),
               'filter_genes_min_cells': MIN_CELLS,
               'n_top_genes': n_top_genes or None,
               'pca_n_comps': N_COMPS, 'pca_svd_solver': 'arpack',
               'n_neighbors': N_NEIGHBORS, 'n_pcs': N_PCS,
               'umap_min_dist': MIN_DIST, 'random_state': RANDOM_STATE,
               'leiden_resolution': LEIDEN_RESOLUTION,
               'rank_genes_method': 'wilcoxon',
           },
           'samples': {}}

    adatas = [load_sample(s, args.barcodes, log['samples']) for s, _ in SAMPLES]
    keys = [k for _, k in SAMPLES]

    adata = sc.concat(adatas, label='batch', join='inner', fill_value=0,
                      keys=keys, index_unique='-')

    # Gene removal and filtering come FIRST, so that normalize_total divides by
    # a total that excludes 16s rRNA rather than being dominated by it.
    genes_to_remove = [g for g in GENES_TO_REMOVE
                       if not (g == 'mCherry' and args.keep_mcherry)]
    adata = adata[:, ~adata.var_names.isin(genes_to_remove)]
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    adata.layers['counts'] = adata.X.copy()   # raw counts, for marker testing
    print(f'{adata.n_obs} cells x {adata.n_vars} genes after gene removal '
          f'(target_sum={target_sum:g})')

    # normalise-then-log-transform on the gene-removed matrix (see docstring
    # point (1)): this is what makes the published UMAP reproducible here
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()  # normalised + log1p, all genes
    adata.raw = adata                          # full gene space, lognorm values
    if n_top_genes:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)
    sc.pp.scale(adata)  # z-score each gene, the standard input to PCA/UMAP
    sc.tl.pca(adata, svd_solver='arpack', n_comps=N_COMPS)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    sc.tl.umap(adata, min_dist=MIN_DIST, random_state=RANDOM_STATE)
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION)
    # marker genes are tested on raw counts (layer='counts'), not the scaled/
    # HVG-subset matrix used for the embedding, so results cover all genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon',
                            use_raw=False, layer='counts')

    log['n_cells_total'] = int(adata.n_obs)
    log['n_hvg'] = int(adata.n_vars)
    log['cells_per_batch'] = {k: int(v) for k, v in
                              adata.obs['batch'].value_counts().items()}
    log['cells_per_cluster'] = {k: int(v) for k, v in
                                adata.obs['leiden'].value_counts().items()}

    # ---- outputs ---------------------------------------------------------
    # X = scaled (z-scored) values; layers['counts'] = raw integer counts;
    # layers['lognorm'] = normalised + log1p; .raw = lognorm over all genes.
    h5ad_path = os.path.join(OUT_DIR, f'adata_shx_{src}_{TS}.h5ad')
    adata.write_h5ad(h5ad_path)

    umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'],
                           index=adata.obs_names)
    umap_df['batch'] = adata.obs['batch'].values
    umap_df['cluster'] = adata.obs['leiden'].values
    coord_path = os.path.join(OUT_DIR, f'umap_coordinates_shx_{src}_{TS}.csv')
    umap_df.to_csv(coord_path)

    # Small semi-transparent points, and a shuffled draw order so that no single
    # group is painted entirely on top of the others.
    shuffled = adata[np.random.default_rng(0).permutation(adata.n_obs)]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    sc.pl.umap(shuffled, color='leiden', ax=axes[0], show=False, frameon=False,
               title='Leiden clusters', size=POINT_SIZE, alpha=POINT_ALPHA)
    sc.pl.umap(shuffled, color='batch', ax=axes[1], show=False, frameon=False,
               title='Sample', size=POINT_SIZE, alpha=POINT_ALPHA)
    fig.tight_layout()
    svg_path = os.path.join(OUT_DIR, f'umap_shx_{src}_{TS}.svg')
    png_path = os.path.join(OUT_DIR, f'umap_shx_{src}_{TS}.png')
    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    markers_path = os.path.join(OUT_DIR, f'marker_genes_shx_{src}_{TS}.xlsx')
    with pd.ExcelWriter(markers_path) as writer:
        for cl in adata.obs['leiden'].cat.categories:
            sc.get.rank_genes_groups_df(adata, group=cl).to_excel(
                writer, sheet_name=f'cluster_{cl}', index=False)

    log['outputs'] = {'adata': h5ad_path, 'coordinates': coord_path, 'figure_svg': svg_path,
                      'figure_png': png_path, 'markers': markers_path}
    log_path = os.path.join(LOG_DIR, f'umap_{src}_{TS}.json')
    with open(log_path, 'w') as fh:
        json.dump(log, fh, indent=2)

    print(f'\n{adata.n_obs} cells, {adata.n_vars} HVGs, '
          f'{len(adata.obs["leiden"].cat.categories)} clusters')
    print('adata       ->', h5ad_path)
    print('coordinates ->', coord_path)
    print('figure      ->', svg_path)
    print('markers     ->', markers_path)
    print('log         ->', log_path)


if __name__ == '__main__':
    main()
