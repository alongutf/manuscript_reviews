"""GMP-Cor per sample x leiden cluster (SHX dataset).

Leiden labels come from data/scanpy_shx.h5ad (built from data_for_umap/).
Count matrices come from data_for_paper/. Cells are matched by barcode.

GMP-Cor = sum_i max(lambda_i - lambda_max^scrambled, 0)

Two gene sets are reported:
  native  - each sample's own genes (as in the per-sample paper analysis)
  common  - intersection of genes across the 5 samples, so that groups are
            compared at identical matrix width p

Groups with fewer than MIN_CELLS matched cells are skipped.
"""
import os
import sys
import json
import datetime

import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import src.analysis_functions as af

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, 'data', 'scanpy_shx.h5ad')
PAPER = os.path.join(ROOT, 'data_for_paper')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')

MIN_CELLS = 100
NORM = True
LOG = False
NORM_METHOD = 'sum'
NORM_SUM = 50          # matches analysis_notebook.ipynb; z-transform makes this scale-free
SEED = 0

# batch key in the h5ad -> per-sample count matrix in data_for_paper
BATCH_TO_FILE = {
    'exp':  'sample_2b_filtered.csv',
    'dis1': 'sample_13a_filtered.csv',
    'dis2': 'sample_15a_filtered.csv',
    'reg1': 'sample_13b_filtered.csv',
    'reg2': 'sample_15b_filtered.csv',
}


def gmp_cor(m):
    """GMP-Cor and the underlying spectrum summary for a cells x genes count matrix."""
    pcs, pcs1, frac_nz = af.get_eig_dist(
        m, norm=NORM, log=LOG, norm_method=NORM_METHOD, norm_sum=NORM_SUM
    )
    thr = float(pcs1.max())
    return {
        'gmp_cor': float(np.sum(np.maximum(pcs - thr, 0))),
        'lambda_max': float(pcs.max()),
        'lambda_max_scrambled': thr,
        'n_modes_above_threshold': int(np.sum(pcs > thr)),
        'fraction_non_zero': float(frac_nz),
    }


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    np.random.seed(SEED)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    adata = ad.read_h5ad(H5AD)
    # h5ad obs names look like "<barcode>-1-<batch>"
    obs = pd.DataFrame({
        'barcode': [n.split('-1-')[0] for n in adata.obs_names],
        'batch': adata.obs['batch'].astype(str).values,
        'leiden': adata.obs['leiden'].astype(str).values,
    })

    matrices = {}
    for batch, fname in BATCH_TO_FILE.items():
        d = pd.read_csv(os.path.join(PAPER, fname), index_col=0)
        d.index = [str(i).split('-')[0] for i in d.index]
        matrices[batch] = d

    common_genes = sorted(set.intersection(*(set(d.columns) for d in matrices.values())))
    print(f'common gene set: {len(common_genes)} genes')

    records = []
    coverage = []
    for batch, d in matrices.items():
        sub = obs[obs.batch == batch]
        labels = dict(zip(sub.barcode, sub.leiden))
        shared = [b for b in d.index if b in labels]
        coverage.append({
            'batch': batch,
            'file': BATCH_TO_FILE[batch],
            'n_cells_matrix': int(d.shape[0]),
            'n_cells_h5ad': int(sub.shape[0]),
            'n_matched': len(shared),
            'match_rate_vs_matrix': round(len(shared) / d.shape[0], 4),
        })
        print(f'\n{batch} ({BATCH_TO_FILE[batch]}): '
              f'{len(shared)}/{d.shape[0]} cells carry a leiden label')

        for cl in sorted(set(labels[b] for b in shared), key=int):
            cells = [b for b in shared if labels[b] == cl]
            if len(cells) < MIN_CELLS:
                print(f'  cluster {cl}: {len(cells)} cells -> skipped (<{MIN_CELLS})')
                records.append({
                    'batch': batch, 'file': BATCH_TO_FILE[batch], 'leiden': cl,
                    'n_cells': len(cells), 'skipped': True,
                })
                continue

            block = d.loc[cells]
            counts = block.values.astype(float)
            rec = {
                'batch': batch,
                'file': BATCH_TO_FILE[batch],
                'leiden': cl,
                'n_cells': len(cells),
                'skipped': False,
                'mean_total_expression': float(counts.sum(axis=1).mean()),
                'median_total_expression': float(np.median(counts.sum(axis=1))),
                'sd_total_expression': float(counts.sum(axis=1).std(ddof=1)),
                'n_genes_native': int(counts.shape[1]),
                'n_genes_common': len(common_genes),
            }
            print(f'  cluster {cl}: {len(cells)} cells, '
                  f'mean total expr {rec["mean_total_expression"]:.1f}')

            for tag, mat in (('native', counts),
                             ('common', block[common_genes].values.astype(float))):
                for k, v in gmp_cor(mat).items():
                    rec[f'{k}_{tag}'] = v
                print(f'    {tag}: GMP-Cor = {rec["gmp_cor_" + tag]:.3f} '
                      f'(lambda_max {rec["lambda_max_" + tag]:.2f} vs '
                      f'scr {rec["lambda_max_scrambled_" + tag]:.2f})')
            records.append(rec)

    df = pd.DataFrame(records)
    cov = pd.DataFrame(coverage)

    csv_path = os.path.join(OUTDIR, f'cluster_gmp_cor_{stamp}.csv')
    df.to_csv(csv_path, index=False)

    txt_path = os.path.join(OUTDIR, f'cluster_gmp_cor_{stamp}.txt')
    kept = df[~df.skipped]
    with open(txt_path, 'w') as fh:
        fh.write('GMP-Cor per sample x leiden cluster (SHX)\n')
        fh.write(f'generated {stamp}\n\n')
        fh.write('Barcode match rate, data_for_paper vs scanpy_shx.h5ad:\n')
        fh.write(cov.to_string(index=False) + '\n\n')
        fh.write(f'Groups with >= {MIN_CELLS} matched cells:\n')
        fh.write(kept[['batch', 'leiden', 'n_cells', 'mean_total_expression',
                       'gmp_cor_native', 'gmp_cor_common']].to_string(index=False) + '\n')
    print('\n' + kept[['batch', 'leiden', 'n_cells', 'mean_total_expression',
                       'gmp_cor_native', 'gmp_cor_common']].to_string(index=False))

    json_path = os.path.join(OUTDIR, f'cluster_gmp_cor_{stamp}.json')
    with open(json_path, 'w') as fh:
        json.dump({
            'params': {
                'h5ad': H5AD, 'data_dir': PAPER, 'min_cells': MIN_CELLS,
                'norm': NORM, 'log': LOG, 'norm_method': NORM_METHOD,
                'norm_sum': NORM_SUM, 'seed': SEED,
                'scramble_reps': 10,
                'batch_to_file': BATCH_TO_FILE,
                'n_common_genes': len(common_genes),
                'gmp_cor_definition': 'sum(max(lambda_i - max_scrambled_lambda, 0))',
            },
            'coverage': coverage,
            'common_genes': common_genes,
            'results': records,
        }, fh, indent=2)

    print(f'\nwrote:\n  {csv_path}\n  {txt_path}\n  {json_path}')


if __name__ == '__main__':
    main()
