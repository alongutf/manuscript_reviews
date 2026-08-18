"""GMP-Cor per sample x leiden cluster, sourced entirely from data/scanpy_shx.h5ad.

Counts come from layers['counts'] (raw integers). The h5ad gene set already has
rRNA / spike-in / reporter probes removed ('16s_mature', '16s_unprocessed',
'LELOBEKK', 'kanR'), so no barcode matching against data_for_* is needed and
every cell carries its own leiden label. All samples share one 4042-gene set,
so GMP-Cor is directly comparable across groups.

GMP-Cor = sum_i max(lambda_i - lambda_max^scrambled, 0)

Groups with fewer than MIN_CELLS cells are skipped.
"""
import os
import sys
import json
import datetime

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import src.analysis_functions as af

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, 'data', 'scanpy_shx.h5ad')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')

MIN_CELLS = 100
NORM = True
LOG = False
NORM_METHOD = 'sum'
NORM_SUM = 50          # matches analysis_notebook.ipynb; z-transform makes this scale-free
SEED = 0

# provenance only - which per-sample matrix each h5ad batch was built from
BATCH_TO_FILE = {
    'exp':  'sample_2b_filtered.csv',
    'dis1': 'sample_13a_filtered.csv',
    'dis2': 'sample_15a_filtered.csv',
    'reg1': 'sample_13b_filtered.csv',
    'reg2': 'sample_15b_filtered.csv',
}
BATCH_ORDER = ['exp', 'dis1', 'dis2', 'reg1', 'reg2']


def gmp_cor(m):
    """GMP-Cor and spectrum summary for a cells x genes count matrix."""
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
    C = adata.layers['counts']
    C = C.toarray() if sp.issparse(C) else np.asarray(C)
    C = C.astype(float)
    batch = adata.obs['batch'].astype(str).values
    leiden = adata.obs['leiden'].astype(str).values
    genes = np.asarray(adata.var_names)
    print(f'h5ad: {C.shape[0]} cells x {C.shape[1]} genes')

    # tmRNA fraction is a known composition confounder - track it per group
    tmrna_ix = int(np.where(genes == 'tmRNA')[0][0]) if 'tmRNA' in set(genes) else None

    records = []
    for b in BATCH_ORDER:
        print(f'\n{b} ({BATCH_TO_FILE[b]}): {int((batch == b).sum())} cells')
        for cl in sorted(set(leiden[batch == b]), key=int):
            sel = (batch == b) & (leiden == cl)
            n = int(sel.sum())
            if n < MIN_CELLS:
                print(f'  cluster {cl}: {n} cells -> skipped (<{MIN_CELLS})')
                records.append({'batch': b, 'file': BATCH_TO_FILE[b], 'leiden': cl,
                                'n_cells': n, 'skipped': True})
                continue

            M = C[sel]
            tot = M.sum(axis=1)
            rec = {
                'batch': b,
                'file': BATCH_TO_FILE[b],
                'leiden': cl,
                'n_cells': n,
                'skipped': False,
                'mean_total_expression': float(tot.mean()),
                'median_total_expression': float(np.median(tot)),
                'sd_total_expression': float(tot.std(ddof=1)),
                'n_genes': int(M.shape[1]),
            }
            if tmrna_ix is not None:
                rec['tmrna_pct_of_counts'] = float(100 * M[:, tmrna_ix].sum() / M.sum())
            rec.update(gmp_cor(M))
            print(f'  cluster {cl}: {n} cells, mean total expr '
                  f'{rec["mean_total_expression"]:.1f}, '
                  f'tmRNA {rec.get("tmrna_pct_of_counts", float("nan")):.1f}%, '
                  f'GMP-Cor = {rec["gmp_cor"]:.3f} '
                  f'(lambda_max {rec["lambda_max"]:.2f} vs scr {rec["lambda_max_scrambled"]:.2f})')
            records.append(rec)

    df = pd.DataFrame(records)
    kept = df[~df.skipped]
    cols = ['batch', 'leiden', 'n_cells', 'mean_total_expression',
            'tmrna_pct_of_counts', 'gmp_cor']

    csv_path = os.path.join(OUTDIR, f'cluster_gmp_cor_h5ad_{stamp}.csv')
    df.to_csv(csv_path, index=False)

    txt_path = os.path.join(OUTDIR, f'cluster_gmp_cor_h5ad_{stamp}.txt')
    with open(txt_path, 'w') as fh:
        fh.write('GMP-Cor per sample x leiden cluster (SHX), counts from scanpy_shx.h5ad\n')
        fh.write(f'generated {stamp}\n')
        fh.write(f'{C.shape[0]} cells x {C.shape[1]} genes; min group size {MIN_CELLS}\n\n')
        fh.write(kept[cols].to_string(index=False) + '\n')
    print('\n' + kept[cols].to_string(index=False))

    json_path = os.path.join(OUTDIR, f'cluster_gmp_cor_h5ad_{stamp}.json')
    with open(json_path, 'w') as fh:
        json.dump({
            'params': {
                'source': H5AD, 'layer': 'counts', 'min_cells': MIN_CELLS,
                'norm': NORM, 'log': LOG, 'norm_method': NORM_METHOD,
                'norm_sum': NORM_SUM, 'seed': SEED, 'scramble_reps': 10,
                'n_cells': int(C.shape[0]), 'n_genes': int(C.shape[1]),
                'batch_to_file': BATCH_TO_FILE,
                'gmp_cor_definition': 'sum(max(lambda_i - max_scrambled_lambda, 0))',
            },
            'results': records,
        }, fh, indent=2)

    print(f'\nwrote:\n  {csv_path}\n  {txt_path}\n  {json_path}')


if __name__ == '__main__':
    main()
