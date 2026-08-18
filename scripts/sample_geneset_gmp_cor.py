"""GMP-Cor per sample under different 2000-gene selections (all 5 SHX samples).

Source: data/scanpy_shx.h5ad, layers['counts'] (raw integers). The h5ad already
excludes '16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR'; here we drop the
remaining high-abundance non-endogenous / reporter probes as well.

All cells of each sample are used (no cluster split). Gene selection is done
independently within each sample.

Gene sets compared (2000 genes each):
  max_expr   - top 2000 by max count across cells
  mean_expr  - top 2000 by mean count across cells (sanity check; max is tie-heavy)
  fano       - top 2000 by Fano factor (var/mean) among genes passing a detection floor

GMP-Cor = sum_i max(lambda_i - lambda_max^scrambled, 0)
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

BATCHES = ['exp', 'dis1', 'dis2', 'reg1', 'reg2']
N_GENES = 2000
DROP_GENES = ['tmRNA', 'tetR', 'mCherry']
# Fano is unstable for genes seen in a handful of cells; require a detection floor
FANO_MIN_CELLS = 10

NORM = True
LOG = False
NORM_METHOD = 'sum'
NORM_SUM = 50          # z-transform makes this scale-free
SEED = 0


def gmp_cor(m):
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


def select_genes(M):
    """Return ({selection_name: gene index array}, per-gene stats) for one sample."""
    mean_expr = M.mean(axis=0)
    max_expr = M.max(axis=0)
    var_expr = M.var(axis=0, ddof=1)
    n_det = (M > 0).sum(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        fano = np.where(mean_expr > 0, var_expr / mean_expr, 0.0)
    fano_rank = np.where(n_det >= FANO_MIN_CELLS, fano, -np.inf)
    sel = {
        'max_expr':  np.sort(np.argsort(-max_expr, kind='stable')[:N_GENES]),
        'mean_expr': np.sort(np.argsort(-mean_expr, kind='stable')[:N_GENES]),
        'fano':      np.sort(np.argsort(-fano_rank, kind='stable')[:N_GENES]),
    }
    return sel, {'mean': mean_expr, 'fano': fano, 'n_det': n_det}


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    np.random.seed(SEED)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    adata = ad.read_h5ad(H5AD)
    C = adata.layers['counts']
    C = (C.toarray() if sp.issparse(C) else np.asarray(C)).astype(float)
    genes = np.asarray(adata.var_names)
    batch = adata.obs['batch'].astype(str).values

    keep = ~np.isin(genes, DROP_GENES)
    dropped = [g for g in DROP_GENES if g in set(genes)]
    C = C[:, keep]
    genes = genes[keep]
    print(f'dropped {dropped}; working gene space: {C.shape[1]} genes')

    records = []
    overlaps_all = {}
    gene_sets_all = {}
    for b in BATCHES:
        M = C[batch == b]
        print(f'\n=== {b}: {M.shape[0]} cells x {M.shape[1]} genes ===')
        sel, st = select_genes(M)
        gene_sets_all[b] = {k: genes[v].tolist() for k, v in sel.items()}

        keys = list(sel)
        ov = {}
        for i, ka in enumerate(keys):
            for kb in keys[i + 1:]:
                ov[f'{ka}|{kb}'] = int(len(set(sel[ka].tolist()) & set(sel[kb].tolist())))
        overlaps_all[b] = ov
        print(f'  gene-set overlaps: {ov}')

        for name, ix in sel.items():
            sub = M[:, ix]
            tot = sub.sum(axis=1)
            rec = {
                'batch': b,
                'selection': name,
                'n_cells': int(sub.shape[0]),
                'n_genes': int(sub.shape[1]),
                'mean_total_expression': float(tot.mean()),
                'median_total_expression': float(np.median(tot)),
                'pct_of_all_counts_retained': float(100 * sub.sum() / M.sum()),
                'mean_gene_mean': float(st['mean'][ix].mean()),
                'mean_gene_fano': float(st['fano'][ix].mean()),
                'median_cells_detected': float(np.median(st['n_det'][ix])),
            }
            rec.update(gmp_cor(sub))
            records.append(rec)
            print(f'  [{name:9s}] mean tot {rec["mean_total_expression"]:7.1f}, '
                  f'{rec["pct_of_all_counts_retained"]:5.1f}% counts, '
                  f'Fano {rec["mean_gene_fano"]:.3f}, '
                  f'GMP-Cor = {rec["gmp_cor"]:8.3f} '
                  f'(lam {rec["lambda_max"]:6.2f} vs scr {rec["lambda_max_scrambled"]:5.2f})')

    df = pd.DataFrame(records)

    pivot = df.pivot(index='batch', columns='selection', values='gmp_cor').loc[BATCHES]
    print('\n=== GMP-Cor by sample x gene selection ===')
    print(pivot.to_string(float_format=lambda v: '%.2f' % v))
    print(f'\n=== gene-set overlap (of {N_GENES}) ===')
    print(pd.DataFrame(overlaps_all).T.loc[BATCHES].to_string())

    csv_path = os.path.join(OUTDIR, f'sample_geneset_gmp_cor_{stamp}.csv')
    df.to_csv(csv_path, index=False)

    cols = ['batch', 'selection', 'n_cells', 'mean_total_expression',
            'pct_of_all_counts_retained', 'mean_gene_fano',
            'lambda_max', 'lambda_max_scrambled', 'gmp_cor']
    txt_path = os.path.join(OUTDIR, f'sample_geneset_gmp_cor_{stamp}.txt')
    with open(txt_path, 'w') as fh:
        fh.write(f'GMP-Cor per sample, {N_GENES}-gene selections\n')
        fh.write(f'generated {stamp}\nsource {H5AD} layers["counts"]\n')
        fh.write(f'dropped genes: {dropped}\n\n')
        fh.write(df[cols].to_string(index=False) + '\n\n')
        fh.write('GMP-Cor by sample x selection:\n')
        fh.write(pivot.to_string(float_format=lambda v: '%.2f' % v) + '\n\n')
        fh.write(f'gene-set overlaps (of {N_GENES}):\n')
        fh.write(pd.DataFrame(overlaps_all).T.loc[BATCHES].to_string() + '\n')

    json_path = os.path.join(OUTDIR, f'sample_geneset_gmp_cor_{stamp}.json')
    with open(json_path, 'w') as fh:
        json.dump({
            'params': {
                'source': H5AD, 'layer': 'counts', 'batches': BATCHES,
                'n_genes': N_GENES, 'dropped_genes': dropped,
                'fano_min_cells': FANO_MIN_CELLS,
                'norm': NORM, 'log': LOG, 'norm_method': NORM_METHOD,
                'norm_sum': NORM_SUM, 'seed': SEED, 'scramble_reps': 10,
                'gmp_cor_definition': 'sum(max(lambda_i - max_scrambled_lambda, 0))',
            },
            'gene_set_overlaps': overlaps_all,
            'gene_sets': gene_sets_all,
            'results': records,
        }, fh, indent=2)

    print(f'\nwrote:\n  {csv_path}\n  {txt_path}\n  {json_path}')


if __name__ == '__main__':
    main()
