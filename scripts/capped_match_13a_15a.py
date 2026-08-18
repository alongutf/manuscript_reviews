"""Cap 15a at 13a's deepest cell, then take the top 1000 below that cap.

Rationale: 13a (dis1) has a much smaller usable pool than 15a (dis2), so a
uniform draw from each (what equate_dims did) leaves 15a shallower. Capping 15a
at 13a's maximum and taking the 1000 cells immediately below removes 15a's
extreme tail before selection.

Selection is done on mRNA counts (rRNA/tRNA dropped first) as the primary
variant, with a total-count variant reported alongside since the original
pipeline thresholded on totals.

Both samples are then put on one gene set (detected in >= DETECTION_FRAC of
retained cells in BOTH samples) so GMP-Cor is comparable between them.
"""
import os
import sys
import json
import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import src.analysis_functions as af

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
DATA = os.path.join(OTHER, 'data')
BIOTYPE = os.path.join(OTHER, 'filtered_data', 'k12_biotype_map.csv')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE_DIR = os.path.join(OUTDIR, 'barcode_stats')

FILES = {'dis1': 'sample_13a_unfiltered.csv', 'dis2': 'sample_15a_unfiltered.csv'}
N_CELLS = 1000
DETECTION_FRAC = 0.05
CHUNK = 20000
NORM_SUM = 50
SEED = 0


def protein_coding_mask(genes):
    bt = pd.read_csv(BIOTYPE)
    pc = bt.gene[(bt.biotype != 'tRNA') & (bt.biotype != 'rRNA')].astype(str)
    pc = set(v.casefold() for v in pc)
    names = [str(v).casefold().replace('lelobekk_', '') for v in genes]
    return np.array([v in pc for v in names])


def load_rows(sample, wanted):
    """Stream the unfiltered CSV, keeping only barcodes in `wanted`."""
    path = os.path.join(DATA, FILES[sample])
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    pc = protein_coding_mask(genes)
    keep_rows, keep_bc = [], []
    for ch in pd.read_csv(path, index_col=0, chunksize=CHUNK, dtype=dtypes):
        bc = np.array([str(b).split('-')[0] for b in ch.index])
        sel = np.isin(bc, list(wanted))
        if sel.any():
            keep_rows.append(ch.values[sel][:, pc])
            keep_bc.append(bc[sel])
    return np.vstack(keep_rows).astype(float), np.concatenate(keep_bc), genes[pc]


def gmp_cor(m):
    pcs, pcs1, _ = af.get_eig_dist(m, norm=True, log=False,
                                   norm_method='sum', norm_sum=NORM_SUM)
    thr = float(pcs1.max())
    return {
        'gmp_cor': float(np.sum(np.maximum(pcs - thr, 0))),
        'lambda_max': float(pcs.max()),
        'lambda_max_scrambled': thr,
        'n_modes_above_threshold': int((pcs > thr).sum()),
    }


def main():
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    np.random.seed(SEED)

    st = {s: dict(np.load(os.path.join(CACHE_DIR, f'{s}.npz'), allow_pickle=True))
          for s in FILES}

    results = {}
    for metric in ['mrna']:
        d1, d2 = st['dis1'][metric], st['dis2'][metric]
        cap = int(d1.max())
        print(f'\n=== selection on {metric} ===')
        print(f'13a deepest cell: {cap} counts  -> cap for 15a')

        # 13a: its own top N
        o1 = np.argsort(-d1)[:N_CELLS]
        # 15a: top N among cells strictly below the cap
        elig = np.where(d2 < cap)[0]
        o2 = elig[np.argsort(-d2[elig])[:N_CELLS]]
        print(f'13a top {N_CELLS}: depth {d1[o1].min():.0f}-{d1[o1].max():.0f}, '
              f'median {np.median(d1[o1]):.0f}')
        print(f'15a top {N_CELLS} below cap: depth {d2[o2].min():.0f}-{d2[o2].max():.0f}, '
              f'median {np.median(d2[o2]):.0f}')
        print(f'15a cells at/above cap excluded: {int((d2 >= cap).sum())}')

        sel = {'dis1': set(st['dis1']['barcode'][o1]),
               'dis2': set(st['dis2']['barcode'][o2])}

        mats = {}
        for s in FILES:
            M, bc, genes = load_rows(s, sel[s])
            mats[s] = (M, genes)
            print(f'  loaded {s}: {M.shape}')

        # common gene set: detected in >=DETECTION_FRAC of cells in BOTH
        keep = None
        for s, (M, genes) in mats.items():
            frac = (M > 0).mean(axis=0)
            g = set(genes[frac >= DETECTION_FRAC])
            keep = g if keep is None else (keep & g)
        keep = sorted(keep)
        print(f'  common gene set: {len(keep)} genes')

        for s, (M, genes) in mats.items():
            ix = pd.Index(genes).get_indexer(keep)
            sub = M[:, ix]
            tot = sub.sum(1)
            rec = {
                'sample': s, 'metric': metric,
                'n_cells': int(sub.shape[0]), 'n_genes': len(keep),
                'mean_total_expression': float(tot.mean()),
                'median_total_expression': float(np.median(tot)),
                'mean_genes_detected': float((sub > 0).sum(1).mean()),
            }
            rec.update(gmp_cor(sub))
            results[f'{s}_{metric}'] = rec
            print(f'  {s}: mean tot {rec["mean_total_expression"]:.1f}, '
                  f'{rec["mean_genes_detected"]:.1f} genes/cell, '
                  f'GMP-Cor = {rec["gmp_cor"]:.3f} '
                  f'(lam {rec["lambda_max"]:.2f} vs scr {rec["lambda_max_scrambled"]:.2f})')

    df = pd.DataFrame(results).T
    cols = ['sample', 'n_cells', 'n_genes', 'mean_total_expression',
            'mean_genes_detected', 'lambda_max', 'lambda_max_scrambled', 'gmp_cor']
    print('\n=== comparison ===')
    print(df[cols].to_string(index=False))

    df.to_csv(os.path.join(OUTDIR, f'capped_match_{stamp}.csv'), index=False)
    with open(os.path.join(OUTDIR, f'capped_match_{stamp}.json'), 'w') as fh:
        json.dump({'params': {'n_cells': N_CELLS, 'detection_frac': DETECTION_FRAC,
                              'norm_sum': NORM_SUM, 'seed': SEED},
                   'results': results}, fh, indent=2)
    print(f'\nwrote results/cluster_gmp_cor/capped_match_{stamp}.*')


if __name__ == '__main__':
    main()
