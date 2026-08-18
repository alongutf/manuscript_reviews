"""equate_dims with depth matching instead of a uniform random draw.

The published matrices came from equate_dims, which takes a uniform random
subsample of target_cells from each pool. Because the eligible pools differ
enormously in size (13a: ~1600 cells, 15a: ~15600), the uniform draw reaches
much further down 15a's depth distribution, leaving the two matrices at very
different depths -- and GMP-Cor tracks depth.

This keeps everything else about the paper pipeline intact:
  * cell pool:   filter_by_umi_count(400, 20000) on TOTAL counts
  * gene panel:  filter_by_gene_dispersion(min_dispersion=1), then the
                 pairwise intersection within each condition pair
  * output size: ~1000 cells x ~2000 genes

and replaces only the cell draw: instead of uniform sampling, cells are drawn so
that every sample ends up with the SAME depth histogram (depth measured as the
row sum on that sample's final gene panel). Bins take min(counts) across all
samples, then are scaled to the target cell number.

Reports how far the resulting barcode sets move from the published files.
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
PAPER = os.path.join(ROOT, 'data_for_paper')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')

UNFILT = {
    'dis1': 'sample_13a_unfiltered.csv',
    'dis2': 'sample_15a_unfiltered.csv',
    'reg1': 'sample_13b_unfiltered.csv',
    'reg2': 'sample_15b_unfiltered.csv',
}
PUBLISHED = {
    'dis1': 'sample_13a_filtered.csv',
    'dis2': 'sample_15a_filtered.csv',
    'reg1': 'sample_13b_filtered.csv',
    'reg2': 'sample_15b_filtered.csv',
}
PAIRS = [('dis1', 'dis2'), ('reg1', 'reg2')]
ORDER = ['dis1', 'dis2', 'reg1', 'reg2']

UMI_MIN, UMI_MAX = 400, 20000     # as in analysis_notebook cell 3
MIN_DISPERSION = 1.0
TARGET_CELLS = 1000
N_BINS = 40
REPS = 5
CHUNK = 20000
NORM_SUM = 50
SEED = 0


def load_eligible(sample):
    """One streaming pass: keep rows with UMI_MIN < total < UMI_MAX."""
    path = os.path.join(DATA, UNFILT[sample])
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    rows, bcs = [], []
    for ch in pd.read_csv(path, index_col=0, chunksize=CHUNK, dtype=dtypes):
        V = ch.values
        tot = V.sum(1)
        sel = (tot > UMI_MIN) & (tot < UMI_MAX)
        if sel.any():
            rows.append(V[sel])
            bcs.append(np.array([str(b).split('-')[0] for b in ch.index])[sel])
    M = np.vstack(rows).astype(float)
    return M, np.concatenate(bcs), genes


def dispersion_genes(M, genes):
    """filter_by_gene_dispersion(min_dispersion=1)."""
    mean = M.mean(axis=0)
    var = M.var(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        disp = np.where(mean > 0, var / mean, 0.0)
    return set(genes[np.nan_to_num(disp) > MIN_DISPERSION])


def matched_draw(depths, target, n_bins, rng):
    """Draw `target` indices per sample so all share one depth histogram."""
    alld = np.concatenate(list(depths.values()))
    edges = np.unique(np.quantile(alld, np.linspace(0, 1, n_bins + 1)))
    binned = {s: np.digitize(d, edges[1:-1]) for s, d in depths.items()}
    nb = len(edges) - 1

    avail = np.array([min(int((binned[s] == b).sum()) for s in depths)
                      for b in range(nb)])
    total_avail = int(avail.sum())
    if total_avail == 0:
        raise RuntimeError('no overlapping depth range across samples')

    take = avail if total_avail <= target else np.floor(
        avail * target / total_avail).astype(int)
    # distribute any remainder to the bins with most headroom
    while take.sum() < min(target, total_avail):
        head = avail - take
        take[int(np.argmax(head))] += 1

    picks = {}
    for s in depths:
        out = []
        for b in range(nb):
            idx = np.where(binned[s] == b)[0]
            k = int(take[b])
            if k > 0:
                out.append(rng.choice(idx, k, replace=False))
        picks[s] = np.concatenate(out) if out else np.array([], dtype=int)
    return picks, int(take.sum())


def gmp_cor(m):
    pcs, pcs1, _ = af.get_eig_dist(m, norm=True, log=False,
                                   norm_method='sum', norm_sum=NORM_SUM)
    thr = float(pcs1.max())
    return (float(np.sum(np.maximum(pcs - thr, 0))), float(pcs.max()), thr)


def main():
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    pool = {}
    for s in ORDER:
        M, bc, genes = load_eligible(s)
        pool[s] = {'M': M, 'bc': bc, 'genes': genes,
                   'disp': dispersion_genes(M, genes)}
        print(f'{s}: eligible {M.shape[0]} cells, '
              f'{len(pool[s]["disp"])} genes pass dispersion')

    # gene panel = pairwise intersection within each condition pair
    panel = {}
    for a, b in PAIRS:
        common = sorted(pool[a]['disp'] & pool[b]['disp'])
        panel[a] = panel[b] = common
        print(f'panel {a}+{b}: {len(common)} genes')

    # depth = row sum on that sample's final panel
    sub, depth = {}, {}
    for s in ORDER:
        ix = pd.Index(pool[s]['genes']).get_indexer(panel[s])
        sub[s] = pool[s]['M'][:, ix]
        depth[s] = sub[s].sum(1)
        print(f'{s}: pool depth on panel  median {np.median(depth[s]):7.1f}  '
              f'mean {depth[s].mean():7.1f}')

    published = {}
    for s in ORDER:
        idx = pd.read_csv(os.path.join(PAPER, PUBLISHED[s]),
                          index_col=0, usecols=[0]).index
        published[s] = set(str(i).split('-')[0] for i in idx)

    records, overlaps = [], []
    for rep in range(REPS):
        rng = np.random.default_rng(SEED + rep)
        picks, n_got = matched_draw(depth, TARGET_CELLS, N_BINS, rng)
        if rep == 0:
            print(f'\nmatched draw yields {n_got} cells per sample '
                  f'(target {TARGET_CELLS})')
        for s in ORDER:
            ix = picks[s]
            m = sub[s][ix]
            bcs = set(pool[s]['bc'][ix])
            g, lam, thr = gmp_cor(m)
            rec = {
                'sample': s, 'rep': rep, 'n_cells': int(m.shape[0]),
                'n_genes': len(panel[s]),
                'mean_total_expression': float(m.sum(1).mean()),
                'median_total_expression': float(np.median(m.sum(1))),
                'mean_genes_detected': float((m > 0).sum(1).mean()),
                'lambda_max': lam, 'lambda_max_scrambled': thr, 'gmp_cor': g,
            }
            records.append(rec)
            overlaps.append({
                'sample': s, 'rep': rep,
                'n_new': len(bcs),
                'n_published': len(published[s]),
                'shared_with_published': len(bcs & published[s]),
                'pct_of_published_retained':
                    round(100 * len(bcs & published[s]) / len(published[s]), 1),
            })
            if rep == 0:
                print(f'  {s}: n={m.shape[0]} p={len(panel[s])} '
                      f'depth mean {rec["mean_total_expression"]:.1f} '
                      f'GMP-Cor {g:.3f}  '
                      f'| shares {overlaps[-1]["shared_with_published"]} '
                      f'barcodes with published '
                      f'({overlaps[-1]["pct_of_published_retained"]}%)')

    df = pd.DataFrame(records)
    ov = pd.DataFrame(overlaps)

    print('\n=== depth after matching (mean +/- sd over %d reps) ===' % REPS)
    print(df.groupby('sample')[['mean_total_expression', 'mean_genes_detected']]
            .agg(['mean', 'std']).loc[ORDER]
            .to_string(float_format=lambda v: '%.2f' % v))

    print('\n=== GMP-Cor (mean +/- sd) ===')
    g = df.groupby('sample')['gmp_cor'].agg(['mean', 'std', 'min', 'max']).loc[ORDER]
    print(g.to_string(float_format=lambda v: '%.3f' % v))

    print('\n=== barcode change vs data_for_paper ===')
    o = ov.groupby('sample')[['shared_with_published',
                              'pct_of_published_retained']].mean().loc[ORDER]
    print(o.to_string(float_format=lambda v: '%.1f' % v))

    d1, d2 = g.loc['dis1', 'mean'], g.loc['dis2', 'mean']
    r1, r2 = g.loc['reg1', 'mean'], g.loc['reg2', 'mean']
    print(f'\n  dis mean {np.mean([d1, d2]):.3f} | reg mean {np.mean([r1, r2]):.3f}')
    print(f'  within-dis gap {abs(d1-d2):.3f} | within-reg gap {abs(r1-r2):.3f} '
          f'| between-condition gap {abs(np.mean([d1,d2])-np.mean([r1,r2])):.3f}')

    df.to_csv(os.path.join(OUTDIR, f'equate_depth_matched_{stamp}.csv'), index=False)
    ov.to_csv(os.path.join(OUTDIR, f'equate_depth_overlap_{stamp}.csv'), index=False)
    with open(os.path.join(OUTDIR, f'equate_depth_matched_{stamp}.json'), 'w') as fh:
        json.dump({'params': {'umi_min': UMI_MIN, 'umi_max': UMI_MAX,
                              'min_dispersion': MIN_DISPERSION,
                              'target_cells': TARGET_CELLS, 'n_bins': N_BINS,
                              'reps': REPS, 'norm_sum': NORM_SUM, 'seed': SEED,
                              'pairs': PAIRS},
                   'results': records, 'overlaps': overlaps}, fh, indent=2)
    print(f'\nwrote results/cluster_gmp_cor/equate_depth_matched_{stamp}.*')


if __name__ == '__main__':
    main()
