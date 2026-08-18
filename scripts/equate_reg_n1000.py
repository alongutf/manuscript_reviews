"""Raise the reg pair to n=1000 by lowering umi_min, then compare against dis.

At filter_by_umi_count(400, 20000) the reg pool caps at 313 cells (reg1), which
forced n=285 and left reg ~3x deeper than dis -- an uncontrolled confound in the
cross-condition comparison. Lowering umi_min admits shallower cells, which both
raises n to 1000 and pulls reg's depth and detection down toward dis.

The reg pools are re-scanned with a permissive floor (UMI_FLOOR) and cached, so
any threshold at or above that floor can be explored without another CSV pass.
umi_min is then chosen as the largest value still yielding >= TARGET_CELLS in
both reg samples.

dis is unchanged (its cached pool already uses umi_min=400) and is re-reported
for comparison.
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
POOLDIR = os.path.join(OUTDIR, 'eligible_pools')

UNFILT = {'reg1': 'sample_13b_unfiltered.csv', 'reg2': 'sample_15b_unfiltered.csv'}
PUBLISHED = {'dis1': 'sample_13a_filtered.csv', 'dis2': 'sample_15a_filtered.csv',
             'reg1': 'sample_13b_filtered.csv', 'reg2': 'sample_15b_filtered.csv'}
ORDER = ['dis1', 'dis2', 'reg1', 'reg2']

UMI_FLOOR = 50          # permissive floor for the re-scan / cache
UMI_MAX = 20000
MIN_DISPERSION = 1.0
TARGET_CELLS = 1000
N_BINS = 40
REPS = 5
CHUNK = 20000
NORM_SUM = 50
SEED = 0
DROP_GENES = ['16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR', 'mCherry']


def scan_low(sample):
    """Re-scan with the permissive floor and cache."""
    cache = os.path.join(POOLDIR, f'{sample}_low.npz')
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True)
        print(f'{sample}: low-floor cache hit {d["M"].shape}')
        return d['M'].astype(float), d['bc'], d['genes']
    path = os.path.join(DATA, UNFILT[sample])
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    rows, bcs = [], []
    for ch in pd.read_csv(path, index_col=0, chunksize=CHUNK, dtype=dtypes):
        V = ch.values
        tot = V.sum(1)
        sel = (tot > UMI_FLOOR) & (tot < UMI_MAX)
        if sel.any():
            rows.append(V[sel])
            bcs.append(np.array([str(b).split('-')[0] for b in ch.index])[sel])
    M = np.vstack(rows)
    bc = np.concatenate(bcs)
    np.savez_compressed(cache, M=M, bc=bc, genes=genes)
    print(f'{sample}: cached {M.shape} at floor {UMI_FLOOR}')
    return M.astype(float), bc, genes


def load_cached(sample):
    d = np.load(os.path.join(POOLDIR, f'{sample}.npz'), allow_pickle=True)
    return d['M'].astype(float), d['bc'], d['genes']


def dispersion_genes(M, genes):
    mean = M.mean(axis=0)
    var = M.var(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        disp = np.where(mean > 0, var / mean, 0.0)
    return set(genes[np.nan_to_num(disp) > MIN_DISPERSION])


def matched_draw(vals, target, n_bins, rng):
    allv = np.concatenate(list(vals.values()))
    edges = np.unique(np.quantile(allv, np.linspace(0, 1, n_bins + 1)))
    binned = {s: np.digitize(v, edges[1:-1]) for s, v in vals.items()}
    nb = len(edges) - 1
    avail = np.array([min(int((binned[s] == b).sum()) for s in vals)
                      for b in range(nb)])
    total = int(avail.sum())
    take = avail if total <= target else np.floor(avail * target / total).astype(int)
    while take.sum() < min(target, total):
        take[int(np.argmax(avail - take))] += 1
    picks = {}
    for s in vals:
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

    raw = {}
    for s in ['reg1', 'reg2']:
        raw[s] = scan_low(s)
    for s in ['dis1', 'dis2']:
        raw[s] = load_cached(s)

    # choose the largest umi_min still giving TARGET_CELLS in both reg samples
    tot = {s: raw[s][0].sum(1) for s in ['reg1', 'reg2']}
    chosen = UMI_FLOOR
    for cand in range(UMI_FLOOR, 401, 5):
        if all(int((tot[s] > cand).sum()) >= TARGET_CELLS for s in ['reg1', 'reg2']):
            chosen = cand
    print(f'\nreg umi_min chosen: {chosen} '
          f'(reg1 {int((tot["reg1"] > chosen).sum())} cells, '
          f'reg2 {int((tot["reg2"] > chosen).sum())} cells)')
    print(f'  for reference at umi_min=400: '
          f'reg1 {int((tot["reg1"] > 400).sum())}, reg2 {int((tot["reg2"] > 400).sum())}')

    pool = {}
    for s in ORDER:
        M, bc, genes = raw[s]
        if s.startswith('reg'):
            sel = M.sum(1) > chosen
            M, bc = M[sel], bc[sel]
        pool[s] = {'M': M, 'bc': bc, 'genes': genes,
                   'disp': dispersion_genes(M, genes)}
        print(f'  {s}: pool {M.shape[0]} cells')

    panel = {}
    for a, b in [('dis1', 'dis2'), ('reg1', 'reg2')]:
        common = sorted((pool[a]['disp'] & pool[b]['disp']) - set(DROP_GENES))
        panel[a] = panel[b] = common
        print(f'  panel {a}+{b}: {len(common)} genes')

    sub, det = {}, {}
    for s in ORDER:
        ix = pd.Index(pool[s]['genes']).get_indexer(panel[s])
        sub[s] = pool[s]['M'][:, ix]
        det[s] = (sub[s] > 0).sum(1)

    published = {}
    for s in ORDER:
        idx = pd.read_csv(os.path.join(PAPER, PUBLISHED[s]),
                          index_col=0, usecols=[0]).index
        published[s] = set(str(i).split('-')[0] for i in idx)

    records, overlaps = [], []
    for rep in range(REPS):
        rng = np.random.default_rng(SEED + rep)
        for a, b in [('dis1', 'dis2'), ('reg1', 'reg2')]:
            picks, n_got = matched_draw({a: det[a], b: det[b]},
                                        TARGET_CELLS, N_BINS, rng)
            if rep == 0:
                print(f'\npair {a}+{b}: {n_got} cells each')
            for s in (a, b):
                ix = picks[s]
                m = sub[s][ix]
                bcs = set(pool[s]['bc'][ix])
                g, lam, thr = gmp_cor(m)
                records.append({
                    'sample': s, 'rep': rep, 'n_cells': int(m.shape[0]),
                    'n_genes': len(panel[s]),
                    'mean_genes_detected': float((m > 0).sum(1).mean()),
                    'mean_depth': float(m.sum(1).mean()),
                    'median_depth': float(np.median(m.sum(1))),
                    'lambda_max': lam, 'lambda_max_scrambled': thr, 'gmp_cor': g,
                })
                overlaps.append({
                    'sample': s, 'rep': rep,
                    'shared_with_published': len(bcs & published[s]),
                    'pct_of_published_retained':
                        round(100 * len(bcs & published[s]) / len(published[s]), 1),
                })
                if rep == 0:
                    r = records[-1]
                    print(f'  {s}: detected {r["mean_genes_detected"]:.1f}  '
                          f'depth {r["mean_depth"]:.1f}  GMP-Cor {g:.3f}')

    df = pd.DataFrame(records)
    ov = pd.DataFrame(overlaps)

    print('\n=== detection / depth / GMP-Cor ===')
    print(df.groupby('sample')[['n_cells', 'n_genes', 'mean_genes_detected',
                                'mean_depth', 'gmp_cor']]
            .agg(['mean', 'std']).loc[ORDER]
            .to_string(float_format=lambda v: '%.2f' % v))

    print('\n=== barcode change vs data_for_paper ===')
    print(ov.groupby('sample')[['shared_with_published',
                                'pct_of_published_retained']]
            .mean().loc[ORDER].to_string(float_format=lambda v: '%.1f' % v))

    g = df.groupby('sample')['gmp_cor'].mean()
    d1, d2, r1, r2 = g['dis1'], g['dis2'], g['reg1'], g['reg2']
    print(f'\n  dis mean {np.mean([d1,d2]):.3f} | reg mean {np.mean([r1,r2]):.3f}')
    print(f'  within-dis gap {abs(d1-d2):.3f} | within-reg gap {abs(r1-r2):.3f} '
          f'| between-condition gap {abs(np.mean([d1,d2])-np.mean([r1,r2])):.3f}')

    df.to_csv(os.path.join(OUTDIR, f'reg_n1000_{stamp}.csv'), index=False)
    ov.to_csv(os.path.join(OUTDIR, f'reg_n1000_overlap_{stamp}.csv'), index=False)
    with open(os.path.join(OUTDIR, f'reg_n1000_{stamp}.json'), 'w') as fh:
        json.dump({'params': {'umi_floor': UMI_FLOOR, 'reg_umi_min': chosen,
                              'dis_umi_min': 400, 'umi_max': UMI_MAX,
                              'min_dispersion': MIN_DISPERSION,
                              'target_cells': TARGET_CELLS, 'reps': REPS,
                              'drop_genes': DROP_GENES, 'norm_sum': NORM_SUM,
                              'seed': SEED},
                   'results': records, 'overlaps': overlaps}, fh, indent=2)
    print(f'\nwrote results/cluster_gmp_cor/reg_n1000_{stamp}.*')


if __name__ == '__main__':
    main()
