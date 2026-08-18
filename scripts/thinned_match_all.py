"""Depth-matched dis vs reg: equal n, equal p, identical per-cell depth.

Same construction as thinned_match_13a_15a.py, extended to all samples so the
condition contrast can be read with the technical confounds removed:

  1. keep cells with mRNA >= T in every sample (rRNA/tRNA dropped first)
  2. subsample all samples to the same n (binding sample sets it)
  3. one gene set: detected in >= DETECTION_FRAC of retained cells in EVERY sample
  4. multinomial-thin every cell to exactly T counts, on that gene set
  5. GMP-Cor, repeated REPS times

Conditions: dis1/dis2 = dysregulated, reg1/reg2 = regulated, exp = exponential
reference. Replicates are reported individually -- with n=2 per condition the
replicate spread is the relevant yardstick for any condition difference.
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

FILES = {
    'exp':  'sample_2b_unfiltered.csv',
    'dis1': 'sample_13a_unfiltered.csv',
    'dis2': 'sample_15a_unfiltered.csv',
    'reg1': 'sample_13b_unfiltered.csv',
    'reg2': 'sample_15b_unfiltered.csv',
}
CONDITION = {'exp': 'exponential', 'dis1': 'dysregulated', 'dis2': 'dysregulated',
             'reg1': 'regulated', 'reg2': 'regulated'}
ORDER = ['exp', 'dis1', 'dis2', 'reg1', 'reg2']

T = 75
DETECTION_FRAC = 0.05
REPS = 10
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
    path = os.path.join(DATA, FILES[sample])
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    pc = protein_coding_mask(genes)
    rows, bcs = [], []
    for ch in pd.read_csv(path, index_col=0, chunksize=CHUNK, dtype=dtypes):
        bc = np.array([str(b).split('-')[0] for b in ch.index])
        sel = np.isin(bc, list(wanted))
        if sel.any():
            rows.append(ch.values[sel][:, pc])
            bcs.append(bc[sel])
    return np.vstack(rows).astype(float), np.concatenate(bcs), genes[pc]


def thin(M, target, rng):
    out = np.zeros_like(M)
    for i in range(M.shape[0]):
        row = M[i]
        tot = row.sum()
        out[i] = row if tot <= target else rng.multinomial(target, row / tot)
    return out


def gmp_cor(m):
    pcs, pcs1, _ = af.get_eig_dist(m, norm=True, log=False,
                                   norm_method='sum', norm_sum=NORM_SUM)
    thr = float(pcs1.max())
    return (float(np.sum(np.maximum(pcs - thr, 0))), float(pcs.max()), thr,
            int((pcs > thr).sum()))


def main():
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    st = {s: dict(np.load(os.path.join(CACHE_DIR, f'{s}.npz'), allow_pickle=True))
          for s in ORDER}

    elig = {s: st[s]['barcode'][st[s]['mrna'] >= T] for s in ORDER}
    n = min(len(v) for v in elig.values())
    binding = min(elig, key=lambda s: len(elig[s]))
    print(f'T={T}: eligible ' + ', '.join(f'{s}={len(elig[s])}' for s in ORDER))
    print(f'binding sample: {binding} -> n={n} per sample')

    mats = {}
    for s in ORDER:
        M, bc, genes = load_rows(s, set(elig[s]))
        mats[s] = (M, genes)
        print(f'  loaded {s}: {M.shape}')

    keep = None
    for s in ORDER:
        M, genes = mats[s]
        frac = (M > 0).mean(axis=0)
        g = set(genes[frac >= DETECTION_FRAC])
        keep = g if keep is None else (keep & g)
    keep = sorted(keep)
    print(f'  common gene set: {len(keep)} genes (gamma = {len(keep)/n:.2f})')

    records = []
    for rep in range(REPS):
        rng = np.random.default_rng(SEED + rep)
        for s in ORDER:
            M, genes = mats[s]
            ix = pd.Index(genes).get_indexer(keep)
            sub = M[:, ix]
            pick = rng.choice(sub.shape[0], n, replace=False)
            sub = thin(sub[pick], T, rng)
            g, lam, thr, nm = gmp_cor(sub)
            records.append({
                'sample': s, 'condition': CONDITION[s], 'rep': rep,
                'n_cells': n, 'n_genes': len(keep), 'target_umi': T,
                'mean_total_expression': float(sub.sum(1).mean()),
                'mean_genes_detected': float((sub > 0).sum(1).mean()),
                'lambda_max': lam, 'lambda_max_scrambled': thr,
                'n_modes_above_threshold': nm, 'gmp_cor': g,
            })
            print(f'  rep{rep} {s:5s}: GMP-Cor = {g:7.3f} '
                  f'({records[-1]["mean_genes_detected"]:.1f} genes/cell)')

    df = pd.DataFrame(records)

    per_sample = df.groupby('sample')[['mean_total_expression',
                                       'mean_genes_detected', 'gmp_cor']] \
                   .agg(['mean', 'std']).loc[ORDER]
    print('\n=== per sample (mean +/- sd over %d reps) ===' % REPS)
    print(per_sample.to_string(float_format=lambda v: '%.3f' % v))

    print('\n=== per sample GMP-Cor range ===')
    for s in ORDER:
        v = df[df['sample'] == s]['gmp_cor'].values
        print(f'  {s:5s} ({CONDITION[s]:13s}): {v.mean():6.3f}  '
              f'range {v.min():6.3f} - {v.max():6.3f}')

    print('\n=== condition contrast ===')
    dis = df[df['condition'] == 'dysregulated']['gmp_cor']
    reg = df[df['condition'] == 'regulated']['gmp_cor']
    print(f'  dysregulated: {dis.mean():.3f} (range {dis.min():.3f}-{dis.max():.3f})')
    print(f'  regulated:    {reg.mean():.3f} (range {reg.min():.3f}-{reg.max():.3f})')
    print(f'  reg / dis ratio: {reg.mean()/dis.mean():.2f}')

    # replicate spread vs condition spread -- the decisive comparison at n=2
    d1 = df[df['sample'] == 'dis1']['gmp_cor'].mean()
    d2 = df[df['sample'] == 'dis2']['gmp_cor'].mean()
    r1 = df[df['sample'] == 'reg1']['gmp_cor'].mean()
    r2 = df[df['sample'] == 'reg2']['gmp_cor'].mean()
    print(f'\n  within-dis replicate gap: {abs(d1-d2):.3f}')
    print(f'  within-reg replicate gap: {abs(r1-r2):.3f}')
    print(f'  between-condition gap:    {abs((d1+d2)/2 - (r1+r2)/2):.3f}')
    print('  -> condition difference is only interpretable if it exceeds '
          'the replicate gaps')

    df.to_csv(os.path.join(OUTDIR, f'thinned_match_all_{stamp}.csv'), index=False)
    with open(os.path.join(OUTDIR, f'thinned_match_all_{stamp}.json'), 'w') as fh:
        json.dump({'params': {'T': T, 'n_cells': n, 'n_genes': len(keep),
                              'reps': REPS, 'detection_frac': DETECTION_FRAC,
                              'norm_sum': NORM_SUM, 'seed': SEED,
                              'binding_sample': binding, 'condition': CONDITION},
                   'results': records}, fh, indent=2)
    print(f'\nwrote results/cluster_gmp_cor/thinned_match_all_{stamp}.*')


if __name__ == '__main__':
    main()
