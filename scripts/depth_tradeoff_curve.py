"""Step 2 feasibility check: how many cells survive at each common depth target?

For a fair GMP-Cor comparison every sample must contribute cells at one fixed
mRNA depth T. This script computes, for each sample and each candidate T, the
number of barcodes with mRNA_counts >= T, and reports the binding (minimum)
count across samples -- i.e. the usable n per sample if T were chosen.

rRNA/tRNA are dropped BEFORE counting, so depth means mRNA depth. 16s is ~92%
of raw counts, so any threshold on total counts is really a threshold on rRNA.

Per-barcode stats are cached to .npz so later steps need no further full passes.
"""
import os
import sys
import datetime

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
DATA = os.path.join(OTHER, 'data')
BIOTYPE = os.path.join(OTHER, 'filtered_data', 'k12_biotype_map.csv')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE_DIR = os.path.join(OUTDIR, 'barcode_stats')

SAMPLES = {
    'exp':  'sample_2b_unfiltered.csv',
    'dis1': 'sample_13a_unfiltered.csv',
    'dis2': 'sample_15a_unfiltered.csv',
    'reg1': 'sample_13b_unfiltered.csv',
    'reg2': 'sample_15b_unfiltered.csv',
}
TARGETS = [25, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000]
CHUNK = 20000
# a gene must be detected in this fraction of retained cells to be usable
DETECTION_FRAC = 0.05


def protein_coding_mask(genes):
    bt = pd.read_csv(BIOTYPE)
    pc = bt.gene[(bt.biotype != 'tRNA') & (bt.biotype != 'rRNA')].astype(str)
    pc = set(v.casefold() for v in pc)
    names = [str(v).casefold().replace('lelobekk_', '') for v in genes]
    return np.array([v in pc for v in names])


def build_stats(sample, fname):
    """Per-barcode mRNA depth + per-gene detection counts, cached."""
    cache = os.path.join(CACHE_DIR, f'{sample}.npz')
    if os.path.exists(cache):
        print(f'  {sample}: cache hit')
        return dict(np.load(cache, allow_pickle=True))

    path = os.path.join(DATA, fname)
    print(f'  {sample}: scanning {fname} ...')
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    pc = protein_coding_mask(genes)
    pc_genes = genes[pc]

    bcs, mrna, det = [], [], []
    # per-gene detection counts, accumulated per depth target
    det_by_T = {T: np.zeros(pc.sum(), dtype=np.int64) for T in TARGETS}
    for i, ch in enumerate(pd.read_csv(path, index_col=0, chunksize=CHUNK,
                                       dtype=dtypes)):
        V = ch.values[:, pc]
        m = V.sum(1)
        bcs.append(np.array([str(b).split('-')[0] for b in ch.index]))
        mrna.append(m)
        det.append((V > 0).sum(1))
        for T in TARGETS:
            sel = m >= T
            if sel.any():
                det_by_T[T] += (V[sel] > 0).sum(0)
    out = {
        'barcode': np.concatenate(bcs),
        'mrna': np.concatenate(mrna),
        'detected': np.concatenate(det),
        'pc_genes': pc_genes,
    }
    for T in TARGETS:
        out[f'genedet_{T}'] = det_by_T[T]
    os.makedirs(CACHE_DIR, exist_ok=True)
    np.savez_compressed(cache, **out)
    return out


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    stats = {}
    for s, f in SAMPLES.items():
        stats[s] = build_stats(s, f)

    print('\n=== per-sample mRNA depth distribution (all barcodes) ===')
    for s in SAMPLES:
        m = stats[s]['mrna']
        print(f'  {s:5s} n={len(m):6d}  median {np.median(m):6.0f}  '
              f'p90 {np.percentile(m, 90):7.0f}  p99 {np.percentile(m, 99):8.0f}  '
              f'max {m.max():8.0f}')

    rows = []
    for T in TARGETS:
        r = {'T': T}
        for s in SAMPLES:
            r[s] = int((stats[s]['mrna'] >= T).sum())
        r['min_across_samples'] = min(r[s] for s in SAMPLES)
        # genes usable at this depth: detected in >=DETECTION_FRAC of retained
        # cells in EVERY sample
        usable = None
        for s in SAMPLES:
            n = max(int((stats[s]['mrna'] >= T).sum()), 1)
            frac = stats[s][f'genedet_{T}'] / n
            keep = set(stats[s]['pc_genes'][frac >= DETECTION_FRAC])
            usable = keep if usable is None else (usable & keep)
        r['genes_usable'] = len(usable)
        rows.append(r)

    df = pd.DataFrame(rows)
    print(f'\n=== cells surviving at each mRNA depth target T ===')
    print(df.to_string(index=False))

    df['n_x_T'] = df['min_across_samples'] * df['T']
    best = df.loc[df['n_x_T'].idxmax()]
    print(f"\nlargest n*T product at T={int(best['T'])}: "
          f"n={int(best['min_across_samples'])} per sample, "
          f"{int(best['genes_usable'])} usable genes")

    csv_path = os.path.join(OUTDIR, f'depth_tradeoff_{stamp}.csv')
    df.to_csv(csv_path, index=False)
    txt_path = os.path.join(OUTDIR, f'depth_tradeoff_{stamp}.txt')
    with open(txt_path, 'w') as fh:
        fh.write('Depth / cell-count tradeoff, mRNA counts (rRNA+tRNA removed)\n')
        fh.write(f'generated {stamp}\ndetection floor: {DETECTION_FRAC}\n\n')
        fh.write(df.to_string(index=False) + '\n')
    print(f'\nwrote:\n  {csv_path}\n  {txt_path}')


if __name__ == '__main__':
    main()
