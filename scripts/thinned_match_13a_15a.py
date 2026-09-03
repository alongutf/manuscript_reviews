"""Depth-matched 13a vs 15a: equal n, equal p, and identical per-cell depth.

The capped-selection variant left 15a ~2x deeper than 13a, and GMP-Cor tracked
that ratio almost 1:1. Here depth is made identical by construction:

  1. keep cells with mRNA >= T in both samples
  2. subsample both to the same n
  3. multinomial-thin every cell to exactly T counts
  4. one gene set, detected in >= DETECTION_FRAC of retained cells in BOTH

Repeated REPS times so the difference can be read against its own noise.

Reads cached per-barcode mRNA totals from results/cluster_gmp_cor/barcode_stats/
(<sample>.npz, produced by an earlier eigenvector/QC pass) to decide which barcodes
are eligible before re-reading the raw per-cell count CSVs. Writes a per-repeat
results table and a JSON log of the run to results/cluster_gmp_cor/. Run directly
with `python thinned_match_13a_15a.py`; no command-line arguments.
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
# hard-coded absolute path to a sibling repository holding the raw per-cell CSVs
# and the biotype map; this script only runs on a machine with that repo checked
# out at this exact location
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
DATA = os.path.join(OTHER, 'data')
BIOTYPE = os.path.join(OTHER, 'filtered_data', 'k12_biotype_map.csv')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE_DIR = os.path.join(OUTDIR, 'barcode_stats')

FILES = {'dis1': 'sample_13a_unfiltered.csv', 'dis2': 'sample_15a_unfiltered.csv'}
T = 75                  # common mRNA depth floor and thinning target
DETECTION_FRAC = 0.05
REPS = 10
CHUNK = 20000            # rows read per pandas chunk when scanning the raw CSVs
NORM_SUM = 50
SEED = 0                 # base RNG seed; each repeat uses SEED + rep for independence


def protein_coding_mask(genes):
    """Boolean mask, in `genes` order, of probes that are protein-coding (not
    tRNA/rRNA) per the K-12 biotype map. Probe names carry a fixed prefix that is
    stripped before matching against the biotype table's gene names."""
    bt = pd.read_csv(BIOTYPE)
    pc = bt.gene[(bt.biotype != 'tRNA') & (bt.biotype != 'rRNA')].astype(str)
    pc = set(v.casefold() for v in pc)
    names = [str(v).casefold().replace('lelobekk_', '') for v in genes]
    return np.array([v in pc for v in names])


def load_rows(sample, wanted):
    """Stream-read the raw per-cell count CSV for `sample`, keeping only rows whose
    barcode (the part of the index before the first '-') is in `wanted`, and only
    protein-coding gene columns. Returns (counts [cells x genes], barcodes, gene
    names) — read in chunks of CHUNK rows so the full unfiltered file, which is
    much larger than the barcodes actually needed, is never held in memory at once.
    """
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
    """Multinomial-downsample each row (cell) of `M` [cells x genes] to exactly
    `target` total counts, redistributing its counts proportionally at random. Rows
    already at or below `target` are left untouched (thinning never adds counts),
    which is why the eligibility filter upstream requires mRNA >= T for every cell:
    it guarantees every retained row actually gets thinned to the same depth."""
    out = np.zeros_like(M)
    for i in range(M.shape[0]):
        row = M[i]
        tot = row.sum()
        out[i] = row if tot <= target else rng.multinomial(target, row / tot)
    return out


def gmp_cor(m):
    """GMP-Cor and its components for one cells x genes count matrix `m`: sum-
    normalise (norm_sum counts/cell), compute the empirical and scrambled
    eigenvalue spectra via af.get_eig_dist, and use the scrambled spectrum's
    largest eigenvalue as the null threshold. Returns
    (GMP-Cor, empirical lambda_max, threshold, number of eigenvalues above
    threshold)."""
    pcs, pcs1, _ = af.get_eig_dist(m, norm=True, log=False,
                                   norm_method='sum', norm_sum=NORM_SUM)
    thr = float(pcs1.max())
    return (float(np.sum(np.maximum(pcs - thr, 0))), float(pcs.max()), thr,
            int((pcs > thr).sum()))


def main():
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    # cached per-barcode total mRNA counts from an earlier QC pass, used only to
    # decide eligibility (mRNA >= T) without re-scanning the full raw CSVs
    st = {s: dict(np.load(os.path.join(CACHE_DIR, f'{s}.npz'), allow_pickle=True))
          for s in FILES}

    elig = {s: st[s]['barcode'][st[s]['mrna'] >= T] for s in FILES}
    n = min(len(v) for v in elig.values())
    print(f'T={T}: eligible cells ' + ', '.join(f'{s}={len(v)}' for s, v in elig.items()))
    print(f'subsampling both to n={n}')

    mats = {}
    for s in FILES:
        M, bc, genes = load_rows(s, set(elig[s]))
        mats[s] = (M, genes)
        print(f'  loaded {s}: {M.shape}')

    # gene set common to both samples: detected (any count) in at least
    # DETECTION_FRAC of eligible cells, intersected across the two samples, so both
    # matrices end up with identical columns
    keep = None
    for s, (M, genes) in mats.items():
        frac = (M > 0).mean(axis=0)
        g = set(genes[frac >= DETECTION_FRAC])
        keep = g if keep is None else (keep & g)
    keep = sorted(keep)
    print(f'  common gene set: {len(keep)} genes (gamma = {len(keep)/n:.2f})')

    # repeat the whole subsample -> thin -> GMP-Cor pipeline REPS times with
    # independent RNGs, so the 13a-vs-15a difference can be judged against its own
    # repeat-to-repeat spread rather than a single draw
    records = []
    for rep in range(REPS):
        rng = np.random.default_rng(SEED + rep)
        for s, (M, genes) in mats.items():
            ix = pd.Index(genes).get_indexer(keep)
            sub = M[:, ix]
            pick = rng.choice(sub.shape[0], n, replace=False)
            sub = thin(sub[pick], T, rng)
            g, lam, thr, nm = gmp_cor(sub)
            records.append({
                'sample': s, 'rep': rep, 'n_cells': n, 'n_genes': len(keep),
                'target_umi': T,
                'mean_total_expression': float(sub.sum(1).mean()),
                'mean_genes_detected': float((sub > 0).sum(1).mean()),
                'lambda_max': lam, 'lambda_max_scrambled': thr,
                'n_modes_above_threshold': nm, 'gmp_cor': g,
            })
            print(f'  rep{rep} {s}: GMP-Cor = {g:7.3f} (lam {lam:5.2f} vs scr {thr:5.2f}, '
                  f'{records[-1]["mean_genes_detected"]:.1f} genes/cell)')

    df = pd.DataFrame(records)
    summ = df.groupby('sample')[['mean_total_expression', 'mean_genes_detected',
                                 'gmp_cor']].agg(['mean', 'std'])
    print('\n=== depth-matched comparison (mean +/- sd over %d reps) ===' % REPS)
    print(summ.to_string(float_format=lambda v: '%.3f' % v))

    a = df[df['sample'] == 'dis1']['gmp_cor'].values
    b = df[df['sample'] == 'dis2']['gmp_cor'].values
    print(f'\ndis2 / dis1 GMP-Cor ratio: {b.mean()/a.mean():.2f}')
    print(f'dis1 range {a.min():.2f}-{a.max():.2f} | dis2 range {b.min():.2f}-{b.max():.2f}')
    print(f'ranges overlap: {not (a.max() < b.min() or b.max() < a.min())}')

    df.to_csv(os.path.join(OUTDIR, f'thinned_match_{stamp}.csv'), index=False)
    with open(os.path.join(OUTDIR, f'thinned_match_{stamp}.json'), 'w') as fh:
        json.dump({'params': {'T': T, 'n_cells': n, 'n_genes': len(keep),
                              'reps': REPS, 'detection_frac': DETECTION_FRAC,
                              'norm_sum': NORM_SUM, 'seed': SEED},
                   'results': records}, fh, indent=2)
    print(f'\nwrote results/cluster_gmp_cor/thinned_match_{stamp}.*')


if __name__ == '__main__':
    main()
