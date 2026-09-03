"""Does the between-sample GMP-Cor ordering survive matching for p, n and depth?

For each source (h5ad, data_for_paper) independently:
  1. restrict all 5 samples to the same gene set (intersection) -> identical p
  2. drop cells below the depth target, subsample to identical n
  3. multinomial-thin every cell to exactly TARGET_UMI counts -> identical depth
  4. compute GMP-Cor

After this, p, n and per-cell depth are identical across samples within a source,
so any remaining GMP-Cor difference reflects correlation structure rather than
matrix geometry or sequencing depth.

Repeated REPS times with different random draws to gauge stability.

Gene exclusion: 16s_mature, 16s_unprocessed, mCherry (tmRNA / tetR / YFP kept,
matching analysis_notebook.ipynb, which used t_genes=[]).
"""
import os
import sys
import json
import functools
import datetime

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import src.analysis_functions as af

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, 'data', 'scanpy_shx.h5ad')
PAPER = os.path.join(ROOT, 'data_for_paper')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')

DROP_GENES = ['16s_mature', '16s_unprocessed', 'mCherry']
BATCH_TO_FILE = {
    'exp':  'sample_2b_filtered.csv',
    'dis1': 'sample_13a_filtered.csv',
    'dis2': 'sample_15a_filtered.csv',
    'reg1': 'sample_13b_filtered.csv',
    'reg2': 'sample_15b_filtered.csv',
}
BATCH_ORDER = ['exp', 'dis1', 'dis2', 'reg1', 'reg2']

REPS = 3
# row-sum target used by af.get_eig_dist's own 'sum' normalization; since every row
# here is separately thinned to exactly target_umi counts first, this only rescales
# the already-uniform row sums and does not itself equalize depth across cells
NORM_SUM = 50
SEED = 0


def gmp_cor(m):
    """Compute GMP-Cor for a cells x genes count matrix m.

    Wraps af.get_eig_dist (empirical eigenvalues `pcs` vs. the mean scrambled/null
    eigenvalues `pcs1`) and reduces it to the same scalar order parameter used
    throughout the paper: the summed excess of empirical eigenvalues above the
    scrambled noise threshold (max eigenvalue of the permuted-gene null).

    Returns (gmp_cor, lambda_max, lambda_max_scrambled, n_modes_above_threshold).
    """
    pcs, pcs1, _ = af.get_eig_dist(m, norm=True, log=False,
                                   norm_method='sum', norm_sum=NORM_SUM)
    thr = float(pcs1.max())
    return (float(np.sum(np.maximum(pcs - thr, 0))), float(pcs.max()), thr,
            int((pcs > thr).sum()))


def thin(M, target, rng):
    """Multinomial-thin each row of M down to exactly `target` counts."""
    out = np.zeros_like(M)
    for i in range(M.shape[0]):
        row = M[i]
        tot = row.sum()
        if tot <= target:
            # row already at/below target: caller is expected to have pre-filtered
            # rows below target_umi, so this branch is a defensive no-op, not the
            # intended path
            out[i] = row
        else:
            out[i] = rng.multinomial(target, row / tot)
    return out


def load_source(name):
    """Return {batch: (matrix, gene_names)} for one source: raw (unfiltered-for-p,
    unthinned) cells x genes count matrices, read either from the combined h5ad
    (split by the 'batch' obs column) or from the per-sample data_for_paper CSVs
    (genes x cells on disk, transposed here via .values on the cells-indexed frame).
    DROP_GENES is removed from both so the two sources start from the same gene
    universe before the common-gene intersection in main().
    """
    if name == 'h5ad':
        adata = ad.read_h5ad(H5AD)
        C = adata.layers['counts']
        C = (C.toarray() if sp.issparse(C) else np.asarray(C)).astype(float)
        genes = np.asarray(adata.var_names)
        keep = ~np.isin(genes, DROP_GENES)
        C, genes = C[:, keep], genes[keep]
        batch = adata.obs['batch'].astype(str).values
        return {b: (C[batch == b], genes) for b in BATCH_ORDER}
    out = {}
    for b, f in BATCH_TO_FILE.items():
        d = pd.read_csv(os.path.join(PAPER, f), index_col=0)
        d = d.drop(columns=[g for g in DROP_GENES if g in d.columns])
        out[b] = (d.values.astype(float), np.asarray(d.columns))
    return out


def main():
    """For each source, match p/n/depth across the 5 samples and compute GMP-Cor
    REPS times per sample; write the per-run results as CSV/JSON and a plain-text
    summary table under results/cluster_gmp_cor/."""
    os.makedirs(OUTDIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    records = []
    setup = {}

    for source in ['h5ad', 'paper']:
        data = load_source(source)

        # 1. common gene set
        common = functools.reduce(lambda a, b: a & b,
                                  (set(g) for _, g in data.values()))
        common = sorted(common)
        mats = {}
        for b, (M, g) in data.items():
            ix = pd.Index(g).get_indexer(common)
            mats[b] = M[:, ix]
        p = len(common)

        # 2/3. depth and cell-count targets: use the tightest sample
        # target_umi = the lowest per-sample median depth: every sample must have a
        # substantial fraction of cells at/above this to be thinned down to it, and
        # taking the min guarantees no sample is asked to be thinned *up*
        med = {b: float(np.median(M.sum(1))) for b, M in mats.items()}
        target_umi = int(np.floor(min(med.values())))
        n_avail = {b: int((M.sum(1) >= target_umi).sum()) for b, M in mats.items()}
        # target_n = the scarcest sample's eligible-cell count, so every sample can
        # supply exactly target_n cells at/above target_umi without replacement
        target_n = min(n_avail.values())
        setup[source] = {
            'n_common_genes': p, 'target_umi': target_umi, 'target_n_cells': target_n,
            'median_depth_per_sample': med, 'cells_available_per_sample': n_avail,
        }
        print(f'\n=== {source}: p={p}, target depth={target_umi} UMI, '
              f'target n={target_n} cells ===')
        print(f'  median depth per sample: '
              f'{ {k: round(v, 1) for k, v in med.items()} }')
        print(f'  cells at/above target:   {n_avail}')

        for rep in range(REPS):
            # SEED + rep: each rep gets its own reproducible draw (cell subsample +
            # multinomial thinning), so repeated runs of this script are deterministic
            # while still sampling REPS independent realizations for the stability check
            rng = np.random.default_rng(SEED + rep)
            for b in BATCH_ORDER:
                M = mats[b]
                elig = np.where(M.sum(1) >= target_umi)[0]
                pick = rng.choice(elig, target_n, replace=False)
                sub = thin(M[pick], target_umi, rng)
                g, lam, thr, nmodes = gmp_cor(sub)
                records.append({
                    'source': source, 'batch': b, 'rep': rep,
                    'n_cells': int(sub.shape[0]), 'n_genes': p,
                    'target_umi': target_umi,
                    'mean_total_expression': float(sub.sum(1).mean()),
                    'mean_genes_detected': float((sub > 0).sum(1).mean()),
                    'lambda_max': lam, 'lambda_max_scrambled': thr,
                    'n_modes_above_threshold': nmodes, 'gmp_cor': g,
                })
                print(f'  rep{rep} {b:5s}: GMP-Cor = {g:8.3f}  '
                      f'(lam {lam:6.2f} vs scr {thr:5.2f}, {nmodes:3d} modes, '
                      f'{records[-1]["mean_genes_detected"]:.1f} genes/cell)')

    df = pd.DataFrame(records)
    summ = (df.groupby(['source', 'batch'])['gmp_cor']
              .agg(['mean', 'std']).reset_index()
              .pivot(index='batch', columns='source').loc[BATCH_ORDER])
    print('\n=== matched GMP-Cor (mean +/- sd over %d reps) ===' % REPS)
    print(summ.to_string(float_format=lambda v: '%.2f' % v))

    det = (df.groupby(['source', 'batch'])['mean_genes_detected'].mean()
             .reset_index().pivot(index='batch', columns='source').loc[BATCH_ORDER])
    print('\n=== mean genes detected per cell after thinning ===')
    print(det.to_string(float_format=lambda v: '%.1f' % v))

    csv_path = os.path.join(OUTDIR, f'matched_downsampling_{stamp}.csv')
    df.to_csv(csv_path, index=False)
    txt_path = os.path.join(OUTDIR, f'matched_downsampling_{stamp}.txt')
    with open(txt_path, 'w') as fh:
        fh.write('Matched downsampling: equal p, n and per-cell depth within source\n')
        fh.write(f'generated {stamp}\ndropped genes: {DROP_GENES}\nreps: {REPS}\n\n')
        for s, v in setup.items():
            fh.write(f'{s}: {json.dumps(v)}\n')
        fh.write('\n' + summ.to_string(float_format=lambda v: '%.2f' % v) + '\n\n')
        fh.write(det.to_string(float_format=lambda v: '%.1f' % v) + '\n')
    json_path = os.path.join(OUTDIR, f'matched_downsampling_{stamp}.json')
    with open(json_path, 'w') as fh:
        json.dump({'params': {'h5ad': H5AD, 'paper_dir': PAPER,
                              'dropped_genes': DROP_GENES, 'reps': REPS,
                              'norm_sum': NORM_SUM, 'seed': SEED,
                              'thinning': 'multinomial to exact target UMI',
                              'gmp_cor_definition':
                                  'sum(max(lambda_i - max_scrambled_lambda, 0))'},
                   'setup': setup, 'results': records}, fh, indent=2)
    print(f'\nwrote:\n  {csv_path}\n  {txt_path}\n  {json_path}')


if __name__ == '__main__':
    main()
