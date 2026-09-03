"""Why does GMP-Cor differ so much between data sources?

Compares, per sample, the h5ad-derived count matrix against the corresponding
data_for_paper matrix, with a matched gene-exclusion rule:
    remove 16s_mature, 16s_unprocessed, mCherry   (tmRNA / tetR / YFP retained)

data_for_paper already lacks 16s_* / mCherry / kanR, so in practice only the
h5ad loses mCherry. Both sources then retain tmRNA, tetR, YFP.

For each matrix we report the structural statistics that drive an eigenvalue
spectrum -- shape, sparsity, depth, dominance of individual genes -- alongside
GMP-Cor, so the spread in GMP-Cor can be attributed.

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

NORM = True
LOG = False
NORM_METHOD = 'sum'
NORM_SUM = 50    # per-cell target total after normalization, before eigen-analysis
SEED = 0         # fixes the RNG used by get_eig_dist's scrambled-matrix repeats


def spectrum(m):
    """GMP-Cor plus the spectral quantities needed to explain it.

    m : cells x genes count/expression matrix. Returns a dict with GMP-Cor
    itself (`gmp_cor`, the total eigenvalue excess above the scrambled
    noise ceiling) and diagnostics for *why* GMP-Cor came out the way it
    did: whether it is dominated by one mode (`top_mode_share_of_gmp`) and
    how far the bulk (non-signal) eigenvalues sit below threshold
    (`bulk_mean_eigenvalue`).
    """
    pcs, pcs1, frac_nz = af.get_eig_dist(
        m, norm=NORM, log=LOG, norm_method=NORM_METHOD, norm_sum=NORM_SUM
    )
    thr = float(pcs1.max())                  # scrambled noise ceiling (lambda_max^scr)
    excess = np.maximum(pcs - thr, 0)         # per-mode signal above that ceiling
    total = float(excess.sum())               # GMP-Cor itself
    return {
        'gmp_cor': total,
        'lambda_max': float(pcs.max()),
        'lambda_max_scrambled': thr,
        'n_modes_above_threshold': int((pcs > thr).sum()),
        # fraction of GMP-Cor carried by the single largest mode - a value near 1
        # means GMP-Cor is really "one big eigenvalue", not a broad spectrum of signal
        'top_mode_share_of_gmp': float(excess.max() / total) if total > 0 else np.nan,
        # mean of the sub-threshold ("noise") eigenvalues, for reference against thr
        'bulk_mean_eigenvalue': float(pcs[pcs <= thr].mean()) if (pcs <= thr).any() else np.nan,
        'fraction_non_zero_after_filter': float(frac_nz),
    }


def describe(M, genes, label):
    """Structural statistics of a cells x genes count matrix, independent of
    the eigenvalue analysis, used to explain differences in GMP-Cor across
    matrices built from different sources/pipelines (depth, sparsity, gene
    dominance, etc. all shape the resulting spectrum)."""
    tot = M.sum(axis=1)
    det = (M > 0).sum(axis=1)
    gmean = M.mean(axis=0)
    gvar = M.var(axis=0, ddof=1)   # ddof=1: sample variance across cells, per gene
    with np.errstate(divide='ignore', invalid='ignore'):
        # Fano factor (variance/mean) per gene, a standard overdispersion measure;
        # guarded against divide-by-zero for genes with zero mean (set to 0 instead)
        fano = np.where(gmean > 0, gvar / gmean, 0.0)
    gsum = M.sum(axis=0)
    order = np.argsort(-gsum)   # genes ranked by total counts, descending
    rec = {
        'dataset': label,
        'n_cells': int(M.shape[0]),
        'n_genes': int(M.shape[1]),
        'aspect_p_over_n': float(M.shape[1] / M.shape[0]),
        'mean_total_expression': float(tot.mean()),
        'median_total_expression': float(np.median(tot)),
        'cv_total_expression': float(tot.std(ddof=1) / tot.mean()),
        'mean_genes_detected': float(det.mean()),
        'frac_nonzero': float((M > 0).mean()),
        'mean_gene_mean': float(gmean.mean()),
        'mean_gene_fano': float(fano.mean()),
        'genes_all_zero': int((gsum == 0).sum()),
        'top1_gene': str(genes[order[0]]),
        'top1_pct_counts': float(100 * gsum[order[0]] / M.sum()),
        'top5_pct_counts': float(100 * gsum[order[:5]].sum() / M.sum()),
        'top20_pct_counts': float(100 * gsum[order[:20]].sum() / M.sum()),
    }
    return rec


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    np.random.seed(SEED)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    adata = ad.read_h5ad(H5AD)
    C = adata.layers['counts']
    C = (C.toarray() if sp.issparse(C) else np.asarray(C)).astype(float)   # densify
    hgenes = np.asarray(adata.var_names)
    batch = adata.obs['batch'].astype(str).values
    keepc = ~np.isin(hgenes, DROP_GENES)
    print(f'h5ad: dropped {[g for g in DROP_GENES if g in set(hgenes)]}')
    C, hgenes = C[:, keepc], hgenes[keepc]

    records = []
    for b in BATCH_ORDER:
        fname = BATCH_TO_FILE[b]

        # --- h5ad source (scanpy_shx.h5ad, counts layer, this batch's cells) ---
        M = C[batch == b]
        rec = describe(M, hgenes, f'{b}:h5ad')
        rec.update({'batch': b, 'source': 'h5ad'})
        rec.update(spectrum(M))
        records.append(rec)

        # --- data_for_paper source (published, already-filtered CSV) ---
        d = pd.read_csv(os.path.join(PAPER, fname), index_col=0)
        d = d.drop(columns=[g for g in DROP_GENES if g in d.columns])
        P = d.values.astype(float)
        rec = describe(P, np.asarray(d.columns), f'{b}:paper')
        rec.update({'batch': b, 'source': 'paper'})
        rec.update(spectrum(P))
        records.append(rec)

        # shared context
        h_only = len(set(hgenes) - set(d.columns))
        p_only = len(set(d.columns) - set(hgenes))
        shared = len(set(hgenes) & set(d.columns))
        print(f'{b}: genes h5ad-only {h_only}, paper-only {p_only}, shared {shared}')

    df = pd.DataFrame(records)

    struct = ['dataset', 'n_cells', 'n_genes', 'aspect_p_over_n', 'mean_total_expression',
              'cv_total_expression', 'mean_genes_detected', 'frac_nonzero',
              'top1_gene', 'top1_pct_counts', 'top5_pct_counts']
    spec = ['dataset', 'frac_nonzero', 'lambda_max', 'lambda_max_scrambled',
            'n_modes_above_threshold', 'top_mode_share_of_gmp', 'gmp_cor']

    fmt = lambda v: '%.3f' % v
    print('\n=== structural statistics ===')
    print(df[struct].to_string(index=False, float_format=fmt))
    print('\n=== spectrum / GMP-Cor ===')
    print(df[spec].to_string(index=False, float_format=fmt))

    piv = df.pivot(index='batch', columns='source', values='gmp_cor').loc[BATCH_ORDER]
    piv['ratio_h5ad_over_paper'] = piv['h5ad'] / piv['paper']
    print('\n=== GMP-Cor by source ===')
    print(piv.to_string(float_format=fmt))

    # rank correlation (not Pearson, since these statistics have very different
    # scales/distributions) between every numeric statistic and GMP-Cor, across all
    # 2*len(BATCH_ORDER) matrices - a cheap way to see which structural property
    # tracks the GMP-Cor spread across sources
    num = df.select_dtypes(include=[np.number])
    corr = num.corr(method='spearman')['gmp_cor'].drop('gmp_cor').sort_values(key=abs, ascending=False)
    print('\n=== Spearman corr of each statistic with GMP-Cor (n=%d matrices) ===' % len(df))
    print(corr.to_string(float_format=fmt))

    csv_path = os.path.join(OUTDIR, f'source_comparison_{stamp}.csv')
    df.to_csv(csv_path, index=False)
    txt_path = os.path.join(OUTDIR, f'source_comparison_{stamp}.txt')
    with open(txt_path, 'w') as fh:
        fh.write('Source comparison: h5ad vs data_for_paper\n')
        fh.write(f'generated {stamp}\ndropped genes: {DROP_GENES}\n\n')
        fh.write(df[struct].to_string(index=False, float_format=fmt) + '\n\n')
        fh.write(df[spec].to_string(index=False, float_format=fmt) + '\n\n')
        fh.write(piv.to_string(float_format=fmt) + '\n\n')
        fh.write('Spearman corr with GMP-Cor:\n')
        fh.write(corr.to_string(float_format=fmt) + '\n')
    json_path = os.path.join(OUTDIR, f'source_comparison_{stamp}.json')
    with open(json_path, 'w') as fh:
        json.dump({
            'params': {
                'h5ad': H5AD, 'paper_dir': PAPER, 'dropped_genes': DROP_GENES,
                'batch_to_file': BATCH_TO_FILE, 'norm': NORM, 'log': LOG,
                'norm_method': NORM_METHOD, 'norm_sum': NORM_SUM, 'seed': SEED,
                'scramble_reps': 10,
                'gmp_cor_definition': 'sum(max(lambda_i - max_scrambled_lambda, 0))',
            },
            'results': records,
            'spearman_with_gmp_cor': corr.to_dict(),
        }, fh, indent=2)
    print(f'\nwrote:\n  {csv_path}\n  {txt_path}\n  {json_path}')


if __name__ == '__main__':
    main()
