"""Total probe counts vs mRNA content: how large is the gap, and do they track?

Section 8 of documents/gmp_cor_provenance_analysis.md shows that rRNA fraction
differs between conditions, so thresholding on TOTAL counts does not mean the
same thing across samples. This script quantifies that directly for all five
samples: the rRNA fraction, the total->mRNA ratio, and the correlation between
total and mRNA per barcode (Pearson on raw and on log1p, Spearman on rank).

Per-barcode caches (total / mRNA / rRNA / detected) are written to
results/cluster_gmp_cor/total_mrna_stats/<sample>.npz so the 840 MB CSVs are
scanned only once. sample_15a (dis2) reuses the existing sample_15a cache.
"""
import os
import datetime

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
DATA = os.path.join(OTHER, 'data')
BIOTYPE = os.path.join(OTHER, 'filtered_data', 'k12_biotype_map.csv')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE_DIR = os.path.join(OUTDIR, 'total_mrna_stats')
LEGACY_15A = os.path.join(OUTDIR, 'sample_15a_barcode_stats.npz')

SAMPLES = {
    'exp':  'sample_2b_unfiltered.csv',
    'dis1': 'sample_13a_unfiltered.csv',
    'dis2': 'sample_15a_unfiltered.csv',
    'reg1': 'sample_13b_unfiltered.csv',
    'reg2': 'sample_15b_unfiltered.csv',
}
CHUNK = 20000
# barcode subsets to report separately: (label, predicate on total counts)
POOLS = [('all', 0), ('total>50', 50), ('total>200', 200), ('total>400', 400)]


def protein_coding_mask(genes):
    bt = pd.read_csv(BIOTYPE)
    pc = bt.gene[(bt.biotype != 'tRNA') & (bt.biotype != 'rRNA')].astype(str)
    pc = set(v.casefold() for v in pc)
    names = [str(v).casefold().replace('lelobekk_', '') for v in genes]
    return np.array([v in pc for v in names])


def build_stats(sample, fname):
    cache = os.path.join(CACHE_DIR, f'{sample}.npz')
    if os.path.exists(cache):
        print(f'  {sample}: cache hit')
        return dict(np.load(cache, allow_pickle=True))
    if sample == 'dis2' and os.path.exists(LEGACY_15A):
        print(f'  {sample}: reusing {os.path.basename(LEGACY_15A)}')
        d = dict(np.load(LEGACY_15A, allow_pickle=True))
        out = dict(barcode=d['barcode'], total=d['total'], mrna=d['mrna'],
                   rrna=d['rrna'], detected_pc=d['detected_pc'])
        np.savez_compressed(cache, **out)
        return out

    path = os.path.join(DATA, fname)
    print(f'  {sample}: scanning {fname} ...')
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    pc = protein_coding_mask(genes)
    rr = np.array([str(g).casefold().startswith('16s') for g in genes])

    bcs, tot, mrna, rrna, det = [], [], [], [], []
    for i, ch in enumerate(pd.read_csv(path, index_col=0, chunksize=CHUNK,
                                       dtype=dtypes)):
        V = ch.values
        bcs.append(np.array([str(b).split('-')[0] for b in ch.index]))
        tot.append(V.sum(1))
        mrna.append(V[:, pc].sum(1))
        rrna.append(V[:, rr].sum(1))
        det.append((V[:, pc] > 0).sum(1))
    out = dict(barcode=np.concatenate(bcs), total=np.concatenate(tot),
               mrna=np.concatenate(mrna), rrna=np.concatenate(rrna),
               detected_pc=np.concatenate(det))
    np.savez_compressed(cache, **out)
    return out


def summarize(sample, d):
    rows = []
    tot, mrna, rrna = (d['total'].astype(float), d['mrna'].astype(float),
                       d['rrna'].astype(float))
    for label, thr in POOLS:
        sel = tot > thr
        n = int(sel.sum())
        if n < 10:
            continue
        t, m, r = tot[sel], mrna[sel], rrna[sel]
        # aggregate (pooled) fractions, then per-barcode means
        row = dict(
            sample=sample, pool=label, n=n,
            mean_total=t.mean(), mean_mrna=m.mean(), mean_rrna=r.mean(),
            median_total=float(np.median(t)), median_mrna=float(np.median(m)),
            rrna_frac_pooled=r.sum() / t.sum(),
            mrna_frac_pooled=m.sum() / t.sum(),
            rrna_frac_percell_mean=float(np.mean(r / np.maximum(t, 1))),
            ratio_total_over_mrna=t.sum() / max(m.sum(), 1),
            pearson_raw=stats.pearsonr(t, m)[0],
            pearson_log=stats.pearsonr(np.log1p(t), np.log1p(m))[0],
            spearman=stats.spearmanr(t, m)[0],
            # how much of the total's variance the mRNA part explains
            r2_log=stats.pearsonr(np.log1p(t), np.log1p(m))[0] ** 2,
            cv_mrna_frac=float(np.std(m / np.maximum(t, 1)) /
                               max(np.mean(m / np.maximum(t, 1)), 1e-12)),
        )
        rows.append(row)
    return rows


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    rows = []
    for s, f in SAMPLES.items():
        d = build_stats(s, f)
        rows += summarize(s, d)
    df = pd.DataFrame(rows)
    out = os.path.join(OUTDIR, f'total_vs_mrna_{stamp}.csv')
    df.to_csv(out, index=False)
    pd.set_option('display.width', 200, 'display.max_columns', 50)
    print(df.round(4).to_string(index=False))
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
