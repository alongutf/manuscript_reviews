"""Per-gene count composition per sample, to see what the analysis panel holds.

The panel used by equate_reg_n1000 drops only ['16s_mature','16s_unprocessed',
'LELOBEKK','kanR','mCherry']. But 574 of the 4184 matrix columns are not in the
K-12 biotype map at all (plasmid/marker/probe columns), and the fraction of
counts they carry differs sharply between conditions (~3-4% in exp/dis vs ~15%
in reg). Anything left in the panel contributes to the reported "depth" and to
the correlation spectrum, so this script reports, per sample, the count share of
each biotype class and the top non-protein-coding contributors.

Computed over all barcodes and over the cell-like pool (total > 200).
"""
import os
import datetime

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
DATA = os.path.join(OTHER, 'data')
BIOTYPE = os.path.join(OTHER, 'filtered_data', 'k12_biotype_map.csv')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE_DIR = os.path.join(OUTDIR, 'gene_sums')

SAMPLES = {
    'exp':  'sample_2b_unfiltered.csv',
    'dis1': 'sample_13a_unfiltered.csv',
    'dis2': 'sample_15a_unfiltered.csv',
    'reg1': 'sample_13b_unfiltered.csv',
    'reg2': 'sample_15b_unfiltered.csv',
}
CHUNK = 20000
POOL_MIN = 200          # cell-like pool, matching the reg umi_min actually used
DROP_GENES = ['16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR', 'mCherry']


def classify(genes):
    bt = pd.read_csv(BIOTYPE)
    m = dict(zip(bt.gene.astype(str).str.casefold(), bt.biotype))
    key = [str(g).casefold().replace('lelobekk_', '') for g in genes]
    return np.array([m.get(k, 'UNMAPPED') for k in key])


def gene_sums(sample, fname):
    cache = os.path.join(CACHE_DIR, f'{sample}.npz')
    if os.path.exists(cache):
        return dict(np.load(cache, allow_pickle=True))
    path = os.path.join(DATA, fname)
    print(f'  {sample}: scanning ...')
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    s_all = np.zeros(len(genes), dtype=np.int64)
    s_pool = np.zeros(len(genes), dtype=np.int64)
    n_pool = 0
    for ch in pd.read_csv(path, index_col=0, chunksize=CHUNK, dtype=dtypes):
        V = ch.values
        s_all += V.sum(0, dtype=np.int64)
        sel = V.sum(1) > POOL_MIN
        if sel.any():
            s_pool += V[sel].sum(0, dtype=np.int64)
            n_pool += int(sel.sum())
    out = dict(genes=genes, sum_all=s_all, sum_pool=s_pool,
               n_pool=np.array(n_pool))
    np.savez_compressed(cache, **out)
    return out


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    shares, tops = [], []
    for s, f in SAMPLES.items():
        d = gene_sums(s, f)
        genes = np.array([str(g) for g in d['genes']])
        cls = classify(genes)
        drop = np.isin(genes, DROP_GENES)
        for scope in ('sum_all', 'sum_pool'):
            v = d[scope].astype(float)
            tot = v.sum()
            row = dict(sample=s, scope=scope, total=tot, n_pool=int(d['n_pool']))
            for c in ['protein_coding', 'pseudogene', 'ncRNA', 'tRNA', 'rRNA',
                      'UNMAPPED']:
                row[c] = v[cls == c].sum() / tot
            # what the panel keeps: everything except DROP_GENES
            row['panel_share'] = v[~drop].sum() / tot
            # of what the panel keeps, how much is NOT protein_coding
            keep = ~drop
            row['panel_nonPC_share'] = (v[keep & (cls != 'protein_coding')].sum()
                                        / v[keep].sum())
            shares.append(row)
        # top non-protein-coding contributors that survive the drop list
        v = d['sum_pool'].astype(float)
        keep = ~drop & (cls != 'protein_coding')
        idx = np.argsort(-v * keep)[:8]
        for i in idx:
            tops.append(dict(sample=s, gene=genes[i], biotype=cls[i],
                             share_of_pool=v[i] / v[~drop].sum()))
    sh = pd.DataFrame(shares)
    tp = pd.DataFrame(tops)
    sh.to_csv(os.path.join(OUTDIR, f'gene_composition_{stamp}.csv'), index=False)
    tp.to_csv(os.path.join(OUTDIR, f'gene_composition_top_{stamp}.csv'), index=False)
    pd.set_option('display.width', 220, 'display.max_columns', 40)
    print(sh.round(4).to_string(index=False))
    print('\ntop non-protein-coding genes surviving DROP_GENES (pool total>200),')
    print('as a share of the panel counts:')
    print(tp.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
