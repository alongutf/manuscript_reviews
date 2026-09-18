"""Does selecting on TOTAL counts systematically elevate reg's depth over dis?

Test: apply ONE identical rule to all five samples (top-N by total counts) and
compare the resulting depth. If total-based selection is what elevates reg, reg
should come out deeper than dis under this uniform rule. If reg is not elevated
here, the elevation in the published runs comes from the asymmetric umi_min
(200 reg / 400 dis) and the equate_dims draw, not from the criterion itself.

Depth is reported three ways, because they are not interchangeable:
  panel  = every column except 16s_mature / 16s_unprocessed / mCherry / kanR
           -- this is what GMP-Cor actually runs on
  pc     = protein_coding + pseudogene + ncRNA only (biotype map)
  tmRNA  = the single largest non-coding column, broken out because its share
           differs ~3x between conditions
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
CACHE_DIR = os.path.join(OUTDIR, 'panel_depth')

SAMPLES = {
    'exp':  'sample_2b_unfiltered.csv',
    'dis1': 'sample_13a_unfiltered.csv',
    'dis2': 'sample_15a_unfiltered.csv',
    'reg1': 'sample_13b_unfiltered.csv',
    'reg2': 'sample_15b_unfiltered.csv',
}
COND = {'exp': 'exp', 'dis1': 'dis', 'dis2': 'dis', 'reg1': 'reg', 'reg2': 'reg'}
CHUNK = 20000
# only these columns leave the panel -- everything else enters GMP-Cor
DROP_GENES = ['16s_mature', '16s_unprocessed', 'mCherry', 'kanR']
TOPN = [500, 1000, 2000]


def masks(genes):
    bt = pd.read_csv(BIOTYPE)
    m = dict(zip(bt.gene.astype(str).str.casefold(), bt.biotype))
    key = [str(g).casefold().replace('lelobekk_', '') for g in genes]
    cls = np.array([m.get(k, 'UNMAPPED') for k in key])
    panel = ~np.isin(genes, DROP_GENES)
    pc = np.isin(cls, ['protein_coding', 'pseudogene', 'ncRNA'])
    rr = np.array([str(g).casefold().startswith('16s') for g in genes])
    tm = np.array([str(g) == 'tmRNA' for g in genes])
    return panel, pc, rr, tm


def build(sample, fname):
    cache = os.path.join(CACHE_DIR, f'{sample}.npz')
    if os.path.exists(cache):
        return dict(np.load(cache, allow_pickle=True))
    path = os.path.join(DATA, fname)
    print(f'  {sample}: scanning ...')
    header = pd.read_csv(path, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    panel, pc, rr, tm = masks(genes)

    acc = {k: [] for k in ('total', 'panel', 'pc', 'rrna', 'tmrna',
                           'det_panel', 'det_pc')}
    for ch in pd.read_csv(path, index_col=0, chunksize=CHUNK, dtype=dtypes):
        V = ch.values
        acc['total'].append(V.sum(1))
        acc['panel'].append(V[:, panel].sum(1))
        acc['pc'].append(V[:, pc].sum(1))
        acc['rrna'].append(V[:, rr].sum(1))
        acc['tmrna'].append(V[:, tm].sum(1).ravel())
        acc['det_panel'].append((V[:, panel] > 0).sum(1))
        acc['det_pc'].append((V[:, pc] > 0).sum(1))
    out = {k: np.concatenate(v) for k, v in acc.items()}
    np.savez_compressed(cache, **out)
    return out


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    rows = []
    for s, f in SAMPLES.items():
        d = build(s, f)
        tot = d['total'].astype(float)
        order = np.argsort(-tot, kind='stable')
        for N in TOPN:
            i = order[:N]
            panel = d['panel'][i].astype(float)
            rows.append(dict(
                sample=s, cond=COND[s], n=N,
                mean_total=tot[i].mean(),
                panel_depth=panel.mean(),
                pc_depth=d['pc'][i].astype(float).mean(),
                tmrna=d['tmrna'][i].astype(float).mean(),
                tmrna_share_of_panel=d['tmrna'][i].sum() / panel.sum(),
                rrna_frac=d['rrna'][i].sum() / tot[i].sum(),
                det_panel=d['det_panel'][i].astype(float).mean(),
                det_pc=d['det_pc'][i].astype(float).mean(),
            ))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUTDIR, f'panel_depth_by_criterion_{stamp}.csv'),
              index=False)
    pd.set_option('display.width', 220, 'display.max_columns', 40)
    print(df.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
