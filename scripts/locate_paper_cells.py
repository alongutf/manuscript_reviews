"""Where do the published data_for_paper/sample_15a cells come from?

The notebook's reconstructed paper path reproduces only ~5% of that file, while
the umap path reproduces its file at 96%. This script profiles the published
paper barcodes inside the unfiltered matrix -- their depth, rRNA content and
rank -- to identify what criterion actually selected them.

Also caches per-barcode totals to an .npz so further questions are cheap.
"""
import os
import sys
import datetime

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
UNFILT = os.path.join(OTHER, 'data', 'sample_15a_unfiltered.csv')
BIOTYPE = os.path.join(OTHER, 'filtered_data', 'k12_biotype_map.csv')
PAPER_F = os.path.join(ROOT, 'data_for_paper', 'sample_15a_filtered.csv')
UMAP_F = os.path.join(ROOT, 'data_for_umap', 'sample_15a_filtered.csv')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE = os.path.join(OUTDIR, 'sample_15a_barcode_stats.npz')
CHUNK = 20000


def protein_coding_mask(genes):
    bt = pd.read_csv(BIOTYPE)
    pc = bt.gene[(bt.biotype != 'tRNA') & (bt.biotype != 'rRNA')].astype(str)
    pc = set(v.casefold() for v in pc)
    names = [str(v).casefold().replace('lelobekk_', '') for v in genes]
    return np.array([v in pc for v in names])


def build_cache():
    header = pd.read_csv(UNFILT, nrows=0)
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    pc = protein_coding_mask(genes)
    rr = np.array([str(g).startswith('16s') for g in genes])

    bcs, tot, mrna, rrna, det_pc = [], [], [], [], []
    for i, ch in enumerate(pd.read_csv(UNFILT, index_col=0, chunksize=CHUNK,
                                       dtype=dtypes)):
        V = ch.values
        bcs.append(np.array([str(b).split('-')[0] for b in ch.index]))
        tot.append(V.sum(1))
        mrna.append(V[:, pc].sum(1))
        rrna.append(V[:, rr].sum(1))
        det_pc.append((V[:, pc] > 0).sum(1))
        print(f'  chunk {i}')
    out = dict(
        barcode=np.concatenate(bcs), total=np.concatenate(tot),
        mrna=np.concatenate(mrna), rrna=np.concatenate(rrna),
        detected_pc=np.concatenate(det_pc),
        genes=genes, pc_mask=pc,
    )
    np.savez_compressed(CACHE, **out)
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    if os.path.exists(CACHE):
        print(f'loading cache {CACHE}')
        d = dict(np.load(CACHE, allow_pickle=True))
    else:
        print('building cache (full matrix pass) ...')
        d = build_cache()

    bc = d['barcode']
    pos = {b: i for i, b in enumerate(bc)}
    print(f'unfiltered barcodes: {len(bc)}')

    sets = {}
    for name, path in [('paper', PAPER_F), ('umap', UMAP_F)]:
        idx = pd.read_csv(path, index_col=0, usecols=[0]).index
        sets[name] = np.array([str(i).split('-')[0] for i in idx])

    # rank within the unfiltered pool (1 = highest)
    rank_tot = pd.Series(-d['total']).rank(method='first').values
    rank_mrna = pd.Series(-d['mrna']).rank(method='first').values

    rows = []
    for name, arr in sets.items():
        ix = np.array([pos[b] for b in arr if b in pos])
        rows.append({
            'file': name,
            'n_found_in_unfiltered': len(ix),
            'median_total': float(np.median(d['total'][ix])),
            'median_mrna': float(np.median(d['mrna'][ix])),
            'median_rrna_frac': float(np.median(d['rrna'][ix] / np.maximum(d['total'][ix], 1))),
            'median_detected_pc': float(np.median(d['detected_pc'][ix])),
            'median_rank_by_total': float(np.median(rank_tot[ix])),
            'median_rank_by_mrna': float(np.median(rank_mrna[ix])),
            'pct_in_top1000_by_total': float(100 * (rank_tot[ix] <= 1000).mean()),
            'pct_in_top1000_by_mrna': float(100 * (rank_mrna[ix] <= 1000).mean()),
        })
    df = pd.DataFrame(rows)
    print('\n=== published cell sets profiled in the unfiltered matrix ===')
    print(df.to_string(index=False, float_format=lambda v: '%.1f' % v))

    # what does the whole pool look like, for reference
    print('\n=== unfiltered pool reference ===')
    for k in ['total', 'mrna', 'detected_pc']:
        v = d[k]
        print(f'  {k:12s} median {np.median(v):8.1f}  p90 {np.percentile(v, 90):8.1f}  max {v.max():8.0f}')

    # is the paper set a contiguous rank band?
    ixp = np.array([pos[b] for b in sets['paper'] if b in pos])
    for lbl, r in [('by total', rank_tot[ixp]), ('by mRNA', rank_mrna[ixp])]:
        print(f'\npaper set rank distribution {lbl}: '
              f'min {r.min():.0f} p25 {np.percentile(r, 25):.0f} median {np.median(r):.0f} '
              f'p75 {np.percentile(r, 75):.0f} max {r.max():.0f}')

    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    df.to_csv(os.path.join(OUTDIR, f'locate_paper_cells_{stamp}.csv'), index=False)
    print(f'\nwrote results/cluster_gmp_cor/locate_paper_cells_{stamp}.csv')


if __name__ == '__main__':
    main()
