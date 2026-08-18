"""Reproduce the two notebook filter paths on sample_15a and see which one
generated each published file.

paper path (analysis_notebook cells 4->5->6):
    read unfiltered -> drop non-protein-coding (rRNA/tRNA)
                    -> filter_by_umi_count(100, 2000)      [counted on mRNA]
                    -> filter_by_gene_dispersion(min_dispersion=1)

umap path (analysis_notebook cell 9):
    read unfiltered -> filter_by_umi_count(200, 20000)     [counted on TOTAL, rRNA in]
                    -> filter_by_gene_dispersion()

Each candidate is compared by barcode overlap against the two published files,
data_for_paper/sample_15a_filtered.csv and data_for_umap/sample_15a_filtered.csv.
"""
import os
import sys
import json
import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
UNFILT = os.path.join(OTHER, 'data', 'sample_15a_unfiltered.csv')
BIOTYPE = os.path.join(OTHER, 'filtered_data', 'k12_biotype_map.csv')
PAPER_F = os.path.join(ROOT, 'data_for_paper', 'sample_15a_filtered.csv')
UMAP_F = os.path.join(ROOT, 'data_for_umap', 'sample_15a_filtered.csv')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')

CHUNK = 20000


def load_unfiltered(path):
    print(f'loading {path} ...')
    header = pd.read_csv(path, nrows=0)
    # first column holds barcodes; force str there and int32 everywhere else
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    parts, index = [], []
    for i, ch in enumerate(pd.read_csv(path, index_col=0, chunksize=CHUNK,
                                       dtype=dtypes)):
        parts.append(ch.values)
        index.append(ch.index.values)
        print(f'  chunk {i}: {ch.shape}')
    M = np.vstack(parts)
    idx = np.concatenate(index)
    cols = ch.columns.values
    print(f'loaded {M.shape}')
    return M, np.asarray(idx), np.asarray(cols)


def protein_coding_mask(genes):
    """Replicates analysis_notebook cell 5: drop tRNA and rRNA biotypes."""
    bt = pd.read_csv(BIOTYPE)
    pc = bt.gene[(bt.biotype != 'tRNA') & (bt.biotype != 'rRNA')].astype(str)
    pc = set(v.casefold() for v in pc)
    names = [str(v).casefold().replace('lelobekk_', '') for v in genes]
    return np.array([v in pc for v in names])


def umi_filter(tot, umi_min, umi_max, target_cells=None):
    """Replicates AnnMat.filter_by_umi_count, including the target_cells branch."""
    if target_cells is not None:
        srt = np.flip(np.sort(tot))
        n_above = int((tot > umi_max).sum())
        k = min(len(srt) - 1, n_above + target_cells)
        umi_min = srt[k]
    return (tot > umi_min) & (tot < umi_max)


def dispersion_mask(M, min_dispersion=1.0):
    """Replicates AnnMat.filter_by_gene_dispersion."""
    mean = M.mean(axis=0)
    var = M.var(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        disp = np.where(mean > 0, var / mean, 0.0)
    disp = np.nan_to_num(disp)
    return disp > min_dispersion


def summarize(M, bcs, genes, label, published):
    tot = M.sum(1)
    det = (M > 0).sum(1)
    rec = {
        'candidate': label,
        'n_cells': int(M.shape[0]),
        'n_genes': int(M.shape[1]),
        'mean_total': float(tot.mean()) if len(tot) else np.nan,
        'median_total': float(np.median(tot)) if len(tot) else np.nan,
        'min_total': float(tot.min()) if len(tot) else np.nan,
        'mean_detected': float(det.mean()) if len(det) else np.nan,
    }
    s = set(bcs)
    for pname, pbc in published.items():
        rec[f'overlap_{pname}'] = len(s & pbc)
        rec[f'pct_of_{pname}'] = round(100 * len(s & pbc) / len(pbc), 1)
    return rec


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    published = {}
    for name, path in [('paper', PAPER_F), ('umap', UMAP_F)]:
        d = pd.read_csv(path, index_col=0, usecols=[0])
        published[name] = set(str(i).split('-')[0] for i in d.index)
        print(f'published {name}: {len(published[name])} barcodes')

    M, bcs, genes = load_unfiltered(UNFILT)
    bcs = np.array([str(b).split('-')[0] for b in bcs])
    rrna_cols = [i for i, g in enumerate(genes) if str(g).startswith('16s')]
    print(f'rRNA columns: {[genes[i] for i in rrna_cols]}')

    pc = protein_coding_mask(genes)
    print(f'protein-coding-ish genes kept: {pc.sum()} / {len(genes)}')

    tot_all = M.sum(1)
    tot_mrna = M[:, pc].sum(1)
    print(f'total counts/barcode: median {np.median(tot_all):.0f}, max {tot_all.max()}')
    print(f'mRNA  counts/barcode: median {np.median(tot_mrna):.0f}, max {tot_mrna.max()}')

    records = []

    # ---- paper path: gene filter first, UMI counted on mRNA ----
    for (lo, hi, tgt) in [(100, 2000, None), (100, 2000, 1000), (400, 20000, None)]:
        cells = umi_filter(tot_mrna, lo, hi, tgt)
        sub = M[np.ix_(cells, pc)]
        keep_g = dispersion_mask(sub)
        final = sub[:, keep_g]
        lbl = f'paper-path umi({lo},{hi},target={tgt}) on mRNA'
        records.append(summarize(final, bcs[cells], None, lbl, published))
        print(f'{lbl}: {final.shape}')

    # ---- umap path: UMI counted on TOTAL incl. rRNA, no biotype filter ----
    for (lo, hi, tgt) in [(200, 20000, None), (200, 20000, 1000), (400, 20000, 1000)]:
        cells = umi_filter(tot_all, lo, hi, tgt)
        sub = M[cells]
        keep_g = dispersion_mask(sub)
        final = sub[:, keep_g]
        lbl = f'umap-path  umi({lo},{hi},target={tgt}) on TOTAL'
        records.append(summarize(final, bcs[cells], None, lbl, published))
        print(f'{lbl}: {final.shape}')

    df = pd.DataFrame(records)
    cols = ['candidate', 'n_cells', 'n_genes', 'mean_total', 'median_total',
            'min_total', 'mean_detected', 'overlap_paper', 'pct_of_paper',
            'overlap_umap', 'pct_of_umap']
    print('\n=== candidate reconstructions vs published files ===')
    print(df[cols].to_string(index=False, float_format=lambda v: '%.1f' % v))

    csv_path = os.path.join(OUTDIR, f'reproduce_filter_paths_{stamp}.csv')
    df.to_csv(csv_path, index=False)
    txt_path = os.path.join(OUTDIR, f'reproduce_filter_paths_{stamp}.txt')
    with open(txt_path, 'w') as fh:
        fh.write(f'Reproduction of notebook filter paths on sample_15a\n')
        fh.write(f'generated {stamp}\nunfiltered: {UNFILT}\n\n')
        fh.write(df[cols].to_string(index=False, float_format=lambda v: '%.1f' % v) + '\n')
    print(f'\nwrote:\n  {csv_path}\n  {txt_path}')


if __name__ == '__main__':
    main()
