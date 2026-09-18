"""Filtering on total counts vs on mRNA counts: do they select the same cells?

Section 8 of documents/gmp_cor_provenance_analysis.md notes that `umi_min`
thresholds TOTAL counts, which are 70-92% rRNA depending on the sample. This
script asks the practical question: if you swap the selection variable from
total to mRNA, how much does the selected barcode set change, and what happens
to the depth/detection of the cells you end up with?

Two matched comparisons, both at equal n so the sets are directly comparable:
  1. top-N by total   vs  top-N by mRNA          (N = 1000, the published n)
  2. threshold pool   vs  top-n_pool by mRNA     (the published umi_min pools)

Reads the per-barcode caches written by scripts/total_vs_mrna_stats.py.
"""
import os
import datetime

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE_DIR = os.path.join(OUTDIR, 'total_mrna_stats')

SAMPLES = ['exp', 'dis1', 'dis2', 'reg1', 'reg2']
TOPN = [1000, 2000]
# published-style total-count thresholds, per condition (see section 7 of the doc)
UMI_MIN = {'exp': 400, 'dis1': 400, 'dis2': 400, 'reg1': 200, 'reg2': 200}
UMI_MAX = 20000


def describe(d, idx):
    t = d['total'][idx].astype(float)
    m = d['mrna'][idx].astype(float)
    r = d['rrna'][idx].astype(float)
    det = d['detected_pc'][idx].astype(float)
    return dict(mean_total=t.mean(), mean_mrna=m.mean(), median_mrna=float(np.median(m)),
                mean_detected=det.mean(), rrna_frac=r.sum() / t.sum())


def compare(sample, d, sel_a, sel_b, label_a, label_b, scheme):
    a, b = set(sel_a.tolist()), set(sel_b.tolist())
    inter = len(a & b)
    row = dict(sample=sample, scheme=scheme, n_a=len(a), n_b=len(b),
               overlap=inter, overlap_frac=inter / max(len(a), 1),
               jaccard=inter / max(len(a | b), 1))
    for pre, sel in ((label_a, sel_a), (label_b, sel_b)):
        for k, v in describe(d, sel).items():
            row[f'{pre}_{k}'] = v
    return row


def main():
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    rows = []
    for s in SAMPLES:
        d = dict(np.load(os.path.join(CACHE_DIR, f'{s}.npz'), allow_pickle=True))
        tot, mrna = d['total'], d['mrna']
        # ties in the count vectors are broken deterministically by argsort order
        order_t = np.argsort(-tot, kind='stable')
        order_m = np.argsort(-mrna, kind='stable')
        for N in TOPN:
            rows.append(compare(s, d, order_t[:N], order_m[:N],
                                'tot', 'mrna', f'top{N}'))
        # published-style pool: threshold on total, vs the same number of cells
        # taken by mRNA rank instead
        pool = np.where((tot > UMI_MIN[s]) & (tot < UMI_MAX))[0]
        n = len(pool)
        rows.append(compare(s, d, pool, order_m[:n], 'tot', 'mrna',
                            f'pool_umi>{UMI_MIN[s]} (n={n})'))
    df = pd.DataFrame(rows)
    out = os.path.join(OUTDIR, f'filter_criterion_overlap_{stamp}.csv')
    df.to_csv(out, index=False)
    pd.set_option('display.width', 250, 'display.max_columns', 60)
    cols = ['sample', 'scheme', 'n_a', 'overlap', 'overlap_frac', 'jaccard',
            'tot_mean_mrna', 'mrna_mean_mrna', 'tot_mean_detected',
            'mrna_mean_detected', 'tot_rrna_frac', 'mrna_rrna_frac']
    print(df[cols].round(3).to_string(index=False))
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
