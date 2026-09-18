"""Does selecting barcodes on TOTAL counts give more or less mRNA, and does the
penalty differ by condition?

Selecting on mRNA is not a neutral fix -- it conditions on the quantity the
downstream analysis measures. So the question is narrower: relative to what is
available, does a total-count criterion enrich or deplete mRNA content, and is
that bias the same in every sample?

Three measurements, all at matched n so nothing is confounded by set size:
  A. mRNA recovery   = mean mRNA of the top-n by TOTAL
                     / mean mRNA of the top-n by mRNA        (<1 = a penalty)
  B. enrichment vs random = mean mRNA of the top-n by TOTAL
                     / mean mRNA of a random n from the eligible pool
  C. bias at fixed mRNA: within a narrow mRNA band, compare cells the total
     criterion accepts vs rejects -- if it is neutral they should look alike.

Reads the per-barcode caches from scripts/total_vs_mrna_stats.py.
"""
import os
import datetime

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
CACHE_DIR = os.path.join(OUTDIR, 'total_mrna_stats')

SAMPLES = ['exp', 'dis1', 'dis2', 'reg1', 'reg2']
COND = {'exp': 'exp', 'dis1': 'dis', 'dis2': 'dis', 'reg1': 'reg', 'reg2': 'reg'}
TOPN = [500, 1000, 2000, 5000]
# pool of non-empty droplets to draw the random reference from
POOL_MIN = 50
SEED = 0
N_BOOT = 200
# mRNA band used for the fixed-mRNA comparison, and the total cut tested there
MRNA_BAND = (40, 60)


def load(s):
    return dict(np.load(os.path.join(CACHE_DIR, f'{s}.npz'), allow_pickle=True))


def recovery_table():
    rng = np.random.default_rng(SEED)
    rows = []
    for s in SAMPLES:
        d = load(s)
        tot = d['total'].astype(float)
        mrna = d['mrna'].astype(float)
        det = d['detected_pc'].astype(float)
        ot = np.argsort(-tot, kind='stable')
        om = np.argsort(-mrna, kind='stable')
        pool = np.where(tot > POOL_MIN)[0]
        for N in TOPN:
            if N > len(pool):
                continue
            sel_t, sel_m = ot[:N], om[:N]
            # random reference: n draws from the same eligible pool
            boot = [mrna[rng.choice(pool, N, replace=False)].mean()
                    for _ in range(N_BOOT)]
            rand_mrna = float(np.mean(boot))
            rows.append(dict(
                sample=s, cond=COND[s], n=N,
                mrna_by_total=mrna[sel_t].mean(),
                mrna_by_mrna=mrna[sel_m].mean(),
                mrna_random=rand_mrna,
                recovery=mrna[sel_t].mean() / mrna[sel_m].mean(),
                enrichment_vs_random=mrna[sel_t].mean() / rand_mrna,
                det_by_total=det[sel_t].mean(),
                det_recovery=det[sel_t].mean() / det[sel_m].mean(),
            ))
    return pd.DataFrame(rows)


def fixed_mrna_table():
    """Within a narrow mRNA band every cell has ~the same mRNA content. Ask what
    a total-count cut does to that homogeneous group."""
    rows = []
    lo, hi = MRNA_BAND
    for s in SAMPLES:
        d = load(s)
        tot = d['total'].astype(float)
        mrna = d['mrna'].astype(float)
        det = d['detected_pc'].astype(float)
        band = (mrna >= lo) & (mrna <= hi)
        if band.sum() < 50:
            rows.append(dict(sample=s, cond=COND[s], n_band=int(band.sum())))
            continue
        # the total cut that keeps half the band -- the cut's own median
        cut = float(np.median(tot[band]))
        keep = band & (tot >= cut)
        drop = band & (tot < cut)
        rows.append(dict(
            sample=s, cond=COND[s], n_band=int(band.sum()), total_cut=cut,
            kept_mrna=mrna[keep].mean(), dropped_mrna=mrna[drop].mean(),
            kept_detected=det[keep].mean(), dropped_detected=det[drop].mean(),
            kept_rrna_frac=float(np.mean(d['rrna'][keep] / tot[keep])),
            dropped_rrna_frac=float(np.mean(d['rrna'][drop] / tot[drop])),
            # spread of total within a fixed mRNA content = how much noise the
            # total criterion injects
            total_iqr_ratio=float(np.percentile(tot[band], 75) /
                                  np.percentile(tot[band], 25)),
        ))
    return pd.DataFrame(rows)


def main():
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    pd.set_option('display.width', 220, 'display.max_columns', 40)
    rec = recovery_table()
    fix = fixed_mrna_table()
    rec.to_csv(os.path.join(OUTDIR, f'total_selection_bias_{stamp}.csv'), index=False)
    fix.to_csv(os.path.join(OUTDIR, f'total_selection_fixedmrna_{stamp}.csv'), index=False)
    print('A/B  mRNA recovery and enrichment, matched n')
    print(rec.round(3).to_string(index=False))
    print('\nby condition (mean recovery over n):')
    print(rec.groupby('cond')[['recovery', 'enrichment_vs_random',
                               'det_recovery']].mean().round(3).to_string())
    print(f'\nC  within mRNA band {MRNA_BAND}, split at the band median total')
    print(fix.round(3).to_string(index=False))
    print(f'\nwrote results to {OUTDIR}')


if __name__ == '__main__':
    main()
