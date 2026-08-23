"""High-B permutation run for the samples near the decision boundary.

Motivation
----------
The parametric tail fit (null_distribution_fit.py) is only identifiable where the
null sample constrains it. For samples with a huge z_1 the fitted families
disagree by hundreds to tens of thousands of orders of magnitude, so no number is
defensible -- but the conclusion is not in doubt either. The samples that
actually need a resolved p-value are the ones near the boundary (small |z_1|),
and for those the answer does not require extrapolation at all: raising B until
the EMPIRICAL permutation p-value resolves gives an assumption-free answer.

This script therefore reruns only the low-|z_1| samples at high B.

Usage:
    python scripts/eigenvalue_permutation_highB.py [B] [sample ...]
"""
import os
import sys
import json
import time
import platform
from datetime import datetime

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import null_distribution_fit as ndf
from eigenvalue_permutation_test import (
    load_matrix, prep, spec, scramble, perm_pvalue, validate,
    DATA_DIR, EV_DIR, OUT, SEED, K, TRACKER_GENES,
    NORM_METHOD, NORM_SUM, LOG, MIN_CELLS, MIN_GENES,
)

# samples within a few sigma of the null: these are the ones where a resolved
# p-value changes the reading, and where high B actually delivers one
DEFAULT_SAMPLES = ['SHX_biorep2B', 'sample_13a_filtered', 'sample_15a_filtered']
DEFAULT_B = 2000


def run(B=DEFAULT_B, samples=None):
    samples = samples or DEFAULT_SAMPLES
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    tag = 'eigenvalue_permutation_highB_' + stamp
    titles = pd.read_excel(os.path.join(EV_DIR, 'titles.xlsx')).set_index('file_name')
    rng = np.random.default_rng(SEED + 1)

    rows, per_sample, draws_store = [], {}, {}
    t_start = time.time()

    for name in samples:
        t0 = time.time()
        stored = np.load(os.path.join(EV_DIR, name + '.npy'))[0]
        m_raw, dropped = load_matrix(os.path.join(DATA_DIR, name + '.csv'))
        M = prep(m_raw)
        n, p = M.shape
        obs = spec(M)

        val = validate(m_raw, M, stored)
        assert val['max_abs_diff_repo_vs_stored'] == 0.0, name + ': repo pipeline mismatch'

        null = np.empty((B, K))
        for b in range(B):
            null[b] = spec(scramble(M, rng))[:K]
            if (b + 1) % 250 == 0:
                print('  {} {}/{}'.format(name, b + 1, B), flush=True)
        draws_store[name] = null

        mu, sd = null.mean(0), null.std(0)
        z = (obs[:K] - mu) / sd
        p_emp = np.array([perm_pvalue(null[:, k], obs[k]) for k in range(K)])
        n_exceed = int((null[:, 0] >= obs[0]).sum())

        fit = ndf.parametric_pvalue(null[:, 0], obs[0], B=B)
        gpd = ndf.gpd_tail_p(null[:, 0], obs[0])
        calib = ndf.calibration_check(null[:, 0], seed=SEED)

        # Wilson interval on the empirical p, so the Monte Carlo uncertainty of the
        # count-based answer is explicit
        from scipy.stats import beta as _beta
        lo = _beta.ppf(0.025, n_exceed + 1, B - n_exceed) if n_exceed >= 0 else 0.0
        hi = _beta.ppf(0.975, n_exceed + 1, B - n_exceed)

        title = titles.title.get(name + '.npy', name)
        cat = titles.category.get(name + '.npy', '?')
        el = time.time() - t0
        rows.append(dict(sample=name, title=title, cat=cat, n=n, p=p,
                         lambda_1=obs[0], null_mean=mu[0], null_sd=sd[0], z_1=z[0],
                         n_exceed=n_exceed, B=B, p_empirical=p_emp[0],
                         p_emp_ci_lo=float(lo), p_emp_ci_hi=float(hi),
                         p_parametric=fit['p_parametric'],
                         log10_p_parametric=fit['log10_p_parametric'],
                         log10_p_lo=fit['log10_p_min'], log10_p_hi=fit['log10_p_max'],
                         best_family=fit['best_family'],
                         gpd_xi=gpd.get('xi'), gpd_log10_p=gpd.get('log10_p'),
                         gpd_beyond_endpoint=gpd.get('obs_beyond_endpoint'),
                         pit_ks_p=calib.get('pit_ks_p')))
        per_sample[name] = dict(title=title, category=cat, n_cells=int(n), n_genes=int(p),
                                dropped_tracker_genes=dropped, validation=val,
                                lambda_obs_topK=obs[:K].tolist(),
                                null_mean_topK=mu.tolist(), null_sd_topK=sd.tolist(),
                                z_topK=z.tolist(), p_empirical_topK=p_emp.tolist(),
                                n_exceed_rank1=n_exceed,
                                parametric_fit_rank1=fit, gpd_rank1=gpd,
                                calibration_rank1=calib, runtime_sec=round(el, 1))
        print('{:24s} z1={:6.2f}  exceedances={}/{}  p_emp={:.5f} '
              '[{:.5f},{:.5f}]  ({:.0f}s)'.format(name, z[0], n_exceed, B,
                                                  p_emp[0], lo, hi, el), flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'raw', tag + '.csv'), index=False)
    np.savez_compressed(os.path.join(OUT, 'raw', tag + '_null_draws.npz'), **draws_store)

    with open(os.path.join(OUT, 'raw', tag + '.txt'), 'w') as f:
        f.write('High-B permutation run for boundary samples\n' + '=' * 78 + '\n')
        f.write('run: {}   B={}   K={}   seed={}\n\n'.format(stamp, B, K, SEED + 1))
        f.write('p_empirical is the exceedance count, (1+#{null>=obs})/(B+1); it needs no\n'
                'distributional assumption. p_emp_ci_lo/hi is its exact (Clopper-Pearson)\n'
                '95% Monte Carlo interval. Compare against p_parametric and gpd_log10_p:\n'
                'agreement means the fitted tail is trustworthy in this range.\n\n')
        f.write(df.to_string(index=False, float_format=lambda x: '{:.5f}'.format(x)))
        f.write('\n')

    log = dict(experiment=tag, timestamp=stamp,
               description='High-B permutation rerun for samples near the decision '
                           'boundary, where the empirical p resolves without any tail '
                           'extrapolation.',
               parameters=dict(B=B, K=K, seed=SEED + 1, samples=samples,
                               tracker_genes_removed=TRACKER_GENES,
                               normalization=dict(method=NORM_METHOD, target_sum=NORM_SUM,
                                                  log=LOG),
                               filtering=dict(min_cells_per_gene=MIN_CELLS,
                                              min_genes_per_cell=MIN_GENES),
                               min_attainable_empirical_p=1 / (B + 1)),
               environment=dict(python=platform.python_version(), numpy=np.__version__,
                                pandas=pd.__version__, platform=platform.platform()),
               total_runtime_sec=round(time.time() - t_start, 1),
               results=per_sample)
    with open(os.path.join(OUT, 'logs', tag + '.json'), 'w') as f:
        json.dump(log, f, indent=2)

    print('\nwrote results/permutation_test/{{raw,logs}}/' + tag + '.*')
    return df


if __name__ == '__main__':
    b = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_B
    smp = sys.argv[2:] or None
    run(b, smp)
