"""Add permutation-test columns to results/data_metrics/data_metrics.csv.

Columns added
-------------
permutation_p
    p_empirical from the B=2000 permutation test: the exceedance count
    (1 + #{lambda_1^perm >= lambda_1^obs}) / (B + 1), matched per dataset.
    No distributional assumption; censored at 1/(B+1) = 4.9975e-4.

gmp_cor_ci
    Uncertainty on GMP-Cor, propagated from the noise in the scrambled
    threshold:
        sigma = SD of lambda_max^scr across the B permutations
        N     = number of observed eigenvalues above the MEAN scrambled
                threshold, mean(lambda_max^scr)
        gmp_cor_ci = sqrt(N) * sigma
    GMP-Cor sums N terms each carrying an error sigma from the shared
    threshold, so the sum's SD is sqrt(N)*sigma.

Supporting columns (the inputs, so the CI can be checked by hand):
    sigma_lambda_max_scr    sigma above
    mean_lambda_max_scr     the threshold used to count N
    n_above_threshold       N

N is counted against the FULL observed spectrum from ev_data/<sample>.npy,
not a truncated top-K, so it is exact even where N > 50.

Usage:
    python scripts/add_permutation_metrics.py [perm_tag]
"""
import os
import sys

import numpy as np
import pandas as pd

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_TAG = 'eigenvalue_permutation_B2000_20260823_132215'
METRICS = os.path.join(_REPO, 'results', 'data_metrics', 'data_metrics.csv')
EV_DIR = os.path.join(_REPO, 'ev_data')


def main(tag=DEFAULT_TAG):
    """Join a permutation-test run onto data_metrics.csv and rewrite it in place.

    tag : basename (no extension) of the permutation run CSV under
        results/permutation_test/raw/, keyed by dataset in a 'sample' column
        with 'null_sd', 'null_mean' and 'p_empirical_1' fields (produced by the
        B=2000 permutation scripts, e.g. eigenvalue_permutation_full_B2000.py).

    Returns the updated data_metrics DataFrame (also written to disk).
    """
    perm_csv = os.path.join(_REPO, 'results', 'permutation_test', 'raw', tag + '.csv')
    if not os.path.exists(perm_csv):
        raise SystemExit('no permutation run at ' + perm_csv)

    try:                      # fail before touching anything if Excel holds the file
        with open(METRICS, 'r+b'):
            pass
    except OSError:
        raise SystemExit('cannot write {} (open in another program?). Close it and rerun.'
                         .format(METRICS))

    perm = pd.read_csv(perm_csv).set_index('sample')
    dm = pd.read_csv(METRICS, index_col=0)
    # data_metrics rows are keyed by CSV filename; permutation results are keyed by
    # the bare sample name, so strip the extension to line the two tables up.
    key = dm['file_name'].str.replace('.csv', '', regex=False)

    # bail out rather than silently leaving NaNs if a dataset in data_metrics has
    # no matching permutation-test row (e.g. the run hasn't finished for it yet)
    missing = sorted(set(key) - set(perm.index))
    if missing:
        raise SystemExit('no permutation result for: ' + ', '.join(missing))

    sigma, thr, n_above, p_emp = [], [], [], []
    for k in key:
        r = perm.loc[k]
        s = float(r['null_sd'])         # SD of lambda_max^scr over the B permutations
        t = float(r['null_mean'])       # mean scrambled threshold
        ev = np.load(os.path.join(EV_DIR, k + '.npy'))[0]   # full observed spectrum
        N = int((ev > t).sum())
        sigma.append(s)
        thr.append(t)
        n_above.append(N)
        p_emp.append(float(r['p_empirical_1']))

    dm['permutation_p'] = p_emp
    dm['sigma_lambda_max_scr'] = sigma
    dm['mean_lambda_max_scr'] = thr
    dm['n_above_threshold'] = n_above
    dm['gmp_cor_ci'] = np.sqrt(np.asarray(n_above, dtype=float)) * np.asarray(sigma)

    dm.to_csv(METRICS)
    print(dm[['file_name', 'sum_denoised_ev', 'permutation_p', 'sigma_lambda_max_scr',
              'n_above_threshold', 'gmp_cor_ci']].to_string(
        index=False, float_format=lambda x: '{:.6g}'.format(x)))
    print('\nupdated ' + METRICS)
    return dm


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_TAG)
