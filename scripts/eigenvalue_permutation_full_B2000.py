"""Full B=2000 permutation run over all 15 datasets, retaining every null draw.

Reports the EMPIRICAL permutation p-value -- the exceedance count -- which needs
no distributional assumption:

    p_k = (1 + #{lambda_k^perm >= lambda_obs_k}) / (B + 1)

and saves every scrambled eigenvalue so that any distribution can be fitted
later without rerunning the permutations.

WHAT IS SAVED
-------------
raw/<tag>_null_draws.npz
    one (B, K_STORE) float64 array per sample: the scrambled eigenvalues for
    ranks 1..K_STORE of every replicate. Load with
        draws = dict(np.load(path));  draws['SHX_biorep2B'][:, 0]   # lambda_1
raw/<tag>_rank1.csv
    tidy long form of rank 1 -- the 15 x B datapoints -- one row per
    (sample, replicate): sample, title, category, replicate, lambda_1_scrambled,
    plus the observed lambda_1 for reference. Directly fittable.
raw/<tag>_observed.csv
    the observed eigenvalues, ranks 1..K_STORE, per sample.
raw/<tag>_null_moments.csv
    per-rank mean/sd/min/max of the null over the FULL spectrum (all n ranks),
    so bulk-of-spectrum work is still possible even though only the top
    K_STORE ranks are stored draw-by-draw.

Storing the full spectrum for every replicate would be ~15 x 2000 x 1000 floats
(~240 MB); the top K_STORE ranks cover the region where the signal lives, and
the moments file preserves the rest in summary form.

Runs samples in parallel across processes; each worker is pinned to a small
number of BLAS threads to avoid oversubscription.

Usage:
    python scripts/eigenvalue_permutation_full_B2000.py [B] [n_workers]
"""
import os
import sys

# must be set before numpy/BLAS is imported in the workers
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, os.environ.get('EV_BLAS_THREADS', '4'))

import json
import time
import platform
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _REPO)

B_DEFAULT = 2000
K_STORE = 50          # ranks kept draw-by-draw (~2.5% of p; signal lives in the top ~1%)
K_REPORT = 20         # ranks carried in the summary tables
SEED = 20260823
ALPHA = 0.05
OUT = os.path.join(_REPO, 'results', 'permutation_test')


def _work(args):
    """Run B permutations for one sample. Executed in a worker process.

    Returns a dict with the observed spectrum (top k_store ranks), the full (B,
    k_store) draw matrix, and running first/second moments plus min/max/exceedance
    counts over the FULL spectrum (all n ranks), computed online so the full draw
    matrix itself never needs to be kept in memory.
    """
    import numpy as np
    from eigenvalue_permutation_test import (
        load_matrix, prep, spec, scramble, validate, DATA_DIR, EV_DIR)

    name, B, seed = args
    t0 = time.time()
    stored = np.load(os.path.join(EV_DIR, name + '.npy'))[0]
    m_raw, dropped = load_matrix(os.path.join(DATA_DIR, name + '.csv'))
    M = prep(m_raw)
    n, p = M.shape
    obs = spec(M)

    val = validate(m_raw, M, stored)
    # Single-threaded, these agree with ev_data to exactly 0.0. Workers pin BLAS to a
    # different thread count, which changes LAPACK's summation order and perturbs the
    # eigenvalues by ~1e-14 absolute (~5e-16 relative). Assert a tight tolerance rather
    # than bitwise equality; the realised value is kept in `validation` for every sample.
    assert val['max_abs_diff_repo_vs_stored'] < 1e-9, name + ': repo pipeline mismatch'
    assert val['max_abs_diff_this_vs_stored'] < 1e-9, name + ': preprocessing mismatch'

    rng = np.random.default_rng(seed)
    k_store = min(K_STORE, len(obs))
    draws = np.empty((B, k_store))
    # running moments over the FULL spectrum, so the bulk is not lost
    s1 = np.zeros(len(obs))
    s2 = np.zeros(len(obs))
    fmin = np.full(len(obs), np.inf)
    fmax = np.full(len(obs), -np.inf)
    n_above_full = np.zeros(len(obs), dtype=np.int64)

    for b in range(B):
        e = spec(scramble(M, rng))
        draws[b] = e[:k_store]
        s1 += e                                    # running sum -> mean
        s2 += e * e                                # running sum of squares -> variance
        np.minimum(fmin, e, out=fmin)
        np.maximum(fmax, e, out=fmax)
        n_above_full += (e >= obs[:len(e)])        # exceedance count, every rank

    mean_full = s1 / B
    # E[X^2] - E[X]^2, Bessel-corrected (x B/(B-1)) to match np.std(ddof=1) used
    # elsewhere; clipped at 0 because floating-point error can drive it slightly
    # negative when the true variance is ~0
    var_full = np.maximum(s2 / B - mean_full ** 2, 0) * B / (B - 1)
    sd_full = np.sqrt(var_full)

    return dict(
        name=name, n=int(n), p=int(p), dropped=dropped, validation=val,
        obs=obs[:k_store].tolist(), draws=draws,
        mean_full=mean_full, sd_full=sd_full, min_full=fmin, max_full=fmax,
        n_above_full=n_above_full, n_ranks_full=int(len(obs)),
        runtime_sec=round(time.time() - t0, 1),
    )


def main(B=B_DEFAULT, n_workers=None):
    """Run the full B-replicate permutation test over every sample in data_for_paper/,
    in parallel across worker processes, and write all raw/log outputs."""
    from scipy.stats import beta as _beta, mannwhitneyu

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    tag = 'eigenvalue_permutation_B{}_{}'.format(B, stamp)
    titles = pd.read_excel(os.path.join(OUT, '..', '..', 'ev_data',
                                        'titles.xlsx')).set_index('file_name')
    data_dir = os.path.join(_REPO, 'data_for_paper')
    samples = sorted(f[:-4] for f in os.listdir(data_dir) if f.endswith('.csv'))
    n_workers = n_workers or min(len(samples), max(1, (os.cpu_count() or 4) // 2))

    print('B={}  samples={}  workers={}'.format(B, len(samples), n_workers), flush=True)
    t_start = time.time()

    # distinct, reproducible stream per sample: each worker seeds its own
    # np.random.default_rng, so samples running concurrently never share a stream
    seeds = {nm: SEED + 1000 + i for i, nm in enumerate(samples)}
    # cost scales with the matrix size, and the largest sample dominates wall time,
    # so submit longest-first rather than letting it start last
    order = sorted(samples, key=lambda nm: -os.path.getsize(
        os.path.join(data_dir, nm + '.csv')))
    results = {}
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(_work, (nm, B, seeds[nm])): nm for nm in order}
        for fut in as_completed(futs):
            r = fut.result()
            results[r['name']] = r
            print('  done {:30s} n={:5d}  ({:.0f}s)'.format(
                r['name'], r['n'], r['runtime_sec']), flush=True)

    # ------------------------------------------------------------ assemble
    rows, per_sample, draws_store = [], {}, {}
    rank1_long, obs_long, moments_long = [], [], []

    for name in samples:
        r = results[name]
        obs = np.asarray(r['obs'])
        d = r['draws']
        draws_store[name] = d
        title = titles.title.get(name + '.npy', name)
        cat = titles.category.get(name + '.npy', '?')

        mu = d.mean(0)
        sd = d.std(0, ddof=1)
        z = (obs - mu) / sd
        n_exc = (d >= obs[None, :]).sum(0)
        p_emp = (1 + n_exc) / (B + 1)          # assumption-free permutation p-value
        # exact (Clopper-Pearson) 95% interval for the exceedance probability
        ci_lo = np.array([_beta.ppf(0.025, c + 1, B - c) if c < B else 1.0 for c in n_exc])
        ci_hi = np.array([_beta.ppf(0.975, c + 1, B - c) if c < B else 1.0 for c in n_exc])

        T_obs = obs[:K_REPORT].sum()
        T_null = d[:, :K_REPORT].sum(1)
        p_topK = float((1 + (T_null >= T_obs).sum()) / (B + 1))
        maxT_obs = float(z[:K_REPORT].max())
        # max over ranks, taken WITHIN each replicate first -> controls the family-wise
        # error rate across the K_REPORT ranks tested simultaneously
        maxT_null = ((d[:, :K_REPORT] - mu[:K_REPORT]) / sd[:K_REPORT]).max(1)
        p_maxT = float((1 + (maxT_null >= maxT_obs).sum()) / (B + 1))

        # parallel analysis over the full spectrum, from the stored full-rank counts:
        # a rank "survives" if fewer than ALPHA of the B null draws reach the observed
        # value at that rank, i.e. the observation exceeds the null's own 95th
        # percentile at its rank. k_hat is the first rank that fails this (or the full
        # spectrum length if every rank passes). Equivalent in spirit to, but computed
        # differently from, the np.quantile-based k_hat in eigenvalue_permutation_test.py
        # (see the cross-script note in eigenvalue_permutation_test.py's log).
        q95_full = None
        above_full = r['n_above_full'] / B < ALPHA      # obs above null 95th pct at that rank
        k_hat = int(np.argmin(above_full)) if not above_full.all() else int(r['n_ranks_full'])

        rows.append(dict(
            sample=name, title=title, cat=cat, n=r['n'], p=r['p'], B=B,
            lambda_1=float(obs[0]), null_mean=float(mu[0]), null_sd=float(sd[0]),
            z_1=float(z[0]), n_exceed_1=int(n_exc[0]), p_empirical_1=float(p_emp[0]),
            p_ci_lo_1=float(ci_lo[0]), p_ci_hi_1=float(ci_hi[0]),
            null_max_1=float(d[:, 0].max()),
            topK_mass=float(T_obs), p_topK_mass=p_topK,
            maxT=maxT_obs, p_maxT=p_maxT, k_hat=k_hat,
            runtime_sec=r['runtime_sec']))

        for b in range(B):
            rank1_long.append((name, title, cat, b, float(d[b, 0])))
        for k in range(len(obs)):
            obs_long.append((name, title, cat, k + 1, float(obs[k])))
        for k in range(r['n_ranks_full']):
            moments_long.append((name, k + 1, float(r['mean_full'][k]),
                                 float(r['sd_full'][k]), float(r['min_full'][k]),
                                 float(r['max_full'][k]), int(r['n_above_full'][k])))

        per_sample[name] = dict(
            title=title, category=cat, n_cells=r['n'], n_genes=r['p'],
            dropped_tracker_genes=r['dropped'], validation=r['validation'],
            lambda_obs=obs[:K_REPORT].tolist(),
            null_mean=mu[:K_REPORT].tolist(), null_sd=sd[:K_REPORT].tolist(),
            z=z[:K_REPORT].tolist(), n_exceed=n_exc[:K_REPORT].tolist(),
            p_empirical=p_emp[:K_REPORT].tolist(),
            p_ci_lo=ci_lo[:K_REPORT].tolist(), p_ci_hi=ci_hi[:K_REPORT].tolist(),
            topK_mass_obs=float(T_obs), p_topK_mass=p_topK,
            maxT_obs=maxT_obs, p_maxT=p_maxT,
            k_hat_parallel_analysis=k_hat, n_ranks_full=r['n_ranks_full'],
            runtime_sec=r['runtime_sec'])

    df = pd.DataFrame(rows)

    # ------------------------------------------------------------ save
    np.savez_compressed(os.path.join(OUT, 'raw', tag + '_null_draws.npz'), **draws_store)
    pd.DataFrame(rank1_long, columns=['sample', 'title', 'category', 'replicate',
                                      'lambda_1_scrambled']).to_csv(
        os.path.join(OUT, 'raw', tag + '_rank1.csv'), index=False)
    pd.DataFrame(obs_long, columns=['sample', 'title', 'category', 'rank',
                                    'lambda_observed']).to_csv(
        os.path.join(OUT, 'raw', tag + '_observed.csv'), index=False)
    pd.DataFrame(moments_long, columns=['sample', 'rank', 'null_mean', 'null_sd',
                                        'null_min', 'null_max', 'n_exceed']).to_csv(
        os.path.join(OUT, 'raw', tag + '_null_moments.csv'), index=False)
    df.to_csv(os.path.join(OUT, 'raw', tag + '.csv'), index=False)

    gr = df[df.cat == 'r'].z_1.values
    gd = df[df.cat == 'd'].z_1.values
    u, p_mw = mannwhitneyu(gd, gr, alternative='less')
    group = dict(statistic='z_1', n_regulated=int(len(gr)), n_dysregulated=int(len(gd)),
                 median_z1_regulated=float(np.median(gr)),
                 median_z1_dysregulated=float(np.median(gd)),
                 mannwhitney_U=float(u), p_value=float(p_mw),
                 alternative='dysregulated < regulated (one-sided)')

    with open(os.path.join(OUT, 'raw', tag + '.txt'), 'w') as f:
        f.write('Permutation test, all datasets, B={}\n'.format(B) + '=' * 78 + '\n')
        f.write('run: {}   K_STORE={}   K_REPORT={}   seed base={}\n\n'.format(
            stamp, K_STORE, K_REPORT, SEED + 1000))
        f.write('p_empirical_1 = (1 + #{{lambda_1^perm >= lambda_1^obs}}) / (B + 1).\n'
                'No distributional assumption. Censored at 1/(B+1) = {:.6f}.\n'
                'p_ci_lo_1 / p_ci_hi_1 are the exact (Clopper-Pearson) 95% interval for\n'
                'the underlying exceedance probability given n_exceed_1 out of B.\n\n'
                'Preprocessing asserted identical to ev_data/*.npy for every sample.\n\n'
                .format(1 / (B + 1)))
        f.write(df.to_string(index=False, float_format=lambda x: '{:.6g}'.format(x)))
        f.write('\n\nGroup comparison (Mann-Whitney U on z_1, dysregulated < regulated):\n')
        for k, v in group.items():
            f.write('  {}: {}\n'.format(k, v))
        f.write('\n\nSaved for reuse (no recomputation needed to fit any distribution):\n')
        f.write('  {}_null_draws.npz    per-sample ({}, {}) scrambled eigenvalues\n'
                .format(tag, B, K_STORE))
        f.write('  {}_rank1.csv         {} x {} = {} rank-1 datapoints, tidy long form\n'
                .format(tag, len(samples), B, len(samples) * B))
        f.write('  {}_observed.csv      observed eigenvalues, ranks 1..{}\n'
                .format(tag, K_STORE))
        f.write('  {}_null_moments.csv  per-rank null mean/sd/min/max over the FULL '
                'spectrum\n'.format(tag))

    log = dict(
        experiment=tag, timestamp=stamp,
        description='Permutation test of the empirical eigenvalue spectrum against the '
                    'column-scrambled null, all datasets at B={}, reporting the '
                    'assumption-free empirical p-value. All scrambled eigenvalues for '
                    'ranks 1..{} are saved.'.format(B, K_STORE),
        null_hypothesis='genes are mutually independent (no gene-gene correlation)',
        null_construction='independent permutation of each gene (column) across cells; '
                          'reproduces src/analysis_functions.scramble',
        parameters=dict(B=B, K_STORE=K_STORE, K_REPORT=K_REPORT, alpha=ALPHA,
                        seed_base=SEED + 1000,
                        seed_per_sample={nm: seeds[nm] for nm in samples},
                        rng='numpy default_rng (PCG64), one independent stream per sample',
                        n_workers=n_workers,
                        p_value_formula='(1 + #{null >= obs}) / (B + 1)',
                        min_attainable_p=1 / (B + 1),
                        confidence_interval='exact Clopper-Pearson 95% on n_exceed/B'),
        saved_artifacts=dict(
            null_draws_npz=tag + '_null_draws.npz',
            rank1_csv=tag + '_rank1.csv',
            observed_csv=tag + '_observed.csv',
            null_moments_csv=tag + '_null_moments.csv',
            note='full spectra per replicate were not stored (~240 MB); ranks 1..{} are '
                 'stored draw-by-draw and all remaining ranks are summarised in the '
                 'moments file'.format(K_STORE)),
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         pandas=pd.__version__, platform=platform.platform()),
        datasets_tested=samples,
        group_comparison=group,
        total_runtime_sec=round(time.time() - t_start, 1),
        results=per_sample)
    with open(os.path.join(OUT, 'logs', tag + '.json'), 'w') as f:
        json.dump(log, f, indent=2)

    print('\n' + df[['title', 'cat', 'lambda_1', 'null_sd', 'z_1', 'n_exceed_1',
                     'p_empirical_1', 'p_ci_hi_1']].to_string(
        index=False, float_format=lambda x: '{:.6g}'.format(x)))
    print('\ntotal {:.0f}s -> results/permutation_test/{{raw,logs}}/{}.*'.format(
        time.time() - t_start, tag))
    return df


if __name__ == '__main__':
    b = int(sys.argv[1]) if len(sys.argv) > 1 else B_DEFAULT
    w = int(sys.argv[2]) if len(sys.argv) > 2 else None
    main(b, w)
