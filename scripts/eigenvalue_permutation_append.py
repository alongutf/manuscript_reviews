"""Add samples to an existing B=2000 permutation run, in place.

Runs the identical permutation pipeline on datasets that were not part of an
earlier run and merges them into that run's output files, so the artifacts stay
a single self-describing set rather than fragmenting across runs.

Seeding: the original run assigned seed = SEED + 1000 + i over its own sorted
sample list. Appending shifts those indices, so new samples are seeded from a
separate block (SEED + 2000 + i over the sorted list of ADDED samples). Existing
samples keep the seeds already recorded in the log, so every sample in the merged
run remains independently reproducible.

Files updated in place (same tag):
    raw/<tag>.csv                  summary row per sample
    raw/<tag>.txt                  human-readable report
    raw/<tag>_rank1.csv            15 -> 18 x B rank-1 datapoints
    raw/<tag>_observed.csv         observed eigenvalues, ranks 1..K_STORE
    raw/<tag>_null_moments.csv     per-rank null moments, full spectrum
    raw/<tag>_null_draws.npz       (B, K_STORE) scrambled eigenvalues per sample
    logs/<tag>.json                full log, with an `appended` provenance block

Usage:
    python scripts/eigenvalue_permutation_append.py <tag> [sample ...]
    # with no samples listed, every data_for_paper/*.csv missing from the run is added
"""
import os
import sys

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, os.environ.get('EV_BLAS_THREADS', '4'))

import json
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _REPO)

from eigenvalue_permutation_full_B2000 import _work, K_STORE, K_REPORT, SEED, ALPHA, OUT

DATA_DIR = os.path.join(_REPO, 'data_for_paper')


def summarise(name, r, titles, B):
    """Same summary computation as the original run, for one sample."""
    from scipy.stats import beta as _beta

    obs = np.asarray(r['obs'])
    d = r['draws']
    title = titles.title.get(name + '.npy', name)
    cat = titles.category.get(name + '.npy', '?')

    mu = d.mean(0)
    sd = d.std(0, ddof=1)
    z = (obs - mu) / sd
    n_exc = (d >= obs[None, :]).sum(0)
    p_emp = (1 + n_exc) / (B + 1)
    ci_lo = np.array([_beta.ppf(0.025, c + 1, B - c) if c < B else 1.0 for c in n_exc])
    ci_hi = np.array([_beta.ppf(0.975, c + 1, B - c) if c < B else 1.0 for c in n_exc])

    T_obs = obs[:K_REPORT].sum()
    T_null = d[:, :K_REPORT].sum(1)
    p_topK = float((1 + (T_null >= T_obs).sum()) / (B + 1))
    maxT_obs = float(z[:K_REPORT].max())
    maxT_null = ((d[:, :K_REPORT] - mu[:K_REPORT]) / sd[:K_REPORT]).max(1)
    p_maxT = float((1 + (maxT_null >= maxT_obs).sum()) / (B + 1))

    above_full = r['n_above_full'] / B < ALPHA
    k_hat = int(np.argmin(above_full)) if not above_full.all() else int(r['n_ranks_full'])

    row = dict(sample=name, title=title, cat=cat, n=r['n'], p=r['p'], B=B,
               lambda_1=float(obs[0]), null_mean=float(mu[0]), null_sd=float(sd[0]),
               z_1=float(z[0]), n_exceed_1=int(n_exc[0]), p_empirical_1=float(p_emp[0]),
               p_ci_lo_1=float(ci_lo[0]), p_ci_hi_1=float(ci_hi[0]),
               null_max_1=float(d[:, 0].max()),
               topK_mass=float(T_obs), p_topK_mass=p_topK,
               maxT=maxT_obs, p_maxT=p_maxT, k_hat=k_hat,
               runtime_sec=r['runtime_sec'])
    entry = dict(title=title, category=cat, n_cells=r['n'], n_genes=r['p'],
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
    return row, entry, title, cat, obs, d, r


def main(tag, new_samples=None, n_workers=None):
    """Compute the new samples' permutation draws and merge them into run `tag`'s
    existing output files in place, leaving already-run samples untouched."""
    from scipy.stats import mannwhitneyu

    log_path = os.path.join(OUT, 'logs', tag + '.json')
    if not os.path.exists(log_path):
        raise SystemExit('no such run: ' + log_path)
    with open(log_path) as f:
        log = json.load(f)
    B = log['parameters']['B']
    existing = list(log['results'].keys())

    all_csv = sorted(f[:-4] for f in os.listdir(DATA_DIR) if f.endswith('.csv'))
    new_samples = new_samples or [s for s in all_csv if s not in existing]
    new_samples = [s for s in new_samples if s not in existing]
    if not new_samples:
        print('nothing to add; run already contains all data_for_paper/*.csv')
        return
    missing = [s for s in new_samples if not os.path.exists(
        os.path.join(DATA_DIR, s + '.csv'))]
    if missing:
        raise SystemExit('no CSV for: ' + ', '.join(missing))

    # Fail fast: every target is rewritten in place, and a file held open by Excel
    # raises PermissionError. Check before spending the compute, not after.
    # ('r+b' opens for read/write without truncating, so this is a pure lock probe.)
    targets = [os.path.join(OUT, 'raw', tag + s) for s in
               ('.csv', '.txt', '_rank1.csv', '_observed.csv', '_null_moments.csv',
                '_null_draws.npz')] + [log_path]
    locked = []
    for t in targets:
        try:
            with open(t, 'r+b'):
                pass
        except OSError:
            locked.append(os.path.basename(t))
    if locked:
        raise SystemExit(
            'cannot write (file open in another program?):\n  ' + '\n  '.join(locked) +
            '\nClose it and rerun; nothing was computed.')

    n_workers = n_workers or min(len(new_samples), max(1, (os.cpu_count() or 4) // 4))
    print('run {}  B={}\nexisting: {} samples\nadding  : {}\nworkers : {}'.format(
        tag, B, len(existing), ', '.join(new_samples), n_workers), flush=True)

    # separate seed block (offset +2000, vs +1000 for the original run) so existing
    # samples keep the seeds already recorded in the log; `i` indexes only THIS call's
    # new_samples list, so a second, later append call reuses the same +2000+i seeds
    # for whatever it adds -- see the cross-append seed-collision note in FINDINGS
    seeds = {nm: SEED + 2000 + i for i, nm in enumerate(sorted(new_samples))}
    t_start = time.time()
    order = sorted(new_samples,
                   key=lambda nm: -os.path.getsize(os.path.join(DATA_DIR, nm + '.csv')))
    results = {}
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(_work, (nm, B, seeds[nm])): nm for nm in order}
        for fut in as_completed(futs):
            r = fut.result()
            results[r['name']] = r
            print('  done {:28s} n={:5d}  ({:.0f}s)'.format(
                r['name'], r['n'], r['runtime_sec']), flush=True)

    # Cache the finished permutations before touching any output file, so a write
    # failure costs a retry of the merge and not a rerun of the compute.
    cache = os.path.join(OUT, 'raw', '.append_cache_' + tag + '.npz')
    np.savez_compressed(cache, **{nm: results[nm]['draws'] for nm in results},
                        **{nm + '::obs': np.asarray(results[nm]['obs']) for nm in results},
                        **{nm + '::meanfull': results[nm]['mean_full'] for nm in results},
                        **{nm + '::sdfull': results[nm]['sd_full'] for nm in results},
                        **{nm + '::minfull': results[nm]['min_full'] for nm in results},
                        **{nm + '::maxfull': results[nm]['max_full'] for nm in results},
                        **{nm + '::abovefull': results[nm]['n_above_full'] for nm in results})
    print('  cached permutations -> ' + os.path.basename(cache), flush=True)
    # (the '::'-suffixed keys pack each sample's per-rank summary arrays alongside its
    # draws in one npz, since savez_compressed only takes flat name -> array pairs)

    titles = pd.read_excel(os.path.join(_REPO, 'ev_data',
                                        'titles.xlsx')).set_index('file_name')

    # ---------------------------------------------------------------- merge
    df = pd.read_csv(os.path.join(OUT, 'raw', tag + '.csv'))
    rank1 = pd.read_csv(os.path.join(OUT, 'raw', tag + '_rank1.csv'))
    obs_df = pd.read_csv(os.path.join(OUT, 'raw', tag + '_observed.csv'))
    mom_df = pd.read_csv(os.path.join(OUT, 'raw', tag + '_null_moments.csv'))
    draws_store = dict(np.load(os.path.join(OUT, 'raw', tag + '_null_draws.npz')))

    new_rows, new_r1, new_obs, new_mom = [], [], [], []
    for name in sorted(new_samples):
        row, entry, title, cat, obs, d, r = summarise(name, results[name], titles, B)
        new_rows.append(row)
        log['results'][name] = entry
        draws_store[name] = d
        for b in range(B):
            new_r1.append((name, title, cat, b, float(d[b, 0])))
        for k in range(len(obs)):
            new_obs.append((name, title, cat, k + 1, float(obs[k])))
        for k in range(r['n_ranks_full']):
            new_mom.append((name, k + 1, float(r['mean_full'][k]), float(r['sd_full'][k]),
                            float(r['min_full'][k]), float(r['max_full'][k]),
                            int(r['n_above_full'][k])))

    df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
    df = df.sort_values('sample').reset_index(drop=True)
    rank1 = pd.concat([rank1, pd.DataFrame(
        new_r1, columns=['sample', 'title', 'category', 'replicate',
                         'lambda_1_scrambled'])], ignore_index=True)
    obs_df = pd.concat([obs_df, pd.DataFrame(
        new_obs, columns=['sample', 'title', 'category', 'rank',
                          'lambda_observed'])], ignore_index=True)
    mom_df = pd.concat([mom_df, pd.DataFrame(
        new_mom, columns=['sample', 'rank', 'null_mean', 'null_sd', 'null_min',
                          'null_max', 'n_exceed'])], ignore_index=True)

    gr = df[df.cat == 'r'].z_1.values
    gd = df[df.cat == 'd'].z_1.values
    u, p_mw = mannwhitneyu(gd, gr, alternative='less')
    group = dict(statistic='z_1', n_regulated=int(len(gr)), n_dysregulated=int(len(gd)),
                 median_z1_regulated=float(np.median(gr)),
                 median_z1_dysregulated=float(np.median(gd)),
                 mannwhitney_U=float(u), p_value=float(p_mw),
                 alternative='dysregulated < regulated (one-sided)')

    # ---------------------------------------------------------------- write
    df.to_csv(os.path.join(OUT, 'raw', tag + '.csv'), index=False)
    rank1.to_csv(os.path.join(OUT, 'raw', tag + '_rank1.csv'), index=False)
    obs_df.to_csv(os.path.join(OUT, 'raw', tag + '_observed.csv'), index=False)
    mom_df.to_csv(os.path.join(OUT, 'raw', tag + '_null_moments.csv'), index=False)
    np.savez_compressed(os.path.join(OUT, 'raw', tag + '_null_draws.npz'), **draws_store)

    samples_all = sorted(log['results'].keys())
    with open(os.path.join(OUT, 'raw', tag + '.txt'), 'w') as f:
        f.write('Permutation test, all datasets, B={}\n'.format(B) + '=' * 78 + '\n')
        f.write('run: {}   K_STORE={}   K_REPORT={}\n'.format(tag, K_STORE, K_REPORT))
        f.write('samples: {}\n\n'.format(len(samples_all)))
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
                .format(tag, len(samples_all), B, len(samples_all) * B))
        f.write('  {}_observed.csv      observed eigenvalues, ranks 1..{}\n'
                .format(tag, K_STORE))
        f.write('  {}_null_moments.csv  per-rank null mean/sd/min/max over the FULL '
                'spectrum\n'.format(tag))

    log['datasets_tested'] = samples_all
    log['group_comparison'] = group
    log['parameters']['seed_per_sample'].update(seeds)
    log.setdefault('appended', []).append(dict(
        timestamp=datetime.now().strftime('%Y%m%d_%H%M%S'),
        samples_added=sorted(new_samples), B=B, n_workers=n_workers,
        seed_block='SEED + 2000 + i over the sorted added samples',
        seeds={k: int(v) for k, v in seeds.items()},
        runtime_sec=round(time.time() - t_start, 1),
        note='existing samples were not recomputed; their seeds and results are '
             'unchanged from the original run'))
    with open(log_path, 'w') as f:
        json.dump(log, f, indent=2)

    if os.path.exists(cache):        # merge succeeded; the resume cache is now dead weight
        os.remove(cache)

    print('\n' + df[['title', 'cat', 'lambda_1', 'null_sd', 'z_1', 'n_exceed_1',
                     'p_empirical_1', 'p_ci_hi_1']].to_string(
        index=False, float_format=lambda x: '{:.6g}'.format(x)))
    print('\nmerged {} new samples into {} ({} total) in {:.0f}s'.format(
        len(new_samples), tag, len(samples_all), time.time() - t_start))
    return df


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2:] or None)
