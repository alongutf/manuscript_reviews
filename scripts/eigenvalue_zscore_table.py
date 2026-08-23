"""z-scores and Gaussian tail probabilities for the scrambled-eigenvalue null.

No distribution is fitted. For each sample and each eigenvalue rank k we take the
B scrambled values of lambda_k, use their empirical mean and standard deviation,
and report

    z_k = (lambda_obs_k - mean_k) / sd_k
    p_k = P(X > lambda_obs_k)  with  X ~ Normal(mean_k, sd_k)
        = SF_normal(z_k)

Both moments come straight from the permutation draws; the only assumption is
that the null lambda_k is Gaussian.

NOTE ON THE ASSUMPTION
----------------------
The null largest eigenvalue is asymptotically Tracy-Widom, whose right tail
exp(-(2/3)s^{3/2}) is HEAVIER than the Gaussian exp(-s^2/2). A Gaussian tail
therefore understates the null's chance of producing a large lambda_1, so these
p-values are anti-conservative -- increasingly so the larger z_k is. The
`skewness`, `excess_kurtosis` and `normality_p` columns report how Gaussian the
null draws actually are, so the size of that approximation is visible. z_k itself
is assumption-free.

Reads the null draws saved by eigenvalue_permutation_test.py, so it recomputes
nothing.

Usage:
    python scripts/eigenvalue_zscore_table.py [run_tag]
"""
import os
import sys
import json
import glob
import platform
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(_REPO, 'results', 'permutation_test')
LOG10 = np.log(10)


def gaussian_tail(z):
    """Upper-tail probability of a standard normal, with a log10 that survives underflow."""
    return float(stats.norm.sf(z)), float(stats.norm.logsf(z) / LOG10)


def load_run(tag=None):
    """Load a permutation run's log and its saved null draws."""
    if tag is None:
        cands = sorted(glob.glob(os.path.join(OUT, 'logs',
                                              'eigenvalue_permutation_test_*.json')))
        if not cands:
            raise SystemExit('no permutation run found in ' + os.path.join(OUT, 'logs'))
        tag = os.path.basename(cands[-1])[:-5]
    with open(os.path.join(OUT, 'logs', tag + '.json')) as f:
        log = json.load(f)
    npz = os.path.join(OUT, 'raw', tag + '_null_draws.npz')
    if not os.path.exists(npz):
        raise SystemExit('missing null draws for ' + tag + '; rerun the permutation test')
    return tag, log, dict(np.load(npz))


def build(tag=None, label=''):
    tag, log, draws = load_run(tag)
    B = log['parameters']['B']
    rows = []
    for name, v in log['results'].items():
        d = np.asarray(draws[name])            # (B, K)
        obs = np.asarray(v['lambda_obs_topK'])
        K = d.shape[1]
        for k in range(K):
            col = d[:, k]
            mu, sd = float(col.mean()), float(col.std(ddof=1))
            z = (obs[k] - mu) / sd
            p, l10 = gaussian_tail(z)
            rows.append(dict(
                sample=name, title=v['title'], cat=v['category'], B=B, rank=k + 1,
                lambda_obs=float(obs[k]), null_mean=mu, null_sd=sd, z=float(z),
                p_gaussian=p, log10_p_gaussian=l10,
                null_min=float(col.min()), null_max=float(col.max()),
                skewness=float(stats.skew(col)),
                excess_kurtosis=float(stats.kurtosis(col)),
                normality_p=float(stats.normaltest(col).pvalue),
                # Monte Carlo error on sd propagates into z: sd(sd)/sd ~ 1/sqrt(2(B-1))
                z_mc_rel_err=float(1.0 / np.sqrt(2 * (B - 1))),
            ))
    df = pd.DataFrame(rows)
    df.attrs['tag'] = tag
    df.attrs['label'] = label
    df.attrs['B'] = B
    return df, log


def main(tag=None):
    df, log = build(tag)
    B = df.attrs['B']
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_tag = 'eigenvalue_zscore_gaussian_' + stamp

    # optional high-B rerun of the boundary samples, for the same statistic
    hi = sorted(glob.glob(os.path.join(OUT, 'logs', 'eigenvalue_permutation_highB_*.json')))
    df_hi = None
    if hi:
        hi_tag = os.path.basename(hi[-1])[:-5]
        hi_npz = os.path.join(OUT, 'raw', hi_tag + '_null_draws.npz')
        if os.path.exists(hi_npz):
            with open(os.path.join(OUT, 'logs', hi_tag + '.json')) as f:
                hlog = json.load(f)
            hdraws = dict(np.load(hi_npz))
            hrows = []
            for name, v in hlog['results'].items():
                col = np.asarray(hdraws[name])[:, 0]
                obs = v['lambda_obs_topK'][0]
                mu, sd = float(col.mean()), float(col.std(ddof=1))
                z = (obs - mu) / sd
                p, l10 = gaussian_tail(z)
                hrows.append(dict(sample=name, title=v['title'], B=hlog['parameters']['B'],
                                  lambda_obs=obs, null_mean=mu, null_sd=sd, z=z,
                                  p_gaussian=p, log10_p_gaussian=l10,
                                  p_empirical=v['p_empirical_topK'][0],
                                  normality_p=float(stats.normaltest(col).pvalue)))
            df_hi = pd.DataFrame(hrows)

    r1 = df[df['rank'] == 1].sort_values('z', ascending=False)

    df.to_csv(os.path.join(OUT, 'raw', out_tag + '.csv'), index=False)

    with open(os.path.join(OUT, 'raw', out_tag + '.txt'), 'w') as f:
        f.write('Eigenvalue z-scores and Gaussian tail probabilities\n' + '=' * 78 + '\n')
        f.write('source run: {}   B={}   ranks 1..{}\n\n'.format(
            df.attrs['tag'], B, df['rank'].max()))
        f.write('No distribution is fitted. For each rank k:\n'
                '  z_k = (lambda_obs_k - mean_k) / sd_k,  moments from the B scrambled draws\n'
                '  p_k = SF_normal(z_k) = P(X > lambda_obs_k), X ~ Normal(mean_k, sd_k)\n\n'
                'The null lambda_1 is asymptotically Tracy-Widom, whose right tail\n'
                'exp(-(2/3)s^1.5) is HEAVIER than the Gaussian exp(-s^2/2), so p_gaussian\n'
                'understates the null tail and is anti-conservative -- more so at large z.\n'
                'skewness / excess_kurtosis / normality_p describe how Gaussian the null\n'
                'draws actually are. z is assumption-free; only p_gaussian is not.\n'
                'sd is estimated from B draws, so z carries ~{:.1%} relative Monte Carlo\n'
                'error at B={}.\n\n'.format(1 / np.sqrt(2 * (B - 1)), B))

        f.write('RANK 1 (largest eigenvalue)\n' + '-' * 78 + '\n')
        cols = ['title', 'cat', 'lambda_obs', 'null_mean', 'null_sd', 'z',
                'p_gaussian', 'log10_p_gaussian', 'skewness', 'excess_kurtosis',
                'normality_p']
        f.write(r1[cols].to_string(index=False,
                                   float_format=lambda x: '{:.4g}'.format(x)))

        f.write('\n\n\nALL RANKS 1..{}: z (upper) and log10 p_gaussian (lower)\n'.format(
            df['rank'].max()) + '-' * 78 + '\n')
        for name, g in df.groupby('title', sort=False):
            g = g.sort_values('rank')
            f.write('\n{}  [{}]\n'.format(name, g.cat.iloc[0]))
            f.write('  k      : ' + ' '.join('{:>8d}'.format(k) for k in g['rank']) + '\n')
            f.write('  z      : ' + ' '.join('{:>8.2f}'.format(v) for v in g.z) + '\n')
            f.write('  log10 p: ' + ' '.join('{:>8.1f}'.format(v)
                                             for v in g.log10_p_gaussian) + '\n')

        if df_hi is not None:
            f.write('\n\n\nBOUNDARY SAMPLES AT HIGHER B (rank 1)\n' + '-' * 78 + '\n')
            f.write('p_empirical is the assumption-free exceedance count at that B, shown\n'
                    'so the Gaussian approximation can be checked where it is checkable.\n\n')
            f.write(df_hi.to_string(index=False,
                                    float_format=lambda x: '{:.5g}'.format(x)))
        f.write('\n')

    log_out = dict(
        experiment=out_tag, timestamp=stamp, source_run=df.attrs['tag'],
        description='z-scores and Gaussian tail probabilities of the observed eigenvalues '
                    'against the column-scrambled null. No distribution is fitted; the '
                    'null mean and sd are the empirical moments of the B permutation draws.',
        method=dict(
            z='z_k = (lambda_obs_k - mean_k) / sd_k',
            moments='empirical mean and sd (ddof=1) of the B scrambled lambda_k',
            p='p_k = SF_normal(z_k), i.e. P(X > lambda_obs_k) for X ~ Normal(mean_k, sd_k)',
            assumption='null lambda_k is Gaussian',
            assumption_caveat='the null largest eigenvalue is asymptotically Tracy-Widom, '
                              'right tail exp(-(2/3)s^{3/2}), heavier than the Gaussian '
                              'exp(-s^2/2); p_gaussian is therefore anti-conservative, '
                              'increasingly so at large z',
            diagnostics='skewness, excess_kurtosis and normaltest p on the null draws '
                        'quantify the departure from Gaussian',
            monte_carlo='sd from B draws carries relative error 1/sqrt(2(B-1)) = '
                        '{:.4f}, which propagates directly into z'.format(
                            1 / np.sqrt(2 * (B - 1)))),
        parameters=dict(B=B, ranks=int(df['rank'].max())),
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         pandas=pd.__version__, platform=platform.platform()),
        rank1=r1.to_dict(orient='records'),
        all_ranks=df.to_dict(orient='records'),
        high_B_rank1=(df_hi.to_dict(orient='records') if df_hi is not None else None),
    )
    with open(os.path.join(OUT, 'logs', out_tag + '.json'), 'w') as f:
        json.dump(log_out, f, indent=2)

    make_figure(df, r1, out_tag, B)
    print(r1[['title', 'cat', 'lambda_obs', 'null_mean', 'null_sd', 'z',
              'p_gaussian', 'log10_p_gaussian', 'normality_p']].to_string(
        index=False, float_format=lambda x: '{:.4g}'.format(x)))
    if df_hi is not None:
        print('\nboundary samples at higher B:')
        print(df_hi[['title', 'B', 'z', 'p_gaussian', 'p_empirical',
                     'normality_p']].to_string(index=False,
                                               float_format=lambda x: '{:.5g}'.format(x)))
    print('\nwrote results/permutation_test/{{raw,logs,figures}}/' + out_tag + '.*')
    return df


def make_figure(df, r1, out_tag, B):
    colors = {'r': '#2c6fbb', 'd': '#c8412f'}
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6),
                             gridspec_kw={'width_ratios': [1.4, 1.15, 1]})

    ax = axes[0]
    d = r1.sort_values(['cat', 'z'], ascending=[True, False]).reset_index(drop=True)
    ax.barh(range(len(d)), d.z, color=[colors.get(c, 'grey') for c in d.cat])
    ax.axvline(0, color='k', lw=0.8)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d.title, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xlim(-3, 30000)
    ax.set_xlabel(r'$z_1=(\lambda_1-\mu_{null})/\sigma_{null}$   (symlog)')
    ax.set_title('z-score of $\\lambda_1$ vs. scrambled null (B={})'.format(B), fontsize=10)
    for i, (z, l10) in enumerate(zip(d.z, d.log10_p_gaussian)):
        lab = 'p={:.2f}'.format(10 ** l10) if l10 > -3 else r'$p{=}10^{%.0f}$' % l10
        ax.text(z * 1.3 if z > 0 else 1.2, i, lab, va='center', ha='left', fontsize=6.5)

    ax = axes[1]
    for _, r in r1.iterrows():
        g = df[(df['sample'] == r['sample'])].sort_values('rank')
        ax.plot(g['rank'], g.z, marker='o', ms=2.5, lw=1,
                color=colors.get(r['cat'], 'grey'), alpha=0.85)
    ax.axhline(0, color='k', lw=0.8, ls=':')
    ax.set_yscale('symlog', linthresh=1)
    ax.set_ylim(-3, 3000)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.set_xlabel('eigenvalue rank $k$')
    ax.set_ylabel('$z_k$  (symlog)')
    ax.set_title('$z_k$ across ranks', fontsize=10)

    ax = axes[2]
    ax.scatter(r1.z, -r1.log10_p_gaussian, s=42,
               color=[colors.get(c, 'grey') for c in r1.cat], zorder=3)
    zz = np.linspace(max(r1.z.min(), 0.1), r1.z.max(), 200)
    ax.plot(zz, -stats.norm.logsf(zz) / LOG10, 'k-', lw=1,
            label=r'$-\log_{10}\overline{\Phi}(z)$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$z_1$')
    ax.set_ylabel('$-\\log_{10}\\,p_{gaussian}$')
    ax.set_title('Gaussian tail: $p$ is a deterministic\nfunction of $z$', fontsize=10)
    ax.legend(fontsize=8)

    fig.tight_layout()
    for ext in ('svg', 'png'):
        fig.savefig(os.path.join(OUT, 'figures', out_tag + '.' + ext),
                    dpi=200, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
