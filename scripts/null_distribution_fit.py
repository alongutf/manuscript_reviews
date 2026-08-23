"""Parametric tail probabilities for scrambled-eigenvalue null distributions.

The empirical permutation p-value, (1 + #{null >= obs}) / (B + 1), is censored at
1/(B+1) and discards the shape of the null. Here we instead fit a parametric law
to the B scrambled values of lambda_k and read the tail probability off the fit:

    p_k = SF_fitted(lambda_obs_k)

WHICH FAMILY?
-------------
For the largest eigenvalue of a Wishart-type null matrix the limiting law is
Tracy-Widom (beta=1), whose right tail decays as exp(-(2/3) s^{3/2}). That sits
BETWEEN a Gaussian, exp(-s^2/2), which decays faster, and a Gumbel, exp(-s),
which decays slower. Consequently:

    normal fit   -> tail too light  -> p too SMALL  (anti-conservative)
    gumbel fit   -> tail too heavy  -> p too LARGE  (conservative)

so the two bracket the Tracy-Widom answer. We fit a panel of families spanning
that range, select by AIC, and always report the min/max across families as an
honest uncertainty band. The 3-parameter gamma is included because the standard
Chiani (2014) approximation to the Tracy-Widom law is a shifted, scaled gamma,
so a free gamma fit contains a Tracy-Widom-shaped solution without our having to
hard-code its constants.

CAVEAT ON EXTRAPOLATION
-----------------------
A fit to B draws is only observationally constrained out to about the 1/B
quantile. Beyond that the reported p is an extrapolation whose value is set by
the assumed tail shape, not by data. Every result here therefore carries
`extrapolated` (p < 1/B) and the across-family spread; when that spread covers
many orders of magnitude, the correct statement is "p is below anything we can
resolve", not the numeral itself.
"""
import warnings

import numpy as np
from scipy import stats

# Families ordered from lightest to heaviest right tail.
CANDIDATES = {
    'norm': stats.norm,               # exp(-s^2/2)          -- lighter than Tracy-Widom
    'skewnorm': stats.skewnorm,       # Gaussian tail, skewed
    'gamma': stats.gamma,             # contains the Chiani approximation to Tracy-Widom
    'genextreme': stats.genextreme,   # GEV; shape decides the tail
    'gumbel_r': stats.gumbel_r,       # exp(-s)              -- heavier than Tracy-Widom
    'lognorm': stats.lognorm,         # heavier still
}

TAIL_ORDER = ['norm', 'skewnorm', 'gamma', 'genextreme', 'gumbel_r', 'lognorm']


def fit_null(x, candidates=None):
    """Fit each candidate family to the null draws x by MLE.

    Returns {name: dict(params, loglik, aic, ks_stat, ks_p, ad_stat)}, skipping
    families that fail to converge. The KS p-value is computed against the fitted
    parameters and is therefore optimistic (the same data set the parameters);
    it is a diagnostic for gross misfit, not a calibrated test.
    """
    x = np.asarray(x, dtype=float)
    out = {}
    for name in (candidates or TAIL_ORDER):
        dist = CANDIDATES[name]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                params = dist.fit(x)
                ll = float(np.sum(dist.logpdf(x, *params)))
                if not np.isfinite(ll):
                    continue
                ks = stats.kstest(x, dist.cdf, args=params)
            out[name] = dict(
                params=[float(v) for v in params],
                n_params=len(params),
                loglik=ll,
                aic=float(2 * len(params) - 2 * ll),
                ks_stat=float(ks.statistic),
                ks_p=float(ks.pvalue),
            )
        except Exception:
            continue
    return out


def tail_p(name, params, obs):
    """Upper-tail probability SF(obs) under a fitted family, with its log10.

    Uses logsf so that values far below double-precision underflow are still
    reported rather than collapsing to 0.
    """
    dist = CANDIDATES[name]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sf = float(dist.sf(obs, *params))
        logsf = float(dist.logsf(obs, *params))
    log10 = logsf / np.log(10) if np.isfinite(logsf) else -np.inf
    return sf, log10


def parametric_pvalue(null_draws, obs, B=None):
    """Full parametric tail assessment of `obs` against the null sample `null_draws`.

    Returns a dict with the AIC-selected family, its tail probability, the
    across-family band, and flags describing how far the answer is extrapolated
    beyond the support of the null sample.
    """
    x = np.asarray(null_draws, dtype=float)
    B = B or len(x)
    fits = fit_null(x)
    if not fits:
        return dict(error='no family converged')

    best = min(fits, key=lambda k: fits[k]['aic'])
    p_best, log10_best = tail_p(best, fits[best]['params'], obs)

    per_family = {}
    for name, f in fits.items():
        sf, l10 = tail_p(name, f['params'], obs)
        per_family[name] = dict(p=sf, log10_p=l10, aic=f['aic'],
                                delta_aic=f['aic'] - fits[best]['aic'],
                                ks_p=f['ks_p'])

    l10s = [v['log10_p'] for v in per_family.values() if np.isfinite(v['log10_p'])]
    mu, sd = float(x.mean()), float(x.std())
    z = (obs - mu) / sd if sd > 0 else np.nan
    p_emp = float((1 + np.sum(x >= obs)) / (B + 1))

    return dict(
        best_family=best,
        best_params=fits[best]['params'],
        best_aic=fits[best]['aic'],
        best_ks_p=fits[best]['ks_p'],
        p_parametric=p_best,
        log10_p_parametric=log10_best,
        log10_p_min=min(l10s) if l10s else None,      # lightest tail  -> smallest p
        log10_p_max=max(l10s) if l10s else None,      # heaviest tail  -> largest p
        family_spread_orders=(max(l10s) - min(l10s)) if len(l10s) > 1 else 0.0,
        p_empirical=p_emp,
        z=float(z),
        null_max=float(x.max()),
        obs_exceeds_null_max=bool(obs > x.max()),
        extrapolated=bool(p_best < 1.0 / B),
        per_family=per_family,
    )


def calibration_check(null_draws, n_fit=None, seed=0):
    """Validate the fitted CDF where data can actually check it.

    Split the null draws, fit on one half, and map the held-out half through the
    fitted CDF. If the family is right those PIT values are Uniform(0,1). This
    validates the fit only over the observable range (roughly the central
    1 - 2/n_fit of the distribution); it says nothing about the far tail, which
    is exactly the region the reported p-values extrapolate into.
    """
    x = np.asarray(null_draws, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(x))
    n_fit = n_fit or len(x) // 2
    x_fit, x_test = x[idx[:n_fit]], x[idx[n_fit:]]
    fits = fit_null(x_fit)
    if not fits:
        return dict(error='no family converged')
    best = min(fits, key=lambda k: fits[k]['aic'])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = CANDIDATES[best].cdf(x_test, *fits[best]['params'])
    ks = stats.kstest(u, 'uniform')
    return dict(family=best, n_fit=int(n_fit), n_test=int(len(x_test)),
                pit_ks_stat=float(ks.statistic), pit_ks_p=float(ks.pvalue),
                pit_mean=float(np.mean(u)))


# ---------------------------------------------------------------- tail-only fit
def gpd_tail_p(null_draws, obs, q_threshold=0.75):
    """Peaks-over-threshold tail probability via a Generalized Pareto fit.

    Fitting a family to the whole null sample and extrapolating is unreliable:
    the family is chosen by how it describes the BODY, which does not determine
    the tail. The Pickands-Balkema-de Haan theorem instead says that exceedances
    over a high threshold u converge to a Generalized Pareto distribution
    regardless of the parent law, so the tail is modelled on its own terms:

        P(X > x) = P(X > u) * SF_GPD(x - u),   x > u

    The fitted shape xi is what decides the answer:
        xi < 0  -> the null has a FINITE upper endpoint u - sigma/xi.
                   If obs lies beyond it, the null assigns it probability zero
                   and no finite p-value exists.
        xi = 0  -> exponential tail.
        xi > 0  -> polynomial (heavy) tail.
    """
    x = np.asarray(null_draws, dtype=float)
    u = float(np.quantile(x, q_threshold))
    exc = x[x > u] - u
    if len(exc) < 10:
        return dict(error='too few exceedances', n_exceed=int(len(exc)), threshold=u)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        xi, loc, sigma = stats.genpareto.fit(exc, floc=0)
        p_exceed = len(exc) / len(x)
        sf = float(stats.genpareto.sf(obs - u, xi, loc=0, scale=sigma))
        logsf = float(stats.genpareto.logsf(obs - u, xi, loc=0, scale=sigma))
    upper_endpoint = (u - sigma / xi) if xi < 0 else np.inf
    p = p_exceed * sf
    return dict(
        threshold=u, q_threshold=q_threshold, n_exceed=int(len(exc)),
        xi=float(xi), sigma=float(sigma), p_exceed=float(p_exceed),
        upper_endpoint=float(upper_endpoint),
        obs_beyond_endpoint=bool(np.isfinite(upper_endpoint) and obs > upper_endpoint),
        p=float(p),
        log10_p=float((np.log(p_exceed) + logsf) / np.log(10)) if np.isfinite(logsf) else -np.inf,
    )
