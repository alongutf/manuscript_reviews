"""Shared permutation p-value annotation for the eigenvalue-spectrum panels.

Used by figure2.py, figure3.py, figure5.py and supplementary_figures/figure_s5.py
so the indicator is worded and formatted identically everywhere.

Source of truth: results/data_metrics/data_metrics.csv, column `permutation_p`
-- the empirical p-value from the B=2000 column-scramble permutation test
(scripts/eigenvalue_permutation_full_B2000.py):

    p = (1 + #{lambda_1^perm >= lambda_1^obs}) / (B + 1)

With B = 2000 this is censored below at 1/2001 = 4.9975e-4, so a sample where no
permutation reached the observed lambda_1 is reported as "p < 5x10^-4" rather
than as a specific number it cannot resolve.

Simulated datasets (ev_data/simulated_pcs_*.npy) have no permutation test, so
`perm_p` returns None for them and the annotation is silently skipped.
"""
import os

import numpy as np
import pandas as pd
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
METRICS = os.path.join(_REPO, 'results', 'data_metrics', 'data_metrics.csv')

# permutation replicates behind data_metrics.csv:permutation_p
B = 2000
P_FLOOR = 1.0 / (B + 1)

_CACHE = {}


def _stem(name):
    """Dataset key from a bare name or a full path, with any .npy/.csv stripped."""
    name = os.path.basename(name)
    for ext in ('.npy', '.csv'):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def _table():
    if 'p' not in _CACHE:
        d = pd.read_csv(METRICS, index_col=0)
        if 'permutation_p' not in d.columns:
            raise RuntimeError(
                'data_metrics.csv has no `permutation_p` column; run '
                'scripts/add_permutation_metrics.py first')
        _CACHE['p'] = {_stem(f): float(p) for f, p in
                       zip(d['file_name'], d['permutation_p'])}
        ci = (d['gmp_cor_ci'] if 'gmp_cor_ci' in d.columns
              else pd.Series(np.nan, index=d.index))
        _CACHE['ci'] = {_stem(f): float(c) for f, c in zip(d['file_name'], ci)}
    return _CACHE['p']


def perm_p(dataset):
    """Empirical permutation p for a dataset, or None if it was not tested.

    `dataset` may be None (panels drawn from data with no underlying dataset, such
    as figure2's simulation inset); those simply get no p-value.
    """
    if dataset is None:
        return None
    return _table().get(_stem(dataset))


def gmp_cor_ci(dataset):
    """sqrt(N)*sigma uncertainty on GMP-Cor, or None if unavailable."""
    _table()
    v = _CACHE['ci'].get(_stem(dataset))
    return None if v is None or not np.isfinite(v) else v


def p_label(p):
    """Format the p-value.

    At or below the B=2000 resolution floor the p-value is not resolved, so it is
    reported as an upper bound rather than as the floor value itself.
    """
    if p is None:
        return None
    if p <= P_FLOOR * 1.001:
        return r'$p<5\times10^{-4}$'
    if p < 0.01:
        return r'$p={:.3f}$'.format(p)
    return r'$p={:.2f}$'.format(p)


class _ZeroWidthHandle(HandlerBase):
    """Legend handler that draws nothing and occupies no width.

    Collapsing the handle column lets the p-value text start at the left edge of
    the legend box, flush with the series' line handles, instead of being indented
    into the label column with the series names.
    """

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        handlebox.set_width(0)
        return Line2D([], [], linestyle='none')


def legend_with_p(ax, dataset, p_fontsize=None, **legend_kwargs):
    """Redraw the axis legend with the permutation p-value as a final entry.

    The p-value belongs with the signal/scrambled keys it qualifies, so it goes
    inside the legend box rather than floating as a separate annotation: left
    aligned to the box edge and in a smaller font than the series labels.

    Datasets without a permutation test (the simulated spectra) get an ordinary
    legend, unchanged.
    """
    handles, labels = ax.get_legend_handles_labels()
    lab = p_label(perm_p(dataset))
    if lab is None:
        return ax.legend(handles, labels, **legend_kwargs)

    ph = Line2D([], [], linestyle='none', marker='')
    handler_map = dict(legend_kwargs.pop('handler_map', None) or {})
    handler_map[ph] = _ZeroWidthHandle()
    leg = ax.legend(list(handles) + [ph], list(labels) + [lab],
                    handler_map=handler_map, **legend_kwargs)

    if p_fontsize is None:
        base = legend_kwargs.get('fontsize')
        base = leg.get_texts()[0].get_fontsize() if base is None else base
        p_fontsize = max(4.0, float(base) - 2)
    leg.get_texts()[-1].set_fontsize(p_fontsize)
    return leg


def gmp_cor_label(dataset, gmp_cor, decimals=1):
    """'GMP-Cor: 13.8 +/- 0.2' -- the value with its sqrt(N)*sigma uncertainty."""
    ci = gmp_cor_ci(dataset)
    fmt = '{:.' + str(decimals) + 'f}'
    if ci is None:
        return ('GMP-Cor: ' + fmt).format(gmp_cor)
    return ('GMP-Cor: ' + fmt + r' $\pm$ ' + fmt).format(gmp_cor, ci)
