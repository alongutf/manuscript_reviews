"""
Supplementary Figure S8 — PC9 cancer cell-line analysis.

Applies the correlation-spectrum / GMP-Cor framework (developed for the
bacterial datasets) to day-14 PC9 scRNA-seq samples, showing that the same
global-dysregulation readout transfers to a mammalian cancer line.

Panels (2x2)
------------
A  Representative correlation-spectrum CCDF for the Proliferating subpopulation,
   drawn in the same loglog style as figure2.py / figure3.py, annotated GMP-Cor.
B  Same, for the Arrested subpopulation.
C  UMAP of all PC9 cells (precomputed coords in pc9_data/obs_metadata.csv),
   colored by proliferation state (Proliferating / Intermediate / Arrested).
D  GMP-Cor per subpopulation, averaged over the two biological replicates;
   error bars are the standard error of the mean (n = 2).

Note: low/med/high in the data refer to a fluorescent-marker level; proliferation
is the opposite (low marker = Proliferating, high marker = Arrested). Only the
proliferation names are shown to readers.

Inputs (all produced by scripts/cancer_pc9_correlation_analysis.py, except the
metadata which ships with the repo):
  pc9_data/obs_metadata.csv                        -> panel C
  results/cancer_pc9/ev_data/<sample>.npy          -> panels A & B  (shape (2, P))
  results/cancer_pc9/pc9_day14_gmp_cor.csv         -> panels A, B & D

Run from this directory:
    cd scripts/supplementary_figures
    python figure_s8.py
Writes figure_s8.pdf + figure_s8_preview.png next to this script.
"""
import os
import sys

# --- import bootstrap: this script lives in scripts/supplementary_figures/ ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                       # repo root
sys.path.insert(0, _REPO)                                            # import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))        # figure_functions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Global style
# ------------------------------------------------------------------
fsize = 10
plt.close("all")
plt.style.use('default')

OBS_META = os.path.join(_REPO, 'pc9_data', 'obs_metadata.csv')
EV_DATA_DIR = os.path.join(_REPO, 'results', 'cancer_pc9', 'ev_data')
METRICS_CSV = os.path.join(_REPO, 'results', 'cancer_pc9', 'pc9_day14_gmp_cor.csv')

# Subpopulation color scheme (shared across panels A/B/C).
# low/med/high refer to the *fluorescent marker* level; the proliferation state
# is the OPPOSITE (low marker = proliferating/cycling, high marker = arrested).
COND_COLOR = {'low': 'steelblue', 'med': '#fec44f', 'high': 'salmon'}
COND_ORDER = ['low', 'med', 'high']
COND_DISPLAY = {'low': 'Proliferating', 'med': 'Intermediate', 'high': 'Arrested'}

# Leiden cluster -> proliferation annotation (drawn above each cluster in panel A)
CLUSTER_LABEL = {2: 'Cycling', 0: 'Non-cycling'}


def _condition(sample_tag):
    """Map a sample tag to its subpopulation: 'high1'/'14_rep1_high' -> 'high'."""
    for cond in COND_ORDER:
        if cond in sample_tag:
            return cond
    return None


# ==================================================================
# Panel A — UMAP colored by subpopulation, clusters annotated
# ==================================================================
def panel_A(ax):
    obs = pd.read_csv(OBS_META, index_col=0)
    obs['condition'] = obs['sample'].map(_condition)

    # draw high (densest) first so the sparser groups stay visible on top
    for cond in ['high', 'med', 'low']:
        sub = obs[obs['condition'] == cond]
        ax.scatter(sub['UMAP_X'], sub['UMAP_Y'], s=2, linewidths=0,
                   color=COND_COLOR[cond], alpha=0.5, rasterized=True,
                   label=COND_DISPLAY[cond])


    ax.set_xlabel('UMAP 1', fontsize=fsize - 2, labelpad=1)
    ax.set_ylabel('UMAP 2', fontsize=fsize - 2, labelpad=1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
    leg = ax.legend(fontsize=fsize - 2, markerscale=4, loc='lower left',
                    frameon=True, framealpha=0.9)
    leg.get_title().set_fontsize(fsize - 2)


# ==================================================================
# Panel B — representative correlation-spectrum CCDFs
# ==================================================================
def _plot_ccdf(ax, sample, title, gmp_cor, signal_color, show_ylabel=True):
    """loglog CCDF of empirical (row 0) vs scrambled (row 1) eigenvalues.

    Color scheme shared with figure2.py / figure3.py:
      grey  -> spurious correlations (below scrambled max)
      color -> true signal (at/above scrambled max)
      black -> scrambled
    """
    arr = np.load(os.path.join(EV_DATA_DIR, f'{sample}.npy'))
    data1 = arr[0];  data1 = data1[data1 > 0]
    data2 = arr[1];  data2 = data2[data2 > 0]
    x2 = float(np.max(data2))                       # scrambled max = GMP-Cor threshold

    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1, p2 = len(d1s), len(d2s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    noise = d1s < x2

    ax.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='spurious', markersize=3)
    ax.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
              color=signal_color, label='signal', markersize=3)
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled', markersize=3)

    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
    ax.text(x2 * 1.1, 0.8, r'$\lambda_\mathrm{max}^\mathrm{scr}$',
            fontsize=fsize - 2, va='center', ha='left', color='k', alpha=0.7,
            transform=ax.get_xaxis_transform())
    ax.text(0.04, 0.05, f'GMP-Cor: {gmp_cor:.1f}', transform=ax.transAxes,
            fontsize=fsize - 2, fontweight='bold', va='bottom', ha='left')

    ax.set_xlim([0.1, 200])
    ax.set_ylim(top=1.5)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    if show_ylabel:
        ax.set_ylabel(r'CCDF: $P(X>\lambda)$', fontsize=fsize - 2, labelpad=0)
    ax.set_title(title, fontsize=fsize - 1)
    ax.legend(fontsize=fsize - 3, loc='upper right')
    ax.tick_params(labelsize=fsize - 2)


_metrics = pd.read_csv(METRICS_CSV).set_index('sample')
_gmp = _metrics['GMP_Cor']

# representative pair: highest-cycling (high) vs lowest-cycling (low), rep 1
REP_HIGH = '14_rep1_high'
REP_LOW = '14_rep1_low'


def panel_proliferating_ccdf(ax):
    _plot_ccdf(ax, REP_LOW, 'Proliferating', _gmp[REP_LOW],
               COND_COLOR['low'], show_ylabel=True)


def panel_arrested_ccdf(ax):
    _plot_ccdf(ax, REP_HIGH, 'Arrested', _gmp[REP_HIGH],
               COND_COLOR['high'], show_ylabel=False)


# ==================================================================
# Panel C — GMP-Cor per subpopulation, mean +/- SEM over 2 replicates
# ==================================================================
def panel_C(ax):
    m = _metrics.reset_index().copy()
    m['condition'] = m['sample'].map(_condition)

    means, sems, colors = [], [], []
    for cond in COND_ORDER:
        vals = m.loc[m['condition'] == cond, 'GMP_Cor'].to_numpy()
        means.append(vals.mean())
        sems.append(vals.std(ddof=1) / np.sqrt(len(vals)))   # SEM, n = 2
        colors.append(COND_COLOR[cond])

    x = np.arange(len(COND_ORDER))
    ax.bar(x, means, yerr=sems, capsize=4, color=colors, alpha=0.85,
           edgecolor='black', linewidth=0.8, width=0.3, zorder=2)

    # overlay the individual replicate points
    for i, cond in enumerate(COND_ORDER):
        vals = m.loc[m['condition'] == cond, 'GMP_Cor'].to_numpy()
        ax.scatter(np.full(len(vals), i), vals, color='black', s=14,
                   zorder=3, alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels([COND_DISPLAY[c] for c in COND_ORDER],
                       fontsize=fsize - 2)
    ax.set_ylabel('GMP-Cor', fontsize=fsize - 2)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax.set_ylim(bottom=50)
    ax.grid(axis='y', linestyle='--', alpha=0.3)


# ------------------------------------------------------------------
# Assemble
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 5), label_offset=(-0.05, 0.03))

# Row 1 — representative CCDFs (Proliferating | Arrested)
pf.add_panel([0.08, 0.6, 0.36, 0.32], label='A', draw_func=panel_proliferating_ccdf)
pf.add_panel([0.56, 0.6, 0.36, 0.32], label='B', draw_func=panel_arrested_ccdf)
# Row 2 — UMAP (smaller) | GMP-Cor bar plot
pf.add_panel([0.08, 0.08, 0.42, 0.4], label='C', draw_func=panel_A)
pf.add_panel([0.62, 0.08, 0.3, 0.4], label='D', draw_func=panel_C)

pf.save('figure_s8.pdf', dpi=300, transparent=True)
pf.save('figure_s8_preview.png', dpi=200)
print('Saved figure_s8.pdf and figure_s8_preview.png')
plt.show()
