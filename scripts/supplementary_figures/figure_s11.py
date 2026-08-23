import os
import sys

# --- import bootstrap: this script lives in scripts/supplementary_figures/ ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                       # repo root
sys.path.insert(0, _REPO)                                             # import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))         # figure_functions
sys.path.insert(0, os.path.join(_REPO, 'scripts'))                    # rpos_regulon_deseq

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figure_functions import PanelFigure
from rpos_regulon_deseq import DESEQ, REGULON, load_regulon

# ------------------------------------------------------------------
# Supplementary Figure S11
#   Volcano plots for three DESeq2 contrasts with the RegulonDB sigma-38
#   (RpoS) sigmulon highlighted, and a Fisher's exact test of regulon
#   over/under-representation in the significant lobes annotated on each panel.
#
#   A. Dis-Arrest  (results/deseq_results/from counts/deseq2_results_disrupted.csv)
#   B. Reg-Arrest  (results/deseq_results/from counts/deseq2_results_regulated.csv)
#   C. Early VapC  (aggregated scRNA-seq, VapC early vs exponential;
#      results/deseq_results/aggregated_sc/deseq2_results_vapc-early_vs_exp.csv)
#
#   Method follows scripts/rpos_regulon_volcano_fisher.py: lobes are
#   |log2FC| > 1 AND padj < ALPHA, and Fisher's exact 2x2 is regulon vs all other
#   genes tested in that contrast (padj not NA), one-sided in the direction the
#   odds ratio points. Unlike that script, ALPHA here is 0.05 for all three
#   panels -- see the note on ALPHA below.
# ------------------------------------------------------------------

fsize = 10
LFC_CUT = 1.0

# One significance threshold for all three panels, so the lobes -- and therefore the
# odds ratios and lobe sizes -- are directly comparable across them. This overrides
# ALPHA_BY_FOLDER in scripts/rpos_regulon_deseq.py, which holds aggregated_sc to 0.01
# on dispersion grounds; that stricter cutoff still governs the 16-panel diagnostic
# figure and documents/rpoS_regulon_deseq_analysis.md. Tightening Early VapC back to
# 0.01 moves its up-lobe odds ratio only 0.60 -> 0.67, so the panel's conclusion does
# not depend on this choice.
ALPHA = 0.05

SIG_COLOR = '#5BC8DC'    # significant genes (light cyan)
NS_COLOR = '#c8c8c8'     # non-significant genes (grey)
REG_COLOR = '#7B3FA0'    # sigma-38 / rpoS regulon genes (purple)

CONTRASTS = [
    ('from counts/deseq2_results_disrupted.csv', 'Dis-Arrest'),
    ('from counts/deseq2_results_regulated.csv', 'Reg-Arrest'),
    ('aggregated_sc/deseq2_results_vapc-early_vs_exp.csv', 'Early VapC'),
]

REGULON_GENES, REGULON_HEADER = load_regulon(REGULON)


def _fisher(in_regulon, lobe):
    """2x2 Fisher for one lobe: regulon vs all other tested genes."""
    a = int((in_regulon & lobe).sum())
    b = int((in_regulon & ~lobe).sum())
    c = int((~in_regulon & lobe).sum())
    d = int((~in_regulon & ~lobe).sum())
    table = [[a, b], [c, d]]
    odds, _ = stats.fisher_exact(table, alternative='two-sided')
    _, p_greater = stats.fisher_exact(table, alternative='greater')
    _, p_less = stats.fisher_exact(table, alternative='less')
    return dict(n=a, lobe_size=a + c, odds_ratio=odds,
                p_greater=p_greater, p_less=p_less)


def _fmt_p(p):
    """Two decimals; one decimal plus an explicit decade once p < 0.01."""
    if p >= 0.01:
        return 'p = %.2f' % p
    exp = int(np.floor(np.log10(p)))
    return r'p = %.1f$\times$10$^{%d}$' % (p / 10 ** exp, exp)


def _lobe_call(res):
    """(one-sided p, box text) for a lobe, phrased in the direction the OR points."""
    if res['lobe_size'] == 0 or not np.isfinite(res['odds_ratio']):
        return np.inf, 'rpoS regulon\nnot testable'
    enriched = res['odds_ratio'] >= 1
    word = 'enriched' if enriched else 'depleted'
    p = res['p_greater'] if enriched else res['p_less']
    return p, r'rpoS regulon $\bf{%s}$' '\nOR = %.2f\n%s' % (word, res['odds_ratio'], _fmt_p(p))


def _best_lobe(up, down):
    """Annotate only the lobe with the stronger one-sided result. Returns (name, text)."""
    p_up, txt_up = _lobe_call(up)
    p_down, txt_down = _lobe_call(down)
    return ('up', txt_up) if p_up <= p_down else ('down', txt_down)


def _volcano(ax, relpath, title):
    path = os.path.join(DESEQ, *relpath.split('/'))
    alpha = ALPHA

    d = pd.read_csv(path, index_col=0)
    d = d[d.padj.notna() & d.log2FoldChange.notna()]

    in_regulon = pd.Series(d.index.isin(REGULON_GENES), index=d.index)
    up = (d.log2FoldChange > LFC_CUT) & (d.padj < alpha)
    down = (d.log2FoldChange < -LFC_CUT) & (d.padj < alpha)
    sig = up | down

    # padj is exactly 0 for a few genes; floor it at the smallest nonzero padj
    # in this contrast so the y axis stays honest instead of hitting 300.
    pad = d.padj.values.copy()
    nz = pad[pad > 0]
    floor = nz.min() if nz.size else 1e-300
    y = pd.Series(-np.log10(np.where(pad > 0, pad, floor)), index=d.index)

    x = d.log2FoldChange

    # non-significant genes are grey whether or not they are in the regulon;
    # purple is reserved for regulon genes that actually clear both cutoffs.
    grey = ~sig
    cyan = sig & ~in_regulon
    purple = sig & in_regulon
    ax.scatter(x[grey], y[grey], s=4, c=NS_COLOR, alpha=0.45,
               linewidths=0, rasterized=True)
    ax.scatter(x[cyan], y[cyan], s=4, c=SIG_COLOR, alpha=0.55,
               linewidths=0, rasterized=True)
    ax.scatter(x[purple], y[purple], s=9, c=REG_COLOR, alpha=0.75,
               linewidths=0, zorder=3)

    ax.axhline(-np.log10(alpha), color='k', lw=0.6, ls='--', alpha=0.5)
    ax.axvline(LFC_CUT, color='k', lw=0.6, ls='--', alpha=0.5)
    ax.axvline(-LFC_CUT, color='k', lw=0.6, ls='--', alpha=0.5)

    res_up = _fisher(in_regulon, up)
    res_down = _fisher(in_regulon, down)

    ax.set_title(title, fontsize=fsize)
    ax.set_xlabel(r'log$_2$ fold change', fontsize=fsize - 2)
    ax.set_ylabel(r'-log$_{10}$ padj', fontsize=fsize - 2)
    ax.set_yticks([0, 100, 200, 300])
    ax.tick_params(axis='both', which='major', labelsize=fsize - 3)

    # headroom for the annotation box
    ax.set_ylim(top=350)
    # one box only, on the side of the lobe with the stronger result
    lobe_name, lobe_txt = _best_lobe(res_up, res_down)
    x_pos, halign = (0.97, 'right') if lobe_name == 'up' else (0.03, 'left')
    ax.text(x_pos, 0.985, lobe_txt, transform=ax.transAxes,
            ha=halign, va='top', multialignment='left', fontsize=fsize - 4,
            linespacing=1.35,
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                      edgecolor='black', linewidth=0.7, alpha=0.92))

    return dict(contrast=relpath, alpha=alpha, genes_tested=len(d),
                regulon_tested=int(in_regulon.sum()),
                up=res_up, down=res_down)


summaries = []


def panel_A(ax):
    summaries.append(_volcano(ax, *CONTRASTS[0]))


def panel_B(ax):
    summaries.append(_volcano(ax, *CONTRASTS[1]))


def panel_C(ax):
    summaries.append(_volcano(ax, *CONTRASTS[2]))


# ------------------------------------------------------------------
# Assemble - single row, three columns
# ------------------------------------------------------------------
plt.close('all')
pf = PanelFigure(figsize=(7, 3.0), label_offset=(-0.045, 0.05))
pf.add_panel([0.070, 0.28, 0.245, 0.62], label='A', draw_func=panel_A)
pf.add_panel([0.395, 0.28, 0.245, 0.62], label='B', draw_func=panel_B)
pf.add_panel([0.720, 0.28, 0.245, 0.62], label='C', draw_func=panel_C)

handles = [
    Line2D([], [], marker='o', ls='', color=SIG_COLOR, markersize=5,
           label='significant'),
    Line2D([], [], marker='o', ls='', color=NS_COLOR, markersize=5,
           label='not significant'),
    Line2D([], [], marker='o', ls='', color=REG_COLOR, markersize=5,
           label='significant rpoS regulon gene (RegulonDB %s)'
                 % REGULON_HEADER.get('RegulonDB release', '').split()[0]),
]
pf.fig.legend(handles=handles, loc='lower center', fontsize=fsize - 2, ncol=3,
              bbox_to_anchor=(0.5, 0.005), frameon=False)

pf.save(os.path.join(_HERE, 'figure_s11.pdf'), dpi=300)
pf.save(os.path.join(_HERE, 'figure_s11_preview.png'), dpi=200)
print('Saved figure_s11.pdf and figure_s11_preview.png')

for s in summaries:
    print('%-52s alpha=%.2g  tested=%d  regulon=%d' %
          (s['contrast'], s['alpha'], s['genes_tested'], s['regulon_tested']))
    for lobe in ('up', 'down'):
        r = s[lobe]
        print('    %-5s lobe n=%4d  regulon=%3d  OR=%.3f  p_enr=%.3g  p_dep=%.3g'
              % (lobe, r['lobe_size'], r['n'], r['odds_ratio'],
                 r['p_greater'], r['p_less']))
