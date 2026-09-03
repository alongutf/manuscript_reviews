"""
Standalone render of Supplementary Figure S6, panel A only.

Cell-area violin plots (from microscopy segmentation) comparing Dis-Arrest
(SHX), VapC 24h and Reg-Arrest, output as a narrower, transparent-background
SVG (the other panels of figure S6 are assembled elsewhere; this script exists
so panel A alone can be re-rendered/resized independently).

Input:  microscopy/all_positions_vapc.csv  (VapC condition, 'kept' filter mask)
        microscopy/all_positions_shx.csv   (SHX/CASP conditions, 'kept' filter mask)
Output: figure_s6_panelA.svg, figure_s6_panelA_preview.png, written next to
        this script.

Run from this directory:
    cd scripts/supplementary_figures
    python figure_s6_panelA.py
Writes figure_s6_panelA.svg next to this script.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Global style
# ------------------------------------------------------------------
fsize = 14
plt.close("all")

MICR_DIR = os.path.join(_REPO, 'microscopy')

# Condition colours (dis-arrest warm, reg-arrest cool)
DIS_COLOR = '#E07B54'
REG_COLOR = 'steelblue'

# ------------------------------------------------------------------
# Load and filter data
# ------------------------------------------------------------------
vapc = pd.read_csv(os.path.join(MICR_DIR, 'all_positions_vapc.csv'))
shx  = pd.read_csv(os.path.join(MICR_DIR, 'all_positions_shx.csv'))

# 'kept' is the upstream segmentation QC mask (per-cell pass/fail); only cells
# that passed are used for the area statistics
vapc_filt = vapc[vapc['kept'] == True].copy()
shx_filt  = shx[shx['kept'] == True].copy()

shx_filt['label'] = shx_filt['condition'].map({'SHX': 'Dis-Arrest\n(SHX)', 'CASP': 'Reg-Arrest'})
vapc_filt['label'] = vapc_filt['condition'].map({'VapC': 'VapC 24h'})

PANEL_A_ORDER = ['Dis-Arrest\n(SHX)', 'VapC 24h', 'Reg-Arrest']
PANEL_A_COLORS = [DIS_COLOR, DIS_COLOR, REG_COLOR]

area_data = {
    'Dis-Arrest\n(SHX)':  shx_filt.loc[shx_filt['label'] == 'Dis-Arrest\n(SHX)',  'area_px'].values,
    'VapC 24h': vapc_filt.loc[vapc_filt['label'] == 'VapC 24h', 'area_px'].values,
    'Reg-Arrest': shx_filt.loc[shx_filt['label'] == 'Reg-Arrest', 'area_px'].values
}


# ------------------------------------------------------------------
# Helper: annotated violin plot with per-violin stats
# ------------------------------------------------------------------
def _violin_with_stats(ax, data_dict, order, colors, ylabel, title):
    """Violin plot of data_dict[k] for each k in order, with a median marker,
    an IQR bar, and a text annotation (n, mean, coefficient of variation)
    placed above each violin."""
    datasets = [data_dict[k] for k in order]

    parts = ax.violinplot(datasets, positions=range(len(order)),
                          showmedians=False, showextrema=False)
    for pc, col in zip(parts['bodies'], colors):
        pc.set_facecolor(col)
        pc.set_edgecolor('0.3')
        pc.set_alpha(0.75)

    for i, vals in enumerate(datasets):
        q25, q75 = np.percentile(vals, [25, 75])
        med = np.median(vals)
        ax.plot([i, i], [q25, q75], color='0.2', lw=2, solid_capstyle='round')
        ax.plot(i, med, 'o', color='white', markeredgecolor='0.2',
                markersize=4, zorder=3)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=fsize - 2)
    ax.set_ylabel(ylabel, fontsize=fsize - 1, labelpad=2)
    ax.tick_params(axis='y', labelsize=fsize - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # extend the y-axis (beyond matplotlib's auto range) to leave headroom for
    # the three-line text annotation above each violin
    ymin, ymax_auto = ax.get_ylim()
    y_pad = (ymax_auto - ymin) * 0.04
    text_top = max(np.max(v) for v in datasets) + (ymax_auto - ymin) * 0.40
    ax.set_ylim(ymin, text_top)

    for i, (key, vals) in enumerate(zip(order, datasets)):
        mean_val = vals.mean()
        cv_val   = vals.std() / mean_val if mean_val != 0 else float('nan')  # coefficient of variation
        ax.text(i, np.max(vals) + y_pad,
                f'n={len(vals)}\nmean={mean_val:.0f}\nCV={cv_val:.2f}',
                ha='center', va='bottom', fontsize=fsize - 4)

    ax.set_title(title, fontsize=fsize, pad=4)


def panel_A(ax):
    _violin_with_stats(ax, area_data, PANEL_A_ORDER, PANEL_A_COLORS,
                       'Cell area (px)', 'Cell area distributions')


# ------------------------------------------------------------------
# Assemble — single narrow panel
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(4.5, 3.2), label_offset=(-0.04, 0.02))
pf.add_panel([0.16, 0.2, 0.8, 0.68], draw_func=panel_A, label=' ')

pf.save("figure_s6_panelA.svg", dpi=300, transparent=True)
pf.fig.savefig("figure_s6_panelA_preview.png", dpi=200)
