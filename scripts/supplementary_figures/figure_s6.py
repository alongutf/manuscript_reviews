"""
Supplementary Figure S6 — Microscopy characterisation of Dis-Arrest vs Reg-Arrest cells.

Panels:
  A. Violin plot of cell area distributions for the four conditions
     (SHX⁺, SHX⁻, VapC⁺ 24h, VapC⁻). Mean, n and CV annotated per condition.
  B. Histogram of constitutive mCherry expression in VapC⁺ 24h cells,
     annotated with mean, n and CV.
  D. Single row of representative phase/fluorescence images, one per condition
     in the same order (and with the same labels) as panel A.

Run from this directory:
    cd scripts/supplementary_figures
    python figure_s6.py
The figure is written next to this script as figure_s6.pdf.
"""
import os
import sys

from networkx.algorithms.bipartite.basic import density

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch

from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Global style
# ------------------------------------------------------------------
fsize = 10
plt.close("all")

MICR_DIR = os.path.join(_REPO, 'microscopy')
IMG_DIR  = os.path.join(MICR_DIR, 'images')

# Condition colours (dis-arrest warm, reg-arrest cool)
DIS_COLOR = '#E07B54'
REG_COLOR = 'steelblue'

# ------------------------------------------------------------------
# Load and filter data
# ------------------------------------------------------------------
vapc = pd.read_csv(os.path.join(MICR_DIR, 'all_positions_vapc.csv'))
shx  = pd.read_csv(os.path.join(MICR_DIR, 'all_positions_shx.csv'))

vapc_filt = vapc[vapc['kept'] == True].copy()
shx_filt  = shx[shx['kept'] == True].copy()

# SHX: 'SHX' -> Dis-Arrest, 'CASP' -> Reg-Arrest
shx_filt = shx_filt.copy()
shx_filt['label'] = shx_filt['condition'].map({'SHX': 'Dis-Arrest (SHX)', 'CASP': 'Reg-Arrest (SHX)'})

# VapC: 'VapC' -> Dis-Arrest, 'Reg-Arrest' -> Reg-Arrest
vapc_filt = vapc_filt.copy()
vapc_filt['label'] = vapc_filt['condition'].map({'VapC': 'VapC 24h', 'Reg-Arrest': 'Reg-Arrest (VapC)'})

# Ordered groups for panel A (data keys) and their display labels
PANEL_A_ORDER = ['Dis-Arrest (SHX)', 'Reg-Arrest (SHX)', 'VapC 24h', 'Reg-Arrest (VapC)']
PANEL_A_LABELS = ['SHX$^+$', 'SHX$^-$', 'VapC$^+$ 24h', 'VapC$^-$']
PANEL_A_COLORS = [DIS_COLOR, REG_COLOR, DIS_COLOR, REG_COLOR]

area_data = {
    'Dis-Arrest (SHX)':  shx_filt.loc[shx_filt['label'] == 'Dis-Arrest (SHX)',  'area_px'].values,
    'Reg-Arrest (SHX)':  shx_filt.loc[shx_filt['label'] == 'Reg-Arrest (SHX)',  'area_px'].values,
    'VapC 24h': vapc_filt.loc[vapc_filt['label'] == 'VapC 24h', 'area_px'].values,
    'Reg-Arrest (VapC)': vapc_filt.loc[vapc_filt['label'] == 'Reg-Arrest (VapC)', 'area_px'].values,
}

mcherry_data = {
    'VapC 24h': vapc_filt.loc[vapc_filt['label'] == 'VapC 24h', 'mcherry_bgsub_median'].values,
    'Reg-Arrest (VapC)': vapc_filt.loc[vapc_filt['label'] == 'Reg-Arrest (VapC)', 'mcherry_bgsub_median'].values,
}

yfp_data = {
    'Dis-Arrest (SHX)': shx_filt.loc[shx_filt['label'] == 'Dis-Arrest (SHX)', 'yfp_bgsub_median'].values,
    'Reg-Arrest (SHX)': shx_filt.loc[shx_filt['label'] == 'Reg-Arrest (SHX)', 'yfp_bgsub_median'].values,
}


# ------------------------------------------------------------------
# Helper: annotated violin plot
# ------------------------------------------------------------------
def _violin(ax, data_dict, order, colors, ylabel, title, fsize):
    datasets = [data_dict[k] for k in order]

    parts = ax.violinplot(datasets, positions=range(len(order)),
                          showmedians=False, showextrema=False)
    for pc, col in zip(parts['bodies'], colors):
        pc.set_facecolor(col)
        pc.set_edgecolor('0.3')
        pc.set_alpha(0.75)

    for i, (key, vals) in enumerate(zip(order, datasets)):
        q25, q75 = np.percentile(vals, [25, 75])
        med = np.median(vals)
        ax.plot([i, i], [q25, q75], color='0.2', lw=2, solid_capstyle='round')
        ax.plot(i, med, 'o', color='white', markeredgecolor='0.2',
                markersize=4, zorder=3)

        mean_val = vals.mean()
        cv_val   = vals.std() / mean_val if mean_val != 0 else float('nan')
        ax.text(i, ax.get_ylim()[1] if ax.get_ylim()[1] != 0 else 1,
                f'μ={mean_val:.0f}\nCV={cv_val:.2f}',
                ha='center', va='bottom', fontsize=fsize - 4,
                transform=ax.get_xaxis_transform())

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=fsize - 2)
    ax.set_ylabel(ylabel, fontsize=fsize - 1, labelpad=2)
    ax.set_title(title, fontsize=fsize)
    ax.tick_params(axis='y', labelsize=fsize - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _violin_with_stats(ax, data_dict, order, colors, ylabel, title, xlabels=None):
    datasets = [data_dict[k] for k in order]

    parts = ax.violinplot(datasets, positions=range(len(order)),
                          showmedians=False, showextrema=False)
    for pc, col in zip(parts['bodies'], colors):
        pc.set_facecolor(col)
        pc.set_edgecolor('0.3')
        pc.set_alpha(0.75)

    # IQR bars + median dot
    for i, vals in enumerate(datasets):
        q25, q75 = np.percentile(vals, [25, 75])
        med = np.median(vals)
        ax.plot([i, i], [q25, q75], color='0.2', lw=2, solid_capstyle='round')
        ax.plot(i, med, 'o', color='white', markeredgecolor='0.2',
                markersize=4, zorder=3)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(xlabels if xlabels is not None else order, fontsize=fsize - 2)
    ax.set_ylabel(ylabel, fontsize=fsize - 1, labelpad=2)
    ax.tick_params(axis='y', labelsize=fsize - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Place stats just above each violin's own maximum
    ymin, ymax_auto = ax.get_ylim()
    y_pad = (ymax_auto - ymin) * 0.04   # small gap above violin top

    # Expand ylim to accommodate two lines of annotation text (~25% extra)
    text_top = max(np.max(v) for v in datasets) + (ymax_auto - ymin) * 0.40
    ax.set_ylim(ymin, text_top)

    for i, (key, vals) in enumerate(zip(order, datasets)):
        mean_val = vals.mean()
        cv_val   = vals.std() / mean_val if mean_val != 0 else float('nan')
        ax.text(i, np.max(vals) + y_pad,
                f'n={len(vals)}\nmean={mean_val:.0f}\nCV={cv_val:.2f}',
                ha='center', va='bottom', fontsize=fsize - 4)

    ax.set_title(title, fontsize=fsize, pad=4)


# ==================================================================
# Panels
# ==================================================================
def panel_A(ax):
    _violin_with_stats(ax, area_data, PANEL_A_ORDER, PANEL_A_COLORS,
                       'Cell area (px)', 'Cell area distributions',
                       xlabels=PANEL_A_LABELS)


def panel_B(ax):
    """Histogram of constitutive mCherry expression in VapC⁺ 24h cells."""
    vals = mcherry_data['VapC 24h']
    ax.hist(vals, bins=30, color=DIS_COLOR, alpha=0.8, edgecolor='0.3', linewidth=0.4, density=True)

    mean_val = vals.mean()
    cv_val   = vals.std() / mean_val if mean_val != 0 else float('nan')
    ax.text(0.95, 0.95,
            f'n={len(vals)}\nmean={mean_val:.0f}\nCV={cv_val:.2f}',
            transform=ax.transAxes, ha='right', va='top', fontsize=fsize - 3,
            bbox=dict(boxstyle='round,pad=0.35', facecolor='none', edgecolor='0.6',
                      linewidth=0.6))

    ax.set_xlabel('mCherry (arbitrary units)', fontsize=fsize - 1, labelpad=2)
    ax.set_ylabel('Density', fontsize=fsize - 1, labelpad=2)
    ax.set_title('VapC expression', fontsize=fsize, pad=2)
    ax.tick_params(labelsize=fsize - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_D(ax):
    """Single row of 4 microscopy images in the same order as panel A."""
    images = [
        ('shx2.png',      'SHX$^+$',       DIS_COLOR),
        ('shx_reg1.png',  'SHX$^-$',       REG_COLOR),
        ('vapc1.png',     'VapC$^+$ 24h',  DIS_COLOR),
        ('vapc_reg1.png', 'VapC$^-$',      REG_COLOR),
    ]

    n = len(images)
    w = 0.22
    gap = 0.025
    margin = (1.0 - n * w - (n - 1) * gap) / 2
    y0 = 0.13
    h = 0.80

    ax.set_axis_off()
    for i, (fname, label, lc) in enumerate(images):
        img_path = os.path.join(IMG_DIR, fname)
        img = mpimg.imread(img_path)
        if img.ndim == 3:
            img_gray = np.dot(img[..., :3], [0.2989, 0.5870, 0.1140])
        else:
            img_gray = img
        x0 = margin + i * (w + gap)
        sub = ax.inset_axes([x0, y0, w, h])
        sub.imshow(img_gray, cmap='gray', aspect='equal')
        sub.set_xticks([]); sub.set_yticks([])
        for spine in sub.spines.values():
            spine.set_edgecolor(lc)
            spine.set_linewidth(2)
        ax.text(x0 + w / 2, y0 - 0.04, label,
                transform=ax.transAxes,
                ha='center', va='top', fontsize=fsize - 2, color=lc)


# ------------------------------------------------------------------
# Assemble
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 5), label_offset=(-0.04, 0.04))

# Row 1: panel A (narrower, area violins) + panel B (mCherry histogram)
pf.add_panel([0.08, 0.62, 0.52, 0.30], draw_func=panel_A, label='A')
pf.add_panel([0.70, 0.62, 0.26, 0.30], draw_func=panel_B, label='B')

# Row 2: panel D — microscopy images (single row)
pf.add_panel([0.1, 0.1, 0.8, 0.4], draw_func=panel_D, hide_axis=True, label='C')

pf.save("figure_s6.pdf", dpi=300, transparent=True)
pf.fig.savefig("figure_s6_preview.png", dpi=200)
plt.show()