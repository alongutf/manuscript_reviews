"""
Supplementary Figure S7 — Microfluidics analysis of the SHX experiment.

Panels:
  A. Kymograph image showing lineage trench over time, with annotation marking
     when SHX was added.
  B. Histogram of cell growth-halt times (frames → minutes), with dashed line
     at frame 18 (SHX addition).
  C. Histogram of division event times (frames → minutes), sharing panel B's
     x-axis.

Both histograms are taken from microscopy/true_events.csv (manually validated
events); halt events with no recorded frame are dropped.

Run from this directory:
    cd scripts/supplementary_figures
    python figure_s7.py
The figure is written next to this script as figure_s7.pdf.
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
import matplotlib.image as mpimg
from matplotlib.patches import FancyArrowPatch

from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Global style
# ------------------------------------------------------------------
fsize = 10
plt.close("all")

MICR_DIR  = os.path.join(_REPO, 'microscopy')
IMG_DIR   = os.path.join(MICR_DIR, 'images')

DIS_COLOR = '#E07B54'
DIV_COLOR = '#4C86A8'

MIN_PER_FRAME = 10        # 1 frame = 10 minutes
SHX_FRAME = 18

BIN_WIDTH_MIN = 1 * MIN_PER_FRAME   # 2 frames per bin

# ------------------------------------------------------------------
# Load data — manually validated events
# ------------------------------------------------------------------
events = pd.read_csv(os.path.join(MICR_DIR, 'true_events.csv'))
events = events.dropna(subset=['frame'])        # drop events with no frame

halt_min = events.loc[events['event_type'] == 'halt', 'frame'].values * MIN_PER_FRAME
div_min  = events.loc[events['event_type'] == 'division', 'frame'].values * MIN_PER_FRAME

# Shared x-axis for both histograms
_max_min = max(halt_min.max(), div_min.max())
BINS = np.arange(0, _max_min + BIN_WIDTH_MIN, BIN_WIDTH_MIN)
XLIM = (0, 500)


# ==================================================================
# Panels
# ==================================================================
def panel_A(ax):
    img = mpimg.imread(os.path.join(IMG_DIR, 'lineage_trench24_timeaxis.png'))
    ax.imshow(img, aspect='auto')
    ax.set_axis_off()

    # Find the x-position of the red dashed line in image pixel coordinates.
    # We add the arrow + label in axes fraction coordinates so it is independent
    # of the exact pixel position: place at the fraction corresponding to the
    # dashed line's location visually (roughly 18/total_frames along the time axis).
    # The image already has the dashed line drawn; we just annotate it.
    ax.annotate(
        'SHX added',
        xy=(0.485, 0.97),          # tip of arrow — near the red dashed line
        xytext=(0.485, 1.06),      # text sits above
        xycoords='axes fraction',
        textcoords='axes fraction',
        ha='center', va='bottom',
        fontsize=fsize - 1,
        color='red',
        arrowprops=dict(arrowstyle='->', color='red', lw=1.2),
    )


def _event_hist(ax, values, color, xlabel, title):
    ax.hist(values, bins=BINS, color=color, edgecolor='0.3',
            linewidth=0.6, alpha=0.85)

    ax.axvline(SHX_FRAME * MIN_PER_FRAME, color='red', linestyle='--',
               linewidth=1.5, label='SHX added')

    ax.set_xlim(*XLIM)
    ax.set_xlabel(xlabel, fontsize=fsize - 1)
    ax.set_ylabel('Number of cells', fontsize=fsize - 1)
    ax.set_title(title, fontsize=fsize)
    ax.tick_params(labelsize=fsize - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=fsize - 2, frameon=False)


def panel_B(ax):
    _event_hist(ax, halt_min, DIS_COLOR,
                'Halt time (minutes)', 'Elongation arrest times')


def panel_C(ax):
    _event_hist(ax, div_min, DIV_COLOR,
                'Division time (minutes)', 'Division times')


# ------------------------------------------------------------------
# Assemble
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 7), label_offset=(-0.02, 0.04))

# Panel A: kymograph image (top, wide)
pf.add_panel([0.05, 0.60, 0.88, 0.36], draw_func=panel_A, label='A')

# Panel B: halt-time histogram (middle) — short and wide
pf.add_panel([0.10, 0.37, 0.85, 0.18], draw_func=panel_B, label='B')

# Panel C: division-time histogram (bottom), matched x-axis with B
pf.add_panel([0.10, 0.06, 0.85, 0.18], draw_func=panel_C, label='C')

pf.save("figure_s7.pdf", dpi=300, transparent=True)
pf.fig.savefig("figure_s7_preview.png", dpi=200)
plt.show()
