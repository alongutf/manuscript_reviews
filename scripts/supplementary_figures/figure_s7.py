"""
Supplementary Figure S7 — Microfluidics analysis of the SHX experiment.

Time in every panel is normalized to the moment SHX is added (frame 18 =
180 min), so t = 0 is the drug-addition time and pre-drug times are negative.

Panels:
  A. Kymograph image showing lineage trench over time, with annotation marking
     when SHX was added. The image's baked-in time axis is cropped away and
     redrawn in SHX-relative time.
  B. Histogram of cell growth-halt times (frames → minutes), with dashed line
     at t = 0 (SHX addition).
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
SHX_MIN = SHX_FRAME * MIN_PER_FRAME   # 180 min — drug addition, defines t = 0

BIN_WIDTH_MIN = 1 * MIN_PER_FRAME   # 2 frames per bin

TIME_LABEL = 'Time from SHX addition (min)'

# Geometry of the baked-in time axis in lineage_trench24_timeaxis.png, measured
# from the tick labels: t = 0 min sits at pixel column 67.5 and one minute spans
# 10.0125 px. Rows 18..1024 / columns 18..3723 are the kymograph box itself;
# everything below is the original (now discarded) axis.
IMG_T0_PX   = 67.5
IMG_PX_PER_MIN = 10.0125
IMG_BOX = dict(row0=18, row1=1025, col0=18, col1=3724)

# ------------------------------------------------------------------
# Load data — manually validated events
# ------------------------------------------------------------------
events = pd.read_csv(os.path.join(MICR_DIR, 'true_events.csv'))
events = events.dropna(subset=['frame'])        # drop events with no frame

# Times are expressed relative to SHX addition (t = 0 at 180 min)
halt_min = events.loc[events['event_type'] == 'halt', 'frame'].values * MIN_PER_FRAME - SHX_MIN
div_min  = events.loc[events['event_type'] == 'division', 'frame'].values * MIN_PER_FRAME - SHX_MIN

# Shared x-axis for both histograms
_max_min = max(halt_min.max(), div_min.max())
BINS = np.arange(-SHX_MIN, _max_min + BIN_WIDTH_MIN, BIN_WIDTH_MIN)
XLIM = (-SHX_MIN, 500 - SHX_MIN)


# ==================================================================
# Panels
# ==================================================================
def panel_A(ax):
    img = mpimg.imread(os.path.join(IMG_DIR, 'lineage_trench24_timeaxis.png'))

    # Crop away the image's own time axis and re-draw it in SHX-relative time.
    box = img[IMG_BOX['row0']:IMG_BOX['row1'], IMG_BOX['col0']:IMG_BOX['col1']].copy()
    left  = (IMG_BOX['col0'] - IMG_T0_PX) / IMG_PX_PER_MIN - SHX_MIN
    right = (IMG_BOX['col1'] - IMG_T0_PX) / IMG_PX_PER_MIN - SHX_MIN

    # The image has a red dashed line baked in at 175 min -- the frame *boundary*
    # before the first post-SHX frame, so it renders half a frame early (t = -5).
    # It sits on the black separator between frames, so erasing it is just a matter
    # of copying a neighbouring separator column over it; the marker is redrawn
    # below at t = 0, centred on the first post-SHX frame itself.
    red = (box[..., 0] > 0.6) & (box[..., 1] < 0.4) & (box[..., 2] < 0.4)
    red_cols = np.where(red.any(axis=0))[0]
    if red_cols.size:
        c0, c1 = red_cols.min(), red_cols.max()
        box[:, c0:c1 + 1] = box[:, c0 - 3][:, None, :]

    ax.imshow(box, aspect='auto', extent=(left, right, 0, 1))
    ax.set_xlim(left, right)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_xlabel(TIME_LABEL, fontsize=fsize - 1)
    ax.tick_params(labelsize=fsize - 2)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)

    # SHX marker, drawn over the middle of the first post-SHX frame (t = 0).
    ax.axvline(0, color='red', linestyle='--', linewidth=1.5)

    ax.annotate(
        'SHX added',
        xy=(0, 1.0),
        xytext=(0, 1.10),
        xycoords=('data', 'axes fraction'),
        textcoords=('data', 'axes fraction'),
        ha='center', va='bottom',
        fontsize=fsize - 1,
        color='red',
        arrowprops=dict(arrowstyle='->', color='red', lw=1.2),
    )


def _event_hist(ax, values, color, xlabel, title):
    ax.hist(values, bins=BINS, color=color, edgecolor='0.3',
            linewidth=0.6, alpha=0.85)

    ax.axvline(0, color='red', linestyle='--',
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
                TIME_LABEL, 'Elongation arrest times')


def panel_C(ax):
    _event_hist(ax, div_min, DIV_COLOR,
                TIME_LABEL, 'Division times')


# ------------------------------------------------------------------
# Assemble
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 7), label_offset=(-0.02, 0.04))

# Panel A: kymograph image (top, wide)
pf.add_panel([0.05, 0.64, 0.88, 0.30], draw_func=panel_A, label='A')

# Panel B: halt-time histogram (middle) — short and wide
pf.add_panel([0.10, 0.37, 0.85, 0.18], draw_func=panel_B, label='B')

# Panel C: division-time histogram (bottom), matched x-axis with B
pf.add_panel([0.10, 0.06, 0.85, 0.18], draw_func=panel_C, label='C')

pf.save("figure_s7.pdf", dpi=300, transparent=True)
pf.fig.savefig("figure_s7_preview.png", dpi=200)
plt.show()
