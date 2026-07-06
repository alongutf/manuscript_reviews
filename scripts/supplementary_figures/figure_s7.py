"""
Supplementary Figure S7 — Microfluidics analysis of the SHX experiment.

Panels:
  A. Kymograph image showing lineage trench over time, with annotation marking
     when SHX was added.
  B. Histogram of cell growth-halt times (frames → minutes), with dashed line
     at frame 18 (SHX addition).

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

FRAMES_PER_MIN = 1 / 10   # 1 frame = 10 minutes
SHX_FRAME = 18

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
halt = pd.read_csv(os.path.join(MICR_DIR, 'halt_times_all.csv'))


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


def panel_B(ax):
    halt_min = halt['halt_frame'].values * 10   # frames → minutes

    bin_width = 2 * 10   # 2 frames in minutes
    bins = np.arange(0, halt_min.max() + bin_width, bin_width)

    ax.hist(halt_min, bins=bins, color=DIS_COLOR, edgecolor='0.3',
            linewidth=0.6, alpha=0.85)

    shx_min = SHX_FRAME * 10
    ax.axvline(shx_min, color='red', linestyle='--', linewidth=1.5,
               label='SHX added')

    ax.set_xlabel('Halt time (minutes)', fontsize=fsize - 1)
    ax.set_ylabel('Number of cells', fontsize=fsize - 1)
    ax.set_title('Distribution of halt times', fontsize=fsize)
    ax.tick_params(labelsize=fsize - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=fsize - 2, frameon=False)


# ------------------------------------------------------------------
# Assemble
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 6), label_offset=(-0.02, 0.04))

# Panel A: kymograph image (top, wide)
pf.add_panel([0.05, 0.50, 0.88, 0.44], draw_func=panel_A, label='A')

# Panel B: halt-time histogram (bottom, centred)
pf.add_panel([0.22, 0.08, 0.56, 0.34], draw_func=panel_B, label='B')

pf.save("figure_s7.pdf", dpi=300, transparent=True)
pf.fig.savefig("figure_s7_preview.png", dpi=200)
plt.show()
