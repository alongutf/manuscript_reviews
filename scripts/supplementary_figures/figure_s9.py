import os
import sys

# --- import bootstrap: this script lives in scripts/supplementary_figures/ ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                       # repo root
sys.path.insert(0, _REPO)                                             # import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))        # figure_functions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Supplementary Figure S9
#   A. Log-log CCDF of lag-time distributions (Reg-Arrest,
#      Reg-Arrest+SHX, Dis-Arrest). Lag times are shifted by t0 (the
#      minimum appearance time of the exponential control) before
#      plotting. Format follows figure1 panel E.
#   B. Summary table of lag-time data across all conditions in this
#      study together with data adapted from Kaplan et al.
# ------------------------------------------------------------------

fsize = 10
REG_COLOR = 'steelblue'  # Reg-Arrest
DIS_COLOR = '#a50f15'    # Dis-Arrest (red), as in figure1
SHX_COLOR = '#6a51a3'    # Reg-Arrest+SHX (purple)
REG_TINT = to_rgba(REG_COLOR, 0.5)
DIS_TINT = to_rgba(DIS_COLOR, 0.5)

SCANLAG_DIR = os.path.join(_REPO, 'scanlag_data', 'CASP+SHX')

# condition -> (lag-time column, OneMinusCDF column, plot color)
CCDF_CONDITIONS = [
    ('Reg-Arrest',     2, 3, REG_COLOR),
    ('Reg-Arrest\n+SHX', 0, 1, SHX_COLOR),
    ('Dis-Arrest',     4, 5, DIS_COLOR),
]


def _load_ccdf_data():
    """Read casp+shx.xlsx and return, per condition, (t0-shifted lag, CCDF)."""
    raw = pd.read_excel(os.path.join(SCANLAG_DIR, 'casp+shx.xlsx'), header=None)
    out = {}
    for name, lag_col, cdf_col, _ in CCDF_CONDITIONS:
        # row 2 holds the t0 string, e.g. 't0=-48'
        t0 = float(str(raw.iloc[2, lag_col]).split('=')[1])
        lag = pd.to_numeric(raw.iloc[4:, lag_col], errors='coerce').dropna().to_numpy()
        ccdf = pd.to_numeric(raw.iloc[4:, cdf_col], errors='coerce').dropna().to_numpy()
        out[name] = (lag - t0, ccdf)
    return out


def panel_A(ax):
    data = _load_ccdf_data()
    for name, lag_col, cdf_col, color in CCDF_CONDITIONS:
        lag, ccdf = data[name]
        order = np.argsort(lag)
        ax.step(lag[order], ccdf[order], where='post', color=color,
                linewidth=1.2, label=name)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Lag time (min)', fontsize=fsize)
    ax.set_ylabel('Fraction of arrested bacteria', fontsize=fsize, labelpad=0)
    ax.set_xlim(0, 2500)
    ax.set_ylim(0.002, 1.5)
    leg = ax.legend(loc='upper right', fontsize=fsize - 2, frameon=True,
                    framealpha=0.6, edgecolor='none')
    leg.get_frame().set_facecolor('white')
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)


def panel_B(ax):
    """Summary table of lag-time statistics. Mean is t0-shifted."""
    s = pd.read_excel(os.path.join(SCANLAG_DIR, 'scanlag_summary.xlsx'), header=None)
    s = s.set_index(0)  # first column = row labels (Study, Strain, ...)

    ax.set_axis_off()
    col_labels = ['Stress', 'Strain', 'Study', 'Lag time\nmean (min)', 'Std (min)']
    rows = []
    cell_colors = []
    for c in s.columns:
        stress = str(s.loc['Stress', c])
        is_reg = 'Reg-Arrest' in stress
        rep = s.loc['Replicate', c]
        if not pd.isna(rep):
            stress = f'{stress} rep{int(rep)}'
        strain = str(s.loc['Strain', c])
        study = str(s.loc['Study', c])
        t0 = float(s.loc['t0', c])
        mean = float(s.loc['mean', c]) - t0  # shift by t0
        std = float(s.loc['std', c])
        rows.append([stress, strain, study, f'{mean:.0f}', f'{std:.0f}'])
        tint = REG_TINT if is_reg else DIS_TINT
        cell_colors.append([tint] * len(col_labels))

    # Size each column to its widest text (header or cell); columns with
    # less text end up narrower. Width = longest line in the column, plus a
    # little padding, normalised so the columns fill the axes width.
    def _col_widths(header, body, pad=2):
        ncol = len(header)
        w = np.array([
            max(max(len(line) for line in str(cell).split('\n'))
                for cell in [header[c]] + [r[c] for r in body]) + pad
            for c in range(ncol)
        ], dtype=float)
        return w / w.sum()

    tbl = ax.table(cellText=rows, colLabels=col_labels,
                   colWidths=_col_widths(col_labels, rows),
                   cellColours=cell_colors, cellLoc='center', loc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fsize - 3)
    tbl.scale(1, 1.4)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor('#bdbdbd')
        if r == 0:
            cell.set_text_props(fontweight='bold')
            cell.set_facecolor('#f0f0f0')


# ------------------------------------------------------------------
# Assemble — single row, two columns
# ------------------------------------------------------------------
plt.close('all')
pf = PanelFigure(figsize=(7, 3.5), label_offset=(-0.03, 0.02))
pf.add_panel([0.08, 0.27, 0.3, 0.63], label='A', draw_func=panel_A)
pf.add_panel([0.43, 0.05, 0.56, 0.85], label='B', draw_func=panel_B)

pf.save('figure_s9.pdf', dpi=300)
pf.save('figure_s9_preview.png', dpi=200)
print('Saved figure_s9.pdf and figure_s9_preview.png')
plt.show()