"""Standalone lag-time histograms for the kaplan_shx dataset, one-off export for a
conference poster/talk (IPS 2026), not part of the main figureN.py figure set.

Reads the two single-column lag-time text files in scanlag_data/kaplan_shx/ (t0:
start of SHX exposure, t1346: after 1346 min) and writes a small histogram PNG for
each, sharing bins/axes so the two can be placed side by side as "before" and
"after" panels.

Input:  scanlag_data/kaplan_shx/t0.txt, scanlag_data/kaplan_shx/t1346.txt
        (each a whitespace/newline-separated list of per-cell lag times, one value
        per cell; loaded with np.loadtxt as a 1-D array)
Output: <output_dir>/panel_B_t0.png, <output_dir>/panel_B_t1346.png
        NOTE: output_dir is a hard-coded local path on the original author's
        machine (see FINDINGS in the accompanying log) -- unlike the figureN.py
        scripts, this one does not write into scripts/figures/.

Usage:
    python export_panel_B_histograms.py
"""
import numpy as np
import matplotlib.pyplot as plt
import os

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
# hard-coded, machine-specific export path (see module docstring / findings log)
output_dir = "/Users/along/Documents/Alon/PhD/documents/IPS 2026"

fsize = 12
bins = np.linspace(0, 1000, 40)

hist_data = np.loadtxt(os.path.join(project_dir, 'scanlag_data', 'kaplan_shx', 't0.txt'))
hist_data2 = np.loadtxt(os.path.join(project_dir, 'scanlag_data', 'kaplan_shx', 't1346.txt'))

# zero the lag-time axis to the first entry of the t0 sample, so both histograms
# (t0 and t1346) share the same origin
x0 = hist_data[0]
hist_data = hist_data - x0
hist_data2 = hist_data2 - x0

def make_hist_fig(data, color, xlabel=False):
    """Small standalone histogram of lag times, styled to match panel B's look
    (no y-axis ticks/labels, top/right spines removed)."""
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.hist(data, bins=bins, color=color, alpha=0.7)
    ax.set_ylabel('# of cells', fontsize=fsize, labelpad=0)
    ax.set_xlim(0, 1000)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    if xlabel:
        ax.set_xlabel('Lag time (min)', fontsize=fsize, labelpad=0)
    else:
        ax.set_xticks([])
    fig.tight_layout()
    return fig

fig1 = make_hist_fig(hist_data, color='steelblue', xlabel=True)
fig1.savefig(os.path.join(output_dir, 'panel_B_t0.png'), dpi=300, bbox_inches='tight')
plt.close(fig1)

fig2 = make_hist_fig(hist_data2, color='#CD5C5C', xlabel=True)
fig2.savefig(os.path.join(output_dir, 'panel_B_t1346.png'), dpi=300, bbox_inches='tight')
plt.close(fig2)

print("Saved panel_B_t0.png and panel_B_t1346.png to", output_dir)