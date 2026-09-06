import os
import sys

# --- import bootstrap: this script lives in scripts/supplementary_figures/ ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                       # repo root
sys.path.insert(0, _REPO)                                             # import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))        # figure_functions

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Supplementary Figure S10
#   A. GMP-Cor calibration curve — GMP-Cor vs. correlation strength (chi),
#      with Reg-Arrest / Dis-Arrest experimental medians as references.
#   B. GMP-Cor box plot — Regulated vs. Dis-Arrest datasets (shares the
#      y-axis with A). Reproduces the former figure3 panels E-F side by side.
# ------------------------------------------------------------------

fsize = 10
root_dir = _REPO
RESULTS_DIR = os.path.join(root_dir, 'results', 'simulation_results')
MODEL_DIR = os.path.join(root_dir,'model fit','gmp_sweep_alpha2')

REG_COLOR = 'steelblue'
DIS_COLOR = '#E07B54'

# Shared style for the median reference lines in panels A and B
# (kept identical so the dashed lines have the same width/opacity in both).
REF_LW = 1
REF_ALPHA = 0.85


def _load_group_medians():
    path = os.path.join(root_dir, 'results', 'data_metrics', 'test8.csv')
    data = pd.read_csv(path, index_col=0)
    ranking_param = 'sum_denoised_ev'
    g1 = data[data['category'] == 'r'][ranking_param]
    g0 = data[data['category'] == 'd'][ranking_param]
    return g1.median(), g0.median()


med_reg, med_dis = _load_group_medians()


def format_p(p):
    if p < 0.0001:
        return '****'
    elif p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'NS'


def panel_A(ax):
    summary = pd.read_csv(os.path.join(RESULTS_DIR, 'raw', 'rho_sweep_summary.txt'),
                          sep=r'\s+', comment='#', index_col=0)
    rho_vals = summary.index.values
    medians  = summary['median'].values
    stds     = summary['std'].values
    analytics = pd.read_csv(os.path.join(MODEL_DIR,'analytical_curve.csv'), index_col=0, header=0)
    ax.plot(analytics['chi'],analytics['gmp_cor'],linestyle='--', linewidth=1.5, color='k', label='analytical')
    ax.errorbar(rho_vals, medians, yerr=stds, fmt='o-', color='steelblue',
                capsize=2, linewidth=1.5, markersize=3, label='simulated')
    ax.fill_between(rho_vals, medians - stds, medians + stds, alpha=0.3, color='steelblue')
    ax.axhline(med_reg, color=REG_COLOR, linestyle='--', linewidth=REF_LW, alpha=REF_ALPHA,
               label='Reg-Arrest median')
    ax.axhline(med_dis, color=DIS_COLOR, linestyle='--', linewidth=REF_LW, alpha=REF_ALPHA,
               label='Dis-Arrest median')
    ax.set_xlabel(r'Correlation strength ($\chi$)', fontsize=fsize - 2)
    ax.set_ylabel('GMP-Cor', fontsize=fsize - 2)
    ax.set_title('GMP-Cor calibration curve', fontsize=fsize)
    ax.set_xticks(np.arange(0, 1.05, 0.2))
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(bottom=-1, top=45)
    ax.grid(True, linestyle='--', alpha=0.2)
    ax.legend(fontsize=fsize - 2)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)


def panel_B(ax):
    path = os.path.join(root_dir, 'results', 'data_metrics', 'data_metrics.csv')
    data = pd.read_csv(path, index_col=0)
    ranking_param = 'sum_denoised_ev'
    to_exclude = ['adam_matrix_filtered.csv','deb_Ec_CDS_untreated.csv','deb_KP_CDS_untreated.csv']
    data = data[~data['file_name'].isin(to_exclude)]
    data['Rank'] = data[ranking_param].rank(method='min').astype(int)
    group1 = data[data['category'] == 'r'][ranking_param]
    group0 = data[data['category'] == 'd'][ranking_param]
    c = "k"
    box = ax.boxplot([group1, group0], meanline=True, showmeans=True, patch_artist=True,
                     boxprops=dict(facecolor="None", color=c), whiskerprops=dict(color=c),
                     capprops=dict(color=c),
                     flierprops=dict(markeredgecolor=c, markersize=2), medianprops=dict(color=c))
    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        for item in box[element]:
            item.set_linewidth(1)
    for mean_line in box['means']:
        mean_line.set_linewidth(1)
        mean_line.set_color(c)
        mean_line.set_linestyle('solid')

    for patch, color in zip(box['boxes'], [REG_COLOR, DIS_COLOR]):
        patch.set_edgecolor(color)
        patch.set_facecolor((*to_rgba(color)[:3], 0.12))

    rng = np.random.default_rng(42)
    for i, (grp, color) in enumerate([(group1, REG_COLOR), (group0, DIS_COLOR)], start=1):
        jitter = rng.uniform(-0.12, 0.12, size=len(grp))
        ax.scatter(i + jitter, grp, color=color, alpha=0.5, s=10, zorder=3, edgecolors='none')

    ax.axhline(med_reg, color=REG_COLOR, linestyle='--', linewidth=REF_LW, alpha=REF_ALPHA)
    ax.axhline(med_dis, color=DIS_COLOR, linestyle='--', linewidth=REF_LW, alpha=REF_ALPHA)

    ax.set_xticklabels(['Regulated', 'Dis-Arrest'], fontsize=fsize - 2, rotation=0, ha='center')
    u_stat, u_p = stats.mannwhitneyu(group1, group0)
    asterisks = format_p(u_p)
    x1, x2 = 1, 2
    y_top = max(group1.max(), group0.max())
    h = 1
    y = y_top + h * 0.5
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.75, color='black')
    ax.text((x1 + x2) * 0.5, y + h * 1.1, asterisks, ha='center', va='bottom',
            color='black', fontsize=fsize - 2)
    ax.set_ylabel('')
    ax.set_ylim(bottom=-1, top=45)
    ax.grid(False)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)


# ------------------------------------------------------------------
# Assemble — single row, two columns (side-by-side, shared y-axis)
# ------------------------------------------------------------------
plt.close('all')
pf = PanelFigure(figsize=(7, 3.2), label_offset=(-0.04, 0.02))
ax_A = pf.add_panel([0.15, 0.18, 0.33, 0.70], label='A', draw_func=panel_A)
ax_B = pf.add_panel([0.55, 0.18, 0.30, 0.70], label='B', draw_func=panel_B)

# Share y-axis between A and B so reference lines align exactly
ax_B.sharey(ax_A)

pf.save('figure_s10.pdf', dpi=300)
pf.save('figure_s10_preview.png', dpi=200)
print('Saved figure_s10.pdf and figure_s10_preview.png')
plt.show()
