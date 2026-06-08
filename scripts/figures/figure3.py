import src.analysis_functions as af
import src.data_functions as df
import matplotlib.pyplot as plt
from figure_functions import PanelFigure
import numpy as np
import pandas as pd
from scipy import stats
import os
import importlib

importlib.reload(af)
importlib.reload(df)
# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------
fsize = 10
plt.close("all")
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
ev_data_dir = os.path.join(root_dir, 'ev_data')
RESULTS_DIR = os.path.join(root_dir, 'results', 'simulation_results')


def _plot_ccdf(ax, npy_file, title):
    arr = np.load(os.path.join(ev_data_dir, npy_file))
    data1 = arr[0, :];  data1 = data1[data1 > 0]
    data2 = arr[1, :];  data2 = data2[data2 > 0]
    x2 = float(np.max(data2))
    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    noise = d1s < x2
    ax.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='noise', markersize=3)
    ax.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
              color='skyblue', label='signal', markersize=3)
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled', markersize=3)
    ax.set_xlim([0.1, np.max(d1s) * 1.5])
    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)
    ax.set_title(title, fontsize=fsize)
    ax.legend(fontsize=fsize - 2)
    ax.tick_params(labelsize=fsize - 2)


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
    _plot_ccdf(ax, 'sample_15b_filtered.npy', 'Reg-Arrest')


def panel_B(ax):
    _plot_ccdf(ax, 'sample_15a_filtered.npy', 'Dis-Arrest')


def panel_C(ax):
    _plot_ccdf(ax, 'simulated_pcs_0.9.npy', r'Simulation ($\rho=0.9$)')


def panel_D(ax):
    _plot_ccdf(ax, 'simulated_pcs_0.5.npy', r'Simulation ($\rho=0.5$)')


def panel_E(ax):
    summary = pd.read_csv(os.path.join(RESULTS_DIR, 'raw', 'rho_sweep_summary.txt'),
                          sep=r'\s+', comment='#', index_col=0)
    rho_vals = summary.index.values
    medians  = summary['median'].values
    stds     = summary['std'].values

    ax.errorbar(rho_vals, medians, yerr=stds, fmt='o-', color='steelblue',
                capsize=4, linewidth=1.8, markersize=6, label='median ± SD')
    ax.fill_between(rho_vals, medians - stds, medians + stds, alpha=0.3, color='steelblue')
    ax.set_xlabel(r'Correlation strength ($\rho$)', fontsize=fsize - 2)
    ax.set_ylabel('GMP-Cor', fontsize=fsize - 2)
    ax.set_title('GMP-Cor calibration curve', fontsize=fsize)
    ax.set_xticks(np.arange(0, 1.05, 0.2))
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(bottom=-1, top=75)
    ax.grid(True, linestyle='--', alpha=0.2)
    #ax.legend(fontsize=fsize - 2)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)


def panel_F(ax):
    path = os.path.join(root_dir, 'results', 'data_metrics', 'test8.csv')
    data = pd.read_csv(path, index_col=0)
    ranking_param = 'sum_denoised_ev'
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

    ax.set_xticklabels(['Regulated', 'Dis-Arrest'], fontsize=fsize - 2, rotation=0, ha='center')
    u_stat, u_p = stats.mannwhitneyu(group1, group0)
    asterisks = format_p(u_p)
    x1, x2 = 1, 2
    y, h = max(max(group1), max(group0)) + 1, 1
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.75, color='black')
    ax.text((x1 + x2) * 0.5, y + h, asterisks, ha='center', va='bottom', color='black', fontsize=fsize - 2)
    ax.set_ylabel('GMP-Cor', fontsize=fsize - 2)
    ax.grid(False)
    ax.set_ylim(bottom=-1, top=40)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)


# Build figure 3:
# Layout: A/C stacked left, B/D stacked middle, E/F stacked right
pf = PanelFigure(figsize=(7, 5.5), label_offset=(-0.04, 0.055))
panel_pos = [
    [0.08, 0.57, 0.23, 0.35],  # A - Reg-Arrest (top-left)
    [0.39, 0.57, 0.23, 0.35],  # B - Dis-Arrest (top-middle)
    [0.08, 0.10, 0.23, 0.35],  # C - Sim ρ=0.9  (bottom-left, below A)
    [0.39, 0.10, 0.23, 0.35],  # D - Sim ρ=0.0  (bottom-middle, below B)
    [0.71, 0.57, 0.27, 0.35],  # E - Calibration curve (top-right)
    [0.71, 0.10, 0.27, 0.35],  # F - Box plot (bottom-right)
]
pf.add_panel(panel_pos[0], draw_func=panel_A)
pf.add_panel(panel_pos[1], draw_func=panel_B)
pf.add_panel(panel_pos[2], draw_func=panel_C)
pf.add_panel(panel_pos[3], draw_func=panel_D)
pf.add_panel(panel_pos[4], draw_func=panel_E)
pf.add_panel(panel_pos[5], draw_func=panel_F)
pf.save("figure3.svg", dpi=300, transparent=True)
plt.show()