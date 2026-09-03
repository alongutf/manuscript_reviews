"""Export the six panels of Figure 3 as individual PNGs, for slide/poster use.

Figure 3 makes the case for GMP-Cor as a correlation-strength metric: panels A-D
compare the empirical eigenvalue CCDF (complementary CDF) against its scrambled
null for two real samples (regulated vs. dis-arrest) and two synthetic ones at
different simulated correlation strengths; panel E is the GMP-Cor calibration
curve against simulated correlation strength chi, with the two real-data group
medians overlaid; panel F is the regulated-vs-dis-arrest GMP-Cor boxplot with a
Mann-Whitney significance bracket.

This is a standalone export utility, not a reusable PanelFigure-based figure
script (contrast scripts/figures/figureN.py, which build one combined SVG via
figure_functions.PanelFigure) -- it saves each panel as its own PNG file to a
presentation folder rather than assembling a single multi-panel manuscript figure.

Reads:
    ev_data/*.npy                              eigenvalue spectra (panels A-D)
    results/simulation_results/raw/rho_sweep_summary.txt   panel E's sweep data
    results/data_metrics/test8.csv             panels E/F's per-sample GMP-Cor

Writes:
    <OUTPUT_DIR>/figure3_panel{A..F}.png        one PNG per panel

Run from anywhere (sys.path is patched to find src/ from this file's location):
    python scripts/figures/export_figure3_panels.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import src.analysis_functions as af
import src.data_functions as df
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
import pandas as pd
from scipy import stats
import importlib

importlib.reload(af)
importlib.reload(df)

fsize = 10
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ev_data_dir = os.path.join(root_dir, 'ev_data')
RESULTS_DIR = os.path.join(root_dir, 'results', 'simulation_results')
OUTPUT_DIR = "/Users/along/Documents/Alon/PhD/documents/IPS 2026"

REG_COLOR = 'steelblue'
DIS_COLOR = '#E07B54'
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


def _plot_ccdf(ax, npy_file, title, signal_color='skyblue'):
    arr = np.load(os.path.join(ev_data_dir, npy_file))
    data1 = arr[0, :]; data1 = data1[data1 > 0]
    data2 = arr[1, :]; data2 = data2[data2 > 0]
    x2 = float(np.max(data2))
    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    noise = d1s < x2
    ax.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='spurious', markersize=3)
    ax.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
              color=signal_color, label='signal', markersize=3)
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled', markersize=3)
    ax.set_xlim([0.1, 30])
    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
    ax.text(x2 * 1.1, 0.8, r'$\lambda_\mathrm{max}^\mathrm{scr}$',
            fontsize=fsize - 2, va='center', ha='left', color='k', alpha=0.7,
            transform=ax.get_xaxis_transform())
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)
    ax.set_title(title, fontsize=fsize)
    ax.legend(fontsize=fsize - 2)
    ax.tick_params(labelsize=fsize - 2)


def panel_A(ax):
    _plot_ccdf(ax, 'sample_15b_filtered.npy', 'Reg-Arrest', signal_color=REG_COLOR)


def panel_B(ax):
    _plot_ccdf(ax, 'sample_15a_filtered.npy', 'Dis-Arrest', signal_color=DIS_COLOR)


def panel_C(ax):
    _plot_ccdf(ax, 'simulated_pcs_0.9.npy', r'Simulation ($\chi=0.9$)', signal_color=REG_COLOR)


def panel_D(ax):
    _plot_ccdf(ax, 'simulated_pcs_0.5.npy', r'Simulation ($\chi=0.5$)', signal_color=DIS_COLOR)


def panel_E(ax):
    summary = pd.read_csv(os.path.join(RESULTS_DIR, 'raw', 'rho_sweep_summary.txt'),
                          sep=r'\s+', comment='#', index_col=0)
    rho_vals = summary.index.values
    medians  = summary['median'].values
    stds     = summary['std'].values

    ax.errorbar(rho_vals, medians, yerr=stds, fmt='o-', color='steelblue',
                capsize=4, linewidth=1.8, markersize=6)
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


panels = [
    ('A', panel_A, (3.5, 2.8)),
    ('B', panel_B, (3.5, 2.8)),
    ('C', panel_C, (3.5, 2.8)),
    ('D', panel_D, (3.5, 2.8)),
    ('E', panel_E, (4.0, 3.0)),
    ('F', panel_F, (3.0, 3.0)),
]

for label, func, figsize in panels:
    fig, ax = plt.subplots(figsize=figsize)
    func(ax)
    fig.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'figure3_panel{label}.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved panel {label} → {out_path}")
