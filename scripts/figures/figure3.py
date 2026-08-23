import src.analysis_functions as af
import src.data_functions as df
import matplotlib.pyplot as plt
from figure_functions import PanelFigure
import permutation_pvalues as pv
from matplotlib.colors import to_rgba
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

REG_COLOR = 'steelblue'
DIS_COLOR = '#E07B54'

# Shared style for the median reference lines in panels E and F
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


def _plot_ccdf(ax, npy_file, title, signal_color='skyblue'):
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
    ax.set_title(title, fontsize=fsize-2, pad=0)
    # p-value rides in the legend box, under the three series keys; a no-op for
    # the simulated spectra, which have no permutation test
    pv.legend_with_p(ax, npy_file, fontsize=fsize - 2)
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


# ── Row 1: Experimental data (Reg-Arrest / Dis-Arrest) ─────────────
def panel_A(ax):
    _plot_ccdf(ax, 'sample_15b_filtered.npy', 'Reg-Arrest', signal_color=REG_COLOR)


def panel_B(ax):
    _plot_ccdf(ax, 'sample_15a_filtered.npy', 'Dis-Arrest', signal_color=DIS_COLOR)


# ── Row 2: Simulations (high / low ρ) ──────────────────────────────
def panel_C(ax):
    _plot_ccdf(ax, 'simulated_pcs_0.9.npy', r'Simulation ($\chi=0.9$)', signal_color=REG_COLOR)


def panel_D(ax):
    _plot_ccdf(ax, 'simulated_pcs_0.5.npy', r'Simulation ($\chi=0.5$)', signal_color=DIS_COLOR)


# ── Panel E: Box plot (Regulated vs Dis-Arrest GMP-Cor) ────────────
def panel_E(ax):
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
    ax.set_ylabel('GMP-Cor', fontsize=fsize - 2)
    ax.set_ylim(bottom=-1, top=45)
    ax.grid(False)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)


# ------------------------------------------------------------------
# Assemble — 2×2 grid of CCDFs (A–D) with the box plot (E) to the right
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 3.5), label_offset=(-0.04, 0.05))
panel_pos = [
    # Left block — 2×2 experimental / simulation CCDFs (short)
    [0.07, 0.61, 0.26, 0.32],  # A — Reg-Arrest CCDF
    [0.42, 0.61, 0.26, 0.32],  # B — Dis-Arrest CCDF
    [0.07, 0.11, 0.26, 0.32],  # C — Sim χ=0.9 CCDF
    [0.42, 0.11, 0.26, 0.32],  # D — Sim χ=0.5 CCDF
    # Right — box plot spanning the height of the A–D block
    [0.76, 0.25, 0.2, 0.5],  # E — Box plot (Regulated vs Dis-Arrest)
]
pf.add_panel(panel_pos[0], draw_func=panel_A)
pf.add_panel(panel_pos[1], draw_func=panel_B)
pf.add_panel(panel_pos[2], draw_func=panel_C)
pf.add_panel(panel_pos[3], draw_func=panel_D)
pf.add_panel(panel_pos[4], draw_func=panel_E)

pf.save("figure3.pdf", dpi=300, transparent=True)
pf.save("figure3_preview.png", dpi=200)
plt.show()
