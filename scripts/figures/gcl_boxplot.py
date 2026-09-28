"""GCL by category (Regulated vs Dis-Arrest), styled like figure 3E (GMP-Cor)."""
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
import pandas as pd
from scipy import stats
import os

fsize = 10
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
REG_COLOR = 'steelblue'
DIS_COLOR = '#E07B54'
REF_LW = 1
REF_ALPHA = 0.85


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


def panel_gcl(ax):
    path = os.path.join(root_dir, 'results', 'data_metrics', 'data_metricsGCL.csv')
    data = pd.read_csv(path, index_col=0)
    ranking_param = 'GCL'
    to_exclude = ['adam_matrix_filtered.csv', 'deb_Ec_CDS_untreated.csv', 'deb_KP_CDS_untreated.csv']
    data = data[~data['file_name'].isin(to_exclude)].dropna(subset=[ranking_param])
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

    ax.axhline(group1.median(), color=REG_COLOR, linestyle='--', linewidth=REF_LW, alpha=REF_ALPHA)
    ax.axhline(group0.median(), color=DIS_COLOR, linestyle='--', linewidth=REF_LW, alpha=REF_ALPHA)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Regulated', 'Dis-Arrest'], fontsize=fsize - 2, rotation=0, ha='center')
    u_stat, u_p = stats.mannwhitneyu(group1, group0)
    print(f'GCL: n_r={len(group1)}, n_d={len(group0)}, Mann-Whitney U={u_stat}, p={u_p:.4f}')
    x1, x2 = 1, 2
    y_top = max(group1.max(), group0.max())
    h = 0.03
    y = y_top + h * 0.5
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.75, color='black')
    ax.text((x1 + x2) * 0.5, y + h * 1.1, format_p(u_p), ha='center', va='bottom',
            color='black', fontsize=fsize - 2)
    ax.set_ylabel('GCL', fontsize=fsize - 2)
    ax.set_ylim(bottom=-0.02, top=0.8)
    ax.grid(False)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)


fig, ax = plt.subplots(figsize=(2.0, 2.8))
panel_gcl(ax)
fig.tight_layout()
fig.savefig('gcl_boxplot.svg')
fig.savefig('gcl_boxplot_preview.png', dpi=300)
