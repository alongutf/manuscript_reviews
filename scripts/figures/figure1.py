import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from scipy.stats import t
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import os
from collections import OrderedDict
from figure_functions import *

# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------

def panel_C(ax):
    # Kill curve panel (matches kill curves/Kill_curve.svg)
    # Load data
    project_dir = os.path.dirname(os.path.dirname(root_dir))
    df = pd.read_csv(os.path.join(project_dir, 'kill curves', 'Kill_curve.csv'))
    # Colors matching the reference SVG
    colors = {'Reg-Arrest': '#0F8554', 'Dis-Arrest': '#CC503E'}
    labels = {'Reg-Arrest': 'Reg-Arrest', 'Dis-Arrest': 'Dis-Arrest'}
    for curve in ['Reg-Arrest', 'Dis-Arrest']:
        sub = df[df['curve'] == curve].sort_values('time_h')
        x = sub['time_h'].values
        y = sub['mean_survival_fraction'].values
        # convert absolute low/high bounds to asymmetric error deltas
        yerr = np.vstack([y - sub['error_low'].values, sub['error_high'].values - y])
        ax.errorbar(x, y, yerr=yerr, fmt='o-', markersize=4, capsize=2,
                    linewidth=1, color=colors[curve], label=labels[curve])
    # Axis and formatting
    ax.set_yscale('log')
    ax.set_xlabel('Time in Ampicillin (h)', fontsize=fsize-1, labelpad=1)
    ax.set_ylabel('Survival fraction', fontsize=fsize-1, labelpad=1)
    ax.set_title('Rotem et. al', style='italic', fontsize=fsize-2, pad=2)
    ax.tick_params(axis='both', labelsize=fsize-2)
    ax.set_xlim(0, 25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=fsize-3, frameon=False, loc='upper right')


def panel_A(axes):
    # Define the sigmoidal function
    def sigmoid(x, a, b, c):
        return c / (1 + a * np.exp(-b * x))

    # Define the inverse gamma distribution
    def inv_gamma(x, a, b):
        return b ** a / sc.special.gamma(a) * x ** (-a - 1) * np.exp(-b / x)

    # Define the gamma distribution
    def gamma(x, a, b):
        return b ** a / sc.special.gamma(a) * x ** (a - 1) * np.exp(-b * x)

    t0 = 0.8
    lw = 1
    green = '#31a354'
    red = '#CD5C5C'
    # Generate x values
    x = np.linspace(0.01, 5, 100)
    # Plot the regulated sigmoid function
    axes[0, 0].plot(x, sigmoid(x, 10, 4, 0.7), color=green, linewidth=lw+1)
    # Plot the disrupted sigmoid function
    axes[1, 0].plot(x[x < t0], sigmoid(x[x < t0], 10, 4, 1), color=red, linewidth=lw+1)
    axes[1, 0].plot(x[x >= t0], sigmoid(x[x >= t0], 10, 4, 1), color=red, linewidth=lw, linestyle='dashed')
    axes[1, 0].plot([t0, 4], [sigmoid(t0, 10, 4, 1), sigmoid(t0, 10, 4, 1)], color=red, linewidth=lw+1)
    axes[1, 0].plot([t0, t0], [0.55, 0.85], color='k', linewidth=lw)
    # add text to the plot
    axes[1, 0].text(0.85, 0.4, 'Acute\nstress', fontsize=fsize-2, color='k')
    # remove the top and right spines
    for ax in axes[:, 0]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # make the spines wider
        ax.spines['left'].set_linewidth(lw)
        ax.spines['bottom'].set_linewidth(lw)
        ax.set_xticks([])  # remove ticks
        ax.set_yticks([])
        ax.set_xlabel('Time', fontsize=fsize, loc='right')
        ax.set_ylabel('# of cells', fontsize=fsize, labelpad=0)
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 1.1)

    # second column:
    # Plot the regulated sigmoid function
    axes[0, 1].plot(x, 0.67 + sigmoid(x - 0.2, 10, 4, 0.7), color=green, linewidth=lw+1)
    # Plot the disrupted sigmoid function
    axes[1, 1].plot(x, 0.705 + sigmoid(x - 0.65, 10, 4, 0.7), color=red, linewidth=lw+1)

    for ax in axes[:, 1]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        # make the spines wider
        ax.spines['left'].set_linewidth(lw)
        ax.spines['bottom'].set_linewidth(lw)
        ax.set_xticks([])  # remove ticks
        ax.set_yticks([])
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 1.1)

    # third column:
    # Plot the gamma distribution for different shape and scale
    axins1 = inset_axes(axes[0, 1], bbox_to_anchor=(0.6, 0.1, 0.8, 0.8), bbox_transform=axes[0, 1].transAxes,
                        width="70%", height="70%", loc="lower left")
    axins1.plot(x, gamma(x, 2, 4.2), label='a=1, b=3', color=green, linewidth=lw+1)
    # Plot the inverse gamma distribution for different shape and scale
    axins2 = inset_axes(axes[1, 1], bbox_to_anchor=(0.6, 0.1, 0.8, 0.8), bbox_transform=axes[1, 1].transAxes,
                        width="70%", height="70%", loc="lower left")
    axins2.plot(x, inv_gamma(x, 2, 2), label='a=2, b=1', color=red, linewidth=lw+1)
    # add labels and legend

    # remove the top and right spines
    for ax in [axins1, axins2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # make the spines wider
        ax.spines['left'].set_linewidth(lw)
        ax.spines['bottom'].set_linewidth(lw)
        ax.set_xlabel('Lag time', fontsize=fsize - 3, labelpad=3)
        ax.set_ylabel('Probability', fontsize=fsize - 3, labelpad=2)
        ax.set_xticks([])  # remove ticks
        ax.set_yticks([])
        ax.set_xlim(0, 5)
        ax.set_ylim(0, 1.6)
        ax.set_title('Single-cell lag', fontsize=fsize - 2, pad=2)


def panel_B(axes):
    project_dir = os.path.dirname(os.path.dirname(root_dir))
    bins = np.linspace(0, 1000, 40)
    hist_data = np.loadtxt(os.path.join(project_dir, 'scanlag_data', 'kaplan_shx', 't0.txt'))
    hist_data2 = np.loadtxt(os.path.join(project_dir, 'scanlag_data', 'kaplan_shx', 't1346.txt'))
    x0 = hist_data[0]
    hist_data = hist_data - x0
    hist_data2 = hist_data2 - x0
    green = '#31a354'
    red = '#CD5C5C'
    axes[0,0].hist(hist_data, bins=bins, color=green, alpha=0.7)
    axes[1,0].hist(hist_data2, bins=bins, color=red, alpha=0.7)
    for ax in axes[:, 0]:
        ax.set_ylabel('Probability', fontsize=fsize-2, labelpad=0)
        ax.set_xlim(0, 1000)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    axes[0,0].set_xticks([])
    axes[1,0].set_xlabel('Lag time (min)', fontsize=fsize, labelpad=0)
def panel_E(ax):
    # scanlag plot:
    data_dir = os.path.join(os.path.dirname(os.path.dirname(root_dir)), 'scanlag_data', 'exp2')
    x_min = 0
    x_max = 4000
    n_points = 400
    n_reps = 3
    interp_data = {'Exponential': [], 'Reg-Arrest': [], 'Dis-Arrest': []}
    colors = {'Exponential': '#9ecae1', 'Reg-Arrest': '#31a354', 'Dis-Arrest': '#a50f15'}
    linestyles = {'Exponential': '-', 'Reg-Arrest': '--', 'Dis-Arrest': '-'}
    common_x = np.linspace(x_min, x_max, num=n_points)
    exp_data = pd.read_csv(os.path.join(data_dir, 'REP3EXP_t00Min_ax1.csv'),header=0)
    t0 = np.min(exp_data['X'])
    for file in os.listdir(data_dir):
        data = pd.read_csv(os.path.join(data_dir, file), header=0)
        y_interpolated = np.interp(common_x, data['X'], data['Y'])
        if 'EXP' in file:
            interp_data['Exponential'].append(y_interpolated)
        elif 'CASP' in file:
            interp_data['Reg-Arrest'].append(y_interpolated)
        elif 'SHX' in file:
            interp_data['Dis-Arrest'].append(y_interpolated)


    for key, value in interp_data.items():
        y_mean = np.mean(value, axis=0)
        y_std = np.std(value, axis=0, ddof=1)
        t_crit = t.ppf(0.84, df=n_reps - 1)  # for 68% CI: t(0.84, df = n_reps-1)
        ci = (y_std / np.sqrt(n_reps))
        existing_values = y_mean > y_mean[-1]
        ax.plot(common_x[existing_values]-t0, y_mean[existing_values], label=key, color=colors[key],
                linestyle=linestyles[key], linewidth=1)
        plt.fill_between(common_x[existing_values]-t0,
                         (y_mean - ci)[existing_values],
                         (y_mean + ci)[existing_values],
                         alpha=0.3,
                         color=colors[key])

    ax.set_xlabel('Lag time (min)', fontsize=fsize)
    ax.set_ylabel('Fraction of arrested bacteria', fontsize=fsize-1, labelpad=0)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0, 2500)
    #ax.set_xticks([10, 100, 1000])
    #ax.set_xticklabels([10, 100, 1000])
    # set tick fontsize
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax.set_ylim(0.002, 2)
    handles, labels = ax.get_legend_handles_labels()

    # use an OrderedDict to remove duplicates while preserving order
    by_label = OrderedDict(zip(labels, handles))

    # re-draw the legend with only the unique labels
    ax.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=fsize - 2)###


def panel_G(ax):
    # UMAP panel:
    # Load data
    # project directory
    project_dir = os.path.dirname(os.path.dirname(root_dir))
    data = pd.read_csv(os.path.join(project_dir, 'scanpy', 'umap_coordinates_shx_paper_barcodes_20260818_142615.csv'), index_col=0,
                       header=0)
    # scatter plot
    exp_data = data[data['batch'] == 'exp']

    reg_data = data[np.logical_or(data['batch'] == 'reg1', data['batch'] == 'reg2')]
    dis_data = data[np.logical_or(data['batch'] == 'dis1', data['batch'] == 'dis2')]
    colors = ['#4393c3', '#a6dba0', '#d6604d']
    ax.scatter(dis_data.UMAP_1, dis_data.UMAP_2, color=colors[2], alpha=.6, s=.5, label='Dis-Arrest (SHX)')
    ax.scatter(reg_data.UMAP_1, reg_data.UMAP_2, color=colors[1], alpha=.6, s=.5, label='Reg-Arrest')
    #ax.scatter(exp_data.UMAP_1, exp_data.UMAP_2, color=colors[0], alpha=.6, s=.5, label='Exponential')
    ax.legend(fontsize=fsize-2, loc='upper left', markerscale=4, frameon=False)
    ax.grid(False)
    ax.set_xlabel('UMAP1', fontsize=fsize-2)
    ax.set_ylabel('UMAP2', fontsize=fsize-2)
    #ax.set_xlim([-9,7])
    # remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def panel_H(ax):
    # UMAP panel:
    # Load data
    # project directory
    project_dir = os.path.dirname(os.path.dirname(root_dir))
    data = pd.read_csv(os.path.join(project_dir, 'scanpy', 'umap_coordinates_shx_paper_barcodes_20260818_142615.csv'), index_col=0,
                       header=0)
    # scatter plot
    colors = ['#8073ac', '#b2182b', '#d6604d', '#4393c3', '#92c5de']
    data = data[data['batch'] != 'exp']
    # color by cluster
    for i in range(max(data['cluster']) + 1):
        ax.scatter(data[data['cluster'] == i].UMAP_1, data[data['cluster'] == i].UMAP_2, color=colors[i], alpha=.6, s=.5,
                   label=f"Cluster {i + 1}")
        ax.text(data[data['cluster'] == i].UMAP_1.mean(), data[data['cluster'] == i].UMAP_2.mean(), str(i),
                fontsize=fsize, color='k', ha='center', va='center')

    ax.grid(False)
    ax.set_xlabel('UMAP1', fontsize=fsize-2)
    ax.set_ylabel('UMAP2', fontsize=fsize-2)
    #ax.set_xlim([-9,7])
    # remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


###
# Build figure 1:
fsize = 10
pf = PanelFigure(figsize=(7, 6.5), label_offset=(-0.03,0.03))
panel_pos = [
    [0.05, 0.68, 0.35, 0.28],  # A
    [0.48, 0.7, 0.22, 0.26],  # B
    [0.8, 0.74, 0.18, 0.22],  # C
    [0.05, 0.37, 0.2, 0.25],  # D
    [0.33, 0.37, 0.22, 0.25],  # E
    [0.05, 0.02, 0.47, 0.28],  # F
    [0.62, 0.37, 0.32, 0.25],  # G
    [0.62, 0.03, 0.32, 0.25],  # H
]
root_dir = os.getcwd()
# panel A:
axes_panel_A = pf.add_grid_panel(panel_pos[0], 2, 2, label="A",
                  sharex=True, sharey=True,
                  wspace=0.15, hspace=0.2)
panel_A(axes_panel_A)
# panel B:
axes_panel_B = pf.add_grid_panel(panel_pos[1], 2, 1, label="B",
                  sharex=True, sharey=True,
                  wspace=0.15, hspace=0.2)
panel_B(axes_panel_B)
# panel C:
pf.add_panel(panel_pos[2], draw_func=panel_C, label="C")
# panel D:
pf.add_panel(panel_pos[3], hide_axis=True, label="D")
# panel E:
pf.add_panel(panel_pos[4], draw_func=panel_E, label="E")
# panel F"
pf.add_panel(panel_pos[5], hide_axis=True, label="F")
# panel G:
pf.add_panel(panel_pos[6], draw_func=panel_G, label="G")
# panel H:
pf.add_panel(panel_pos[7], draw_func=panel_H, label="H")

pf.save("figure1.svg", dpi=300)
plt.show()
