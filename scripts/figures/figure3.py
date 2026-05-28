import src.analysis_functions as af
import src.data_functions as df
import matplotlib.pyplot as plt
from figure_functions import PanelFigure
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ttest_ind
from scipy.stats import t
import os
import importlib
from collections import OrderedDict

importlib.reload(af)
importlib.reload(df)
# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------
def get_data_for_plot(path, norm=True, log=False, norm_method='sum', norm_sum=1):
    # get annotated matrix from file
    amat = df.read_from_csv(path)
    # calculate the eigenvalues and plot:
    m = amat.get_filtered_matrix().m
    pcs, pcs1 = af.get_eig_dist(m, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    return pcs, pcs1, m.shape[0]


def plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, n_bins, x_label=True, y_label=True):
    # plot the eigenvalue distribution of the normalized filtered matrix
    # define limits and bin number
    P = len(pcs)
    scale = 1  # scale factor for the Marchenko-Pastur distribution
    edges = np.linspace(-0.1, x_max, num=n_bins)

    # remove zeros in pcs and pcs1
    # if alpha>1 adjust the scale factor to match theoretical results
    if P / N > 1:
        scale = N / P
        pcs = pcs[pcs != 0]
        pcs1 = pcs1[pcs1 != 0]

    # first plot
    counts, bins = np.histogram(pcs, bins=edges, density=True)
    ax.plot(bins[1:], scale * counts, color='#3182bd', linewidth=1, label='original data')
    ax.fill_between(bins[1:], scale * counts, 0, color='#9ecae1', alpha=.4)
    # second plot
    counts, bins = np.histogram(pcs1, bins=edges, density=True)
    ax.plot(bins[1:], scale * counts, color='#de2d26', linewidth=1, label='scrambled data')
    ax.fill_between(bins[1:], scale * counts, 0, color='#fc9272', alpha=.4)
    # plot analytical Marchenko-Pastur distribution
    x = np.linspace(-0.1, x_max, 100)
    y = [af.mp_distribution(val, P / N) for val in x]
    ax.plot(x, y, color='#756bb1', linestyle='dashed', label='MP')
    # labels and limits
    if x_label:
        ax.set_xlabel("$\lambda$", fontsize=fsize)
    if y_label:
        ax.set_ylabel(r"$\rho(\lambda)$", fontsize=fsize)
    ax.set_ylim(0, y_max)
    ax.set_xlim(0, x_max)
    # set x_ticks with difference of 2
    ax.set_xticks(np.arange(0, (x_max // 2) * 2 + 2, 2))
    # set y_ticks with difference of 0.1
    ax.set_yticks(np.arange(0, (y_max // 0.1) * 0.1 + 0.1, 0.1))
    ax.legend(facecolor='white', framealpha=1, fontsize=fsize-2, loc='upper right')
    # set the font size of the ticks
    ax.tick_params(axis='both', which='major', labelsize=fsize)


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
    nbins = 81
    x_max = 8
    y_max = 0.3
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    plt.rcParams.update({'font.size': fsize})
    file_name = ('sample_2b_filtered.csv')
    path = os.path.join(root_dir, 'data_for_paper', file_name)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    ax.set_title('Exponential', fontsize=fsize)
    plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, nbins)
    ax.set_xlabel(r"Eigenvalue - $\lambda$", fontsize=fsize, labelpad=0)
    ax.set_ylabel(r"Probability Density - $\rho(\lambda)$", fontsize=fsize, labelpad=0)


def panel_B(ax):
    nbins = 81
    x_max = 8
    y_max = 0.3
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    plt.rcParams.update({'font.size': fsize})
    file_name = ('sample_15b_filtered.csv')
    path = os.path.join(root_dir, 'data_for_paper', file_name)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    ax.set_title('Reg-Arrest', fontsize=fsize)
    plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, nbins)


def panel_C(ax):
    nbins = 81
    x_max = 8
    y_max = 0.3
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    plt.rcParams.update({'font.size': fsize})
    file_name = ('sample_15a_filtered.csv')
    path = os.path.join(root_dir, 'data_for_paper', file_name)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    ax.set_title('Dis-Arrest', fontsize=fsize)
    plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, nbins)


def panel_D(ax):
    # set theme
    plt.style.use('default')
    # plot data with model fit:
    # load data
    fit_path = os.path.join(root_dir, 'model fit', 'fit_sample_2b.txt')
    data_path = os.path.join(root_dir, 'model fit', 'pc_data', 'sample_2b_filtered.csv')
    # get fit data
    fit_data = pd.read_csv(fit_path, sep='\t', header=None)
    # get data
    data = pd.read_csv(data_path)
    # remove zeros
    data = data[data['Y'] > 0]
    fit_data = fit_data[fit_data[0] > 0]
    # plot data and fit
    ax.scatter(data['X'], data['Y'], s=5, color='#de2d26', label='exponential data')
    ax.plot(fit_data[1], fit_data[0], color='#3182bd', label='GMP fit', linewidth=1.5)
    ax.set_xlabel('$\lambda$', fontsize=fsize)
    ax.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 0.35)
    # axis ticks
    ax.set_yticks(np.arange(0, 0.4, 0.1))
    # set tick fontsize
    ax.tick_params(axis='both', which='major', labelsize=fsize)
    ax.legend(fontsize=fsize-2, loc='upper right')


def panel_E(ax):
    # Create bar plot graph from model fit:
    # load data
    dataset_values = os.path.join(root_dir, 'scripts', 'figures', 'figure3', 'dataset_summary_no_plasmid_genes.csv')
    scrambled_values = os.path.join(root_dir, 'scripts', 'figures', 'figure3', 'scrambled_GMPcor.csv')
    # read data
    data = pd.read_csv(dataset_values)
    scrambled = pd.read_csv(scrambled_values)
    samples = {'sample_2b', 'sample_15b', 'sample_15a'}
    labels = {'sample_2b': 'Exponential','sample_13b': 'Reg-Arrest', 'sample_15b': 'Reg-Arrest','sample_13a': 'Dis-Arrest', 'sample_15a': 'Dis-Arrest'}
    colors = ["#9ecae1", "#9ecae1", "#a50f15","#bdbdbd"]
    # get data for the samples
    data = data[data['sample'].isin(samples)]
    # sort data
    data = data.sort_values('sigma fit', ascending=False)
    sigma_dict = {'Exponential': [], 'Reg-Arrest': [], 'Dis-Arrest': []}
    error_dict = {'Exponential': [], 'Reg-Arrest': [], 'Dis-Arrest': []}
    # get the sigma fit values
    for sample in samples:
        sigma_dict[labels[sample]].append(data[data['sample'] == sample]['sigma fit'].values[0])
        error_dict[labels[sample]].append(data[data['sample'] == sample]['standard error'].values[0])
    # create bar plot from dictionary
    labels = list(sigma_dict.keys())
    means = [np.mean(values) for values in sigma_dict.values()]
    errors = [np.mean(values) for values in error_dict.values()]  # standard error
    # add the scrambled values
    labels.append('Scrambled')
    means.append(scrambled['GMP-Cor'].mean())
    errors.append(scrambled['GMP-Cor'].std() / np.sqrt(len(scrambled)))  # standard error
    ax.bar(labels, means, yerr=errors, capsize=5,
           color=colors, edgecolor='black', alpha=0.7, width=0.3)

    ax.set_ylim([0, 1])
    ax.set_ylabel('GMP-Cor', fontsize=fsize)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=25)
    ax.tick_params(axis='both', which='major', labelsize=fsize)


def panel_F(ax):
    path = os.path.join(root_dir, 'scripts', 'figures', 'figure3', 'dataset_summary_no_plasmid_genes.csv')
    data = pd.read_csv(path, header=0)
    ranking_param = 'sigma fit'
    # Sort values and compute rank
    data['Rank'] = data[ranking_param].rank(method='min').astype(int)
    # Split data into two groups based on labels
    group1 = data[data['class'] == 'reg'][ranking_param]
    group0 = data[data['class'] == 'dis'][ranking_param]
    c = "k"
    box = ax.boxplot([group1, group0], meanline=True, showmeans=True, patch_artist=True,
                     boxprops=dict(facecolor="None", color=c), whiskerprops=dict(color=c), capprops=dict(color=c),
                     flierprops=dict(markeredgecolor=c, markersize=2), medianprops=dict(color=c))
    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        for item in box[element]:
            item.set_linewidth(1)
    for mean_line in box['means']:
        mean_line.set_linewidth(1)
        mean_line.set_color(c)
        mean_line.set_linestyle('solid')

    ax.set_xticklabels(['Regulated', 'Dis-Arrest'], fontsize=fsize, rotation=0, ha='center')
    ax.set_ylabel('Enrichment score', fontsize=fsize, labelpad=0)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
    ax.set_ylim([0.2, 1])
    ax.set_yticklabels([0.2, 0.4, 0.6, 0.8, 1], fontsize=fsize - 2)
    # add u-test results
    u_stat, u_p = stats.mannwhitneyu(group1, group0)

    # Define the level of significance
    asterisks = format_p(u_p)
    # Add significance annotation
    x1, x2 = 1, 2  # x-coordinates of the box plots
    y, h = max(max(group1), max(group0))+0.02, 0.02  # y-position and height of the annotation
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.75, color='black')
    ax.text((x1 + x2) * 0.5, y + h, asterisks, ha='center', va='bottom', color='black', fontsize=fsize - 2)
    ax.set_ylabel(r'GMP-Cor', fontsize=fsize)


# Build figure 1:
fsize = 10
plt.close("all")
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
pf = PanelFigure(figsize=(7, 5), label_offset=(0, 0.04))
panel_pos = [
    [0.1, 0.6, 0.23, 0.35],  # A
    [0.42, 0.6, 0.23, 0.35],  # B
    [0.74, 0.6, 0.23, 0.35],  # C
    [0.1, 0.12, 0.23, 0.35],  # D
    [0.42, 0.12, 0.23, 0.35],  # E
    [0.74, 0.12, 0.23, 0.35]  # F
]
# panel A:
pf.add_panel(panel_pos[0], draw_func=panel_A)
# panel B:
pf.add_panel(panel_pos[1], draw_func=panel_B)
# panel C:
pf.add_panel(panel_pos[2], draw_func=panel_C)
# panel D:
pf.add_panel(panel_pos[3], draw_func=panel_D)
# panel E:
pf.add_panel(panel_pos[4], draw_func=panel_E)
# panel F:
pf.add_panel(panel_pos[5], draw_func=panel_F)
#pf.save("figure3.svg", dpi=300, transparent=True)
plt.show()