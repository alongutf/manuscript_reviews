import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from figure_functions import PanelFigure
import numpy as np
import pandas as pd
import os
from scipy.stats import t
import src.analysis_functions as af
import src.data_functions as df

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
    ax.plot(bins[1:], scale * counts, color='#3182bd', linewidth=.75, label='original data')
    ax.fill_between(bins[1:], scale * counts, 0, color='#9ecae1', alpha=.4)
    # second plot
    counts, bins = np.histogram(pcs1, bins=edges, density=True)
    ax.plot(bins[1:], scale * counts, color='#de2d26', linewidth=.5, label='scrambled data')
    ax.fill_between(bins[1:], scale * counts, 0, color='#fc9272', alpha=.4)
    # plot analytical Marchenko-Pastur distribution
    x = np.linspace(-0.1, x_max, 100)
    y = [af.mp_distribution(val, P / N) for val in x]
    ax.plot(x, y, color='#756bb1', linestyle='dashed', label='MP')
    # labels and limits
    if x_label:
        ax.set_xlabel("$\lambda$", fontsize=fsize, labelpad=0)
    if y_label:
        ax.set_ylabel(r"$\rho(\lambda)$", fontsize=fsize, labelpad=0)
    ax.set_ylim(0, y_max)
    ax.set_xlim(0, x_max)
    # set x_ticks with difference of 2
    ax.set_xticks(np.arange(0, (x_max // 2) * 2 + 2, 2))
    # set y_ticks with difference of 0.1
    ax.set_yticks([0.1,0.2])
    ax.legend(facecolor='white', framealpha=1, fontsize=fsize-2, loc='upper right')
    # set the font size of the ticks
    ax.tick_params(axis='both', which='major', labelsize=fsize)


def panel_B(ax):
    data = pd.read_csv(os.path.join(root_dir, 'scanpy', 'umap_coordinates_vapc.csv'), index_col=0,
                       header=0)

    # scatter plot
    exp_data = data[data['batch'] == 'exp']
    t2_data = data[data['batch'] == 'T2']
    t5a_data = data[np.logical_or(data['batch'] == 'T5A', data['batch'] == 'T5B')]
    dis_data = data[data['batch'] == 'TON']
    colors = ['#4393c3', '#f4a582', '#d6604d', '#b2182b']
    ax.scatter(exp_data.UMAP_1, exp_data.UMAP_2, color=colors[0], alpha=.8, s=.5, label='Exponential')
    ax.scatter(t2_data.UMAP_1, t2_data.UMAP_2, color=colors[1], alpha=.8, s=.5, label='VapC: 2h')
    ax.scatter(t5a_data.UMAP_1, t5a_data.UMAP_2, color=colors[2], alpha=.8, s=.5, label='VapC: 5h')
    ax.scatter(dis_data.UMAP_1, dis_data.UMAP_2, color=colors[3], alpha=.8, s=.5, label='VapC: 24h')
    ax.legend(fontsize=fsize-2, loc='lower right', bbox_to_anchor=(1.2,-0.2), markerscale=4, frameon=False)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

def panel_C(ax):
    # UMAP panel:
    # Load data
    # project directory
    data = pd.read_csv(os.path.join(root_dir, 'scanpy', 'umap_coordinates_vapc.csv'), index_col=0,
                       header=0)
    # scatter plot
    colors = ['#8073ac', '#b2182b', '#4393c3', '#d6604d', '#92c5de']
    # color by cluster
    for i in range(max(data['cluster']) + 1):
        ax.scatter(data[data['cluster'] == i].UMAP_1, data[data['cluster'] == i].UMAP_2, color=colors[i], alpha=.6, s=.5,
                   label=f"Cluster {i + 1}")
        ax.text(data[data['cluster'] == i].UMAP_1.mean(), data[data['cluster'] == i].UMAP_2.mean(), str(i),
                fontsize=fsize, color='k', ha='center', va='center')
    #remove axis spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

def panel_D(ax):
    file_name = 'VapC_biorep_t2A_filtered.csv'
    nbins = 81
    x_max = 8
    y_max = 0.27
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    path = os.path.join(root_dir, 'data_for_paper', file_name)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    ax.set_title('Early VapC (2h)', fontsize=fsize)
    plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, nbins)
    ax.set_ylabel(r"Probability Density - $\rho(\lambda)$", fontsize=fsize, labelpad=0)
    ax.set_xlabel(r"Eigenvalue - $\lambda$", fontsize=fsize, labelpad=0)

def panel_E(ax):
    file_name = 'VapC_biorep_tONA_filtered.csv'
    nbins = 81
    x_max = 8
    y_max = 0.25
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    path = os.path.join(root_dir, 'data_for_paper', file_name)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    ax.set_title('Late VapC (24h)', fontsize=fsize)
    plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, nbins)

def panel_F(ax):
    # Create bar plot graph from model fit:
    # load data
    dataset_values = os.path.join(root_dir,'scripts','figures','figure5','dataset_summary_no_plasmid_genes.csv')
    scrambled_values = os.path.join(root_dir, 'scripts', 'figures', 'figure3', 'scrambled_GMPcor.csv')
    # read data
    data = pd.read_csv(dataset_values)
    scrambled = pd.read_csv(scrambled_values)
    samples = {'EXP_biorep_t0A','VapC_biorep_t2A','VapC_biorep_t5A', 'VapC_biorep_tONA'}
    labels = ['Exponential','VapC\n2h','VapC\n5h','VapC\n24h']
    # get data for the samples
    data = data[data['sample'].isin(samples)]
    data['time'] = [20, 2, 5, 0]
    # sort data
    colors = ['#4393c3', '#f4a582', '#d6604d', '#b2182b', "#bdbdbd"]
    data = data.sort_values('time', ascending=True)
    bar_width = 0.25
    gap_between_bars = 0.4
    # add the scrambled values
    labels.append('Scrambled')
    means = data['sigma fit'].to_list()
    errors = data['standard error'].to_list()
    means.append(scrambled['GMP-Cor'].mean())
    errors.append(scrambled['GMP-Cor'].std() / np.sqrt(len(scrambled)))  # standard error
    positions = [i * (bar_width + gap_between_bars) for i in range(len(means))]
    # create bar plot of sigma fit values
    positions[0] -= gap_between_bars / 2  # adjust first bar position
    positions[-1] += gap_between_bars / 2 # adjust last bar position
    ax.bar(positions, means, yerr=errors, capsize=2.5, color=colors, edgecolor='black',
           alpha=0.7, width=bar_width)

    # set positions of bars
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=0, fontsize=fsize-2, ha='center')
    ax.set_ylabel('GMP-Cor', fontsize=fsize, labelpad=0)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8])
    ax.tick_params(axis='both', which='major', labelsize=fsize-2)
    ax.set_ylim([0, 1])
    ax.set_xlim([positions[0] - gap_between_bars/2, positions[-1] + gap_between_bars/2])


def panel_G(ax):
    path = os.path.join(root_dir,'scripts','figures','figure5','normalizedOD_at_20h.csv')
    data = pd.read_csv(path, header=0,
                       index_col=0)

    data_dict = {}  # will hold (group_mean, group_SE) for plotting
    bio_means = {}  # stores list of biological-replicate means per group (for tests)

    for col in data.columns:
        X = data[col].dropna().to_numpy()
        means = []
        errors = []
        for i in range(int(len(X) / 2)):
            biorep = X[2 * i:2 * i + 2]
            means.append(np.mean(biorep))
            errors.append(np.std(biorep, ddof=1) / np.sqrt(2))
        bio_means[col] = means[:]  # keep biological-replicate means

        # your original rule to choose SE (biological vs technical)
        if np.std(means, ddof=1) < np.mean(errors):
            data_dict[col] = (np.mean(means), np.mean(errors))
        else:
            data_dict[col] = (np.mean(means), np.std(means, ddof=1) / np.sqrt(len(means)))

    # --- plotting ---
    labels = list(data_dict.keys())
    plot_colors = ['#2166ac', '#92c5de', '#f4a582', '#d6604d', '#b2182b']

    bar_means = [data_dict[col][0] for col in labels]
    bar_errs = [data_dict[col][1] for col in labels]
    bar_labels = ['Reg\n2h', 'Reg\n24h', 'VapC\n2h', 'VapC\n5h', 'VapC\n24h']
    bars = ax.bar(bar_labels, bar_means, yerr=bar_errs, capsize=3,
                  color=plot_colors[:len(labels)], alpha=0.8, edgecolor='black', width=0.4)

    ax.set_ylabel('Normalized OD', fontsize=fsize, labelpad=0)

    ax.set_yscale('log')
    ax.set_yticks([0.1,1])
    ax.set_yticklabels(['0.1','1'])
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax.set_title('SDS added to culture', fontsize=fsize-2, pad=0)
    ax.margins(y=0.2)  # add headroom for significance caps

    # --- significance helpers ---
    def get_pval(means, errors, n1, n2):
        tstat = (means[0] - means[1]) / np.sqrt(errors[0] ** 2 + errors[1] ** 2)
        df = n1 + n2 - 2
        return t.sf(tstat, df)

    def p_to_stars(p):
        if np.isnan(p):   return 'n.s.'
        if p < 1e-4:      return '****'
        if p < 1e-3:      return '***'
        if p < 0.01:      return '**'
        if p < 0.05:      return '*'
        return 'n.s.'

    def add_sig_between(ax, bar1, bar2, text, y_mult=1.2, cap_mult=1.1, text_mult=1.01):
        x1 = bar1.get_x() + bar1.get_width() * 0.6
        x2 = bar2.get_x() + bar2.get_width() * 0.4
        y_base = max(bar1.get_height(), bar2.get_height())
        y = y_base * y_mult
        y_cap = y * cap_mult
        ax.plot([x1, x1, x2, x2], [y, y_cap, y_cap, y], lw=1, c='black')
        ax.text((x1 + x2) / 2, y_cap * text_mult, text, ha='center', va='bottom', fontsize=fsize-2)

    # --- compute Welch t-tests between neighbor bars and annotate ---
    for i in range(len(labels) - 1):
        g1, g2 = labels[i], labels[i + 1]
        means = [data_dict[g1][0], data_dict[g2][0]]
        errors = [data_dict[g1][1], data_dict[g2][1]]

        # Welch t-test on biological-replicate means (needs n>=2 ideally)
        p = get_pval(means, errors, len(bio_means[g1]), len(bio_means[g2]))
        print(p)
        add_sig_between(ax, bars[i], bars[i + 1], p_to_stars(p), y_mult=1.3, cap_mult=1.2)



def panel_H(ax):
    # vapc lag time distribution
    # Load the data
    conditions = ['CTRLt0', 'VAPCt240', 'VAPCt1400']
    labels = ['Reg-Arrest', 'Early VapC', 'Late VapC']
    data = {}
    colors = ['#9ecae1', '#fb6a4a', '#a50f15']
    plt.style.use('default')

    for condition in conditions:
        path = os.path.join(root_dir, 'scripts', 'figures', 'figure5', condition + '.csv')
        data[condition] = pd.read_csv(path, index_col=False, header=None).to_numpy()
    # plot histograms of lag time
    edges = np.linspace(300, 1100, 51)
    for i, condition in enumerate(conditions):
        x = data[condition].flatten() + 400
        ax.hist(x, bins=edges, color=colors[i], histtype='stepfilled', edgecolor='k', alpha=0.5, label=labels[i],
                 density=True)
    ax.set_xlabel('Lag time (min)', fontsize=fsize)
    ax.set_ylabel(r'Frequency', labelpad=0, fontsize=fsize)
    ax.text(310, 0.0175, r'$\times{10}^{-2}$', fontsize=fsize - 3)
    # set axis label size
    ax.legend(fontsize=fsize-2)
    ax.set_xlim([300, 1100])
    ax.set_xticks([300, 700, 1100])
    ax.set_yticks([0, 0.01, 0.02])
    ax.set_yticklabels([0, 1, 2])
    ax.tick_params(axis='both', which='major', labelsize=fsize)

###
# Build figure 5:
fsize = 10
plt.close("all")
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
pf = PanelFigure(figsize=(7, 6), label_offset=(0, 0.03))
panel_pos = [
    [0.075, 0.6, 0.5, 0.34],  # A
    [0.075, 0.45, 0.24, 0.2],  # B
    [0.4, 0.45, 0.24, 0.2],  # C
    [0.075, 0.08, 0.24, 0.28],  # D
    [0.4, 0.08, 0.24, 0.28],  # E
    [0.7, 0.72, 0.275, 0.22],  # F
    [0.7, 0.41, 0.275, 0.22],  # G
    [0.7, 0.08, 0.275, 0.22],  # H
]
# panel A:
pf.add_panel(panel_pos[0], hide_axis=True)
# panel B:
pf.add_panel(panel_pos[1], draw_func=panel_B)
# panel C:
pf.add_panel(panel_pos[2], draw_func=panel_C)
# panel C:
pf.add_panel(panel_pos[3], draw_func=panel_D)
# panel D:
pf.add_panel(panel_pos[4], draw_func=panel_E)
# panel E:
pf.add_panel(panel_pos[5], draw_func=panel_F)
# panel F:
pf.add_panel(panel_pos[6], draw_func=panel_G)
# panel G:
pf.add_panel(panel_pos[7], draw_func=panel_H)
#pf.save("figure5.svg", dpi=300)
plt.show()