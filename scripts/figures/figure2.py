import src.analysis_functions as af
import src.data_functions as df
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.datasets import make_spd_matrix
from figure_functions import PanelFigure
import numpy as np
import os
import importlib

importlib.reload(af)
importlib.reload(df)
# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------
def get_data_for_plot(path, norm=True, log=False, norm_method='sum', norm_sum=1):
    # get annotated matrix from file
    amat = df.read_from_csv(path)
    # calculate the eigenvalues and plot:
    # remove tracker genes from the matrix
    t_genes = ['16s_mature', '16s_unprocessed']
    # amat.reset_filters()
    index = [np.where(amat.var_names == val)[0][0] for val in t_genes if val in amat.var_names]
    # filter the genes
    amat.filtered_var[np.array(index).astype(int)] = False
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


def panel_A(axes):
    # Create a random NXP matrix
    N = 500  # number of rows
    P = 1000  # number of columns
    matrix = np.random.randn(N, P)
    # get correlation matrix
    corr_matrix = matrix.T @ matrix / N
    sns.heatmap(corr_matrix[:20, :20], ax=axes[0,0], cmap='bwr', center=0, cbar=True, vmin=-.5, vmax=.5)
    colorbar = axes[0,0].collections[0].colorbar
    colorbar.set_ticks([-0.5, 0, 0.5])  # Set desired ticks here
    colorbar.set_label('correlation', fontsize=fsize-2, labelpad=0)  # Set desired label here
    colorbar.ax.tick_params(labelsize=fsize-2)  # Set desired fontsize here

    # plot the MP distribution
    eigvals, _ = np.linalg.eig(corr_matrix)
    eigvals = np.real(eigvals[eigvals > 1e-6])  # keep only positive eigenvalues
    bins = np.linspace(0, 6, 40)
    axes[1,0].hist(eigvals, weights=np.ones_like(eigvals)/(P*(bins[1]-bins[0])), bins=bins, color='red', alpha=0.5, density=False, label='simulated\neigenvalues')
    x = np.linspace(0, 6, 100)
    y = np.array([af.mp_distribution(val, P / N) for val in x])
    axes[1,0].plot(x, y, color='red', linewidth=1, label='MP')
    axes[0,0].set_title('Random matrix', fontsize=fsize)
    axes[0,0].set_yticks([])
    axes[0,0].set_xticks([])
    axes[1,0].set_title('Eigenvalue density', fontsize=fsize)
    axes[1,0].set_ylabel(r'Probability density-$\rho(\lambda)$', fontsize=fsize-2,labelpad=0)
    axes[1,0].set_xlabel(r'Eigenvalue-$\lambda$', fontsize=fsize-2, labelpad=0)
    # set the font size of the ticks
    axes[1,0].tick_params(axis='both', which='major', labelsize=fsize-2)
    axes[1,0].grid(False)
    axes[1,0].legend(fontsize=fsize - 3)

def panel_B(ax):
    files = ['model_alpha2_sigma0.txt','model_alpha2_sigma07.txt', 'model_alpha2_sigma08.txt',
             'model_alpha2_sigma09.txt']
    labels = ['$\chi=0$ (MP)','$\chi=0.7$', '$\chi=0.8$', '$\chi=0.9$']
    # Get the colormap
    cmap = plt.cm.RdBu  # Choose your colormap

    # Get 4 evenly spaced colors from the colormap
    colors = [cmap(i) for i in [0.1,0.7, 0.85, 0.95]]
    for i, file in enumerate(files):
        data = np.loadtxt(os.path.join(root_dir,'model fit', file))
        ax.plot(data[:, 0], data[:, 1], label=labels[i], color=colors[i], linewidth=1.5)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize, labelpad=0)
    ax.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize, labelpad=0)
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 0.35)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    # set the font size of the ticks
    ax.tick_params(axis='both', which='major', labelsize=fsize)
    ax.set_title('Generalized MP', fontsize=fsize)
    ax.legend(fontsize=fsize)


def panel_C(axes):
    # Step 1: Generate a correlated matrix
    np.random.seed(42)
    n_samples = 15
    n_features = 10

    # Create a random positive semi-definite covariance matrix
    cov_matrix = make_spd_matrix(n_features, random_state=42)

    # Generate multivariate Gaussian data with this covariance
    data = np.random.multivariate_normal(mean=np.zeros(n_features), cov=cov_matrix, size=n_samples)
    # cluster the data
    linkage_matrix = linkage(data.T, method='average')
    ordered_cols = leaves_list(linkage_matrix)
    # reorder the data
    data = data[:, ordered_cols]
    # Step 2: Compute correlation matrix of the features (columns)
    corr_matrix = np.corrcoef(data, rowvar=False)
    # Step 3: Cluster the columns based on correlation
    linkage_matrix = linkage(corr_matrix, method='average')
    ordered_cols = leaves_list(linkage_matrix)
    reordered_corr = corr_matrix[ordered_cols][:, ordered_cols]
    background = np.full((data.shape), -np.inf)
    background[:n_features, :] = reordered_corr
    # Step 4: Plot
    sns.heatmap(data, linewidths=.5, linecolor='black', ax=axes[0, 0], cmap='Blues', center=0, cbar=False)
    sns.heatmap(background, ax=axes[0, 1], cmap='Blues', center=0, cbar=False, vmin=-1, vmax=1)
    # generate scrambled data
    scrambled_data = data.copy()
    # scramble the data
    for i in range(n_features):
        np.random.shuffle(scrambled_data[:, i])
    scrambled_corr = np.corrcoef(scrambled_data, rowvar=False)
    # Step 3: Cluster the columns based on correlation
    linkage_matrix = linkage(scrambled_corr, method='average')
    ordered_cols = leaves_list(linkage_matrix)
    scrambled_data = scrambled_data[:, ordered_cols]
    sns.heatmap(scrambled_data, linewidths=.5, linecolor='black', ax=axes[1, 0], cmap='Reds', center=0, cbar=False)
    reordered_corr = scrambled_corr[ordered_cols][:, ordered_cols]
    scrambled_background = np.full((data.shape), -np.inf)
    scrambled_background[:n_features, :] = reordered_corr
    # plot
    sns.heatmap(scrambled_background, ax=axes[1, 1], cmap='Reds', center=0, cbar=False, vmin=-1, vmax=1)
    for i in range(axes.shape[0]):
        for j in range(axes.shape[1]):
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])
            if j == 0:
                axes[i, j].set_ylabel('cells', fontsize=fsize-2)
                axes[i, j].set_xlabel('genes', fontsize=fsize-2)
                axes[i, j].xaxis.set_label_position('top')
            else:
                axes[i, j].set_xlabel('correlation matrix', fontsize=fsize-2)
                axes[i, j].xaxis.set_label_position('top')


def panel_D(ax):
    nbins = 81
    x_max = 8
    y_max = 0.3

    # plot individual distributions
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    # set style to default
    plt.style.use('default')
    plt.rcParams.update({'font.size': fsize})
    # create grid of experimental results
    # subplot 1
    # get annotated matrix from file
    file_name = ('sample_2b_filtered.csv')
    path = os.path.join(root_dir,'filtered_data', file_name)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    ax.set_title('Exponential', fontsize=fsize-2)
    plot_eigvals(ax,pcs, pcs1, N, x_max, y_max, nbins)


def panel_E(ax):
    nbins = 81
    x_max = 8
    y_max = 0.3
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    # set style to default
    plt.style.use('default')
    plt.rcParams.update({'font.size': fsize})
    # get annotated matrix from file
    file_name = ('adam_matrix_filtered.csv')
    path = os.path.join(root_dir,'filtered_data', file_name)
    ax.set_title(r'Exponential $\mathit{E. coli}$,'
                 "\n"
                 r'McNulty $\mathit{et. al.}$', fontsize=fsize-2)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, nbins)


def panel_F(ax):
    nbins = 81
    x_max = 8
    y_max = 0.3
    norm = True
    log = False
    norm_method = 'sum'
    norm_sum = 1
    # set style to default
    plt.style.use('default')
    plt.rcParams.update({'font.size': fsize})
    # get annotated matrix from file
    file_name = ('deb_KP_CDS_untreated.csv')
    path = os.path.join(root_dir, 'filtered_data', file_name)
    ax.set_title(r'Exponential $\mathit{K. pneumoniae}$,'
                 '\n'
                 r'Ma $\mathit{et. al.}$', fontsize=fsize-2)
    pcs, pcs1, N = get_data_for_plot(path, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    plot_eigvals(ax, pcs, pcs1, N, x_max, y_max, nbins)

###
# Build figure 1:
fsize = 10
plt.close("all")
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
pf = PanelFigure(figsize=(7, 5.5), label_offset=(-0.04, 0.04))
panel_pos = [
    [0.06,0.52, 0.21, 0.43],  # A
    [0.4, 0.65, 0.26, 0.3],  # B
    [0.72, 0.55, 0.25, 0.4],  # C
    [0.08, 0.1, 0.23, 0.3],  # D
    [0.4, 0.1, 0.23, 0.3],  # E
    [0.72, 0.1, 0.23, 0.3],  # F
]
# panel A:
axes_panel_A = pf.add_grid_panel(panel_pos[0], 2, 1, hspace=0.4)
panel_A(axes_panel_A)
# panel B:
pf.add_panel(panel_pos[1], draw_func=panel_B)
# panel C:
axes_panel_C = pf.add_grid_panel(panel_pos[2], 2, 2,
                  sharex=True, sharey=True,
                  wspace=0.3, hspace=0.2)
panel_C(axes_panel_C)
# panel D:
pf.add_panel(panel_pos[3], draw_func=panel_D)
# panel E:
pf.add_panel(panel_pos[4], draw_func=panel_E)
# panel F:
pf.add_panel(panel_pos[5], draw_func=panel_F)
#pf.save("figure2.svg", dpi=300)
plt.show()