import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
import seaborn as sns
#import sklearn
import numpy as np
import pandas as pd
from scipy import stats
#import statsmodels.api as sm
#from statsmodels.stats.multitest import multipletests

def scramble(m):
    # Scramble the column indices in each row of a matrix m
    #m = np.array([np.random.permutation(row) for row in m])
    # Scramble the row indices in each column of a matrix m
    m = np.array([np.random.permutation(row) for row in m.T]).T
    return m


def normalize(m, method='norm', target_sum=1):
    # Normalize the rows of a matrix m by norm or by sum
    if method == 'sum':
        m = target_sum*m / m.sum(axis=1)[:, None]
    else:
        m = target_sum*m / la.norm(m, axis=1)[:, None]
    m[np.isnan(m)] = 0
    return m


def z_transform(m):
    # Z-transform the columns of a matrix m
    m = (m - m.mean(axis=0)) / m.std(axis=0)
    m[np.isnan(m)] = 0
    return m


def log_transform(m):
    # Log-transform the elements of a matrix m
    return np.log(m + 1)


def get_pcs(m):
    # Get the principal components of a matrix m
    n = m.shape[0]
    p = m.shape[1]
    pcs = np.zeros(p)
    pcs[:min(n, p)] = la.svd(m)[1]**2/n
    return pcs


def get_eig_dist(m, norm=True, log=False, norm_method='sum', norm_sum=1):
    # get the eigenvalue distribution of the normalized matrix m
    # scramble the matrix m, and get the eigenvalue distribution of the normalized matrix
    #m = log_transform(m)  # z-transform the matrix m
    n_reps = 10
    gene_sums = (m>0).sum(axis=0)
    min_cells = 1
    m = m[:,gene_sums >= min_cells]
    min_genes = 1
    cell_sums = (m>0).sum(axis=1)
    m = m[cell_sums >= min_genes,:]
    if norm:
        m = normalize(m, method=norm_method, target_sum=norm_sum)  # normalize the rows of the matrix m
    if log:
        m = log_transform(m)  # log-transform the matrix m
    m = z_transform(m)

    pcs = get_pcs(m)  # get the principal components of the matrix m
    pcs1 = np.zeros(len(pcs))
    for _ in range(n_reps):
        m1 = m.copy()  # copy the matrix m for scrambling
        m1 = scramble(m1)  # scramble the matrix m
        #m1 = z_transform(m1) # z-transform the matrix m1
        pcs1 += get_pcs(m1)/n_reps  # get the principal components of the matrix m1
    return pcs, pcs1


def mp_distribution(x, a):
    # Marchenko-Pastur distribution with ratio a
    l_min = (1-np.sqrt(a))**2
    l_max = (1+np.sqrt(a))**2
    if l_min < x < l_max:
        f = (1/(2*np.pi*x*a))*np.sqrt((x-l_min)*(l_max-x))
    else:
        f = 0
    return f


def plot_eig_dist(pcs, pcs1, N, x_max, y_max, n_bins, ax=None, x_label=True, y_label=True):
    # plot the eigenvalue distribution of the normalized filtered matrix
    # define limits and bin number
    P = len(pcs)
    scale = 1 # scale factor for the Marchenko-Pastur distribution
    edges = np.linspace(-0.1, x_max, num=n_bins)

    # remove zeros in pcs and pcs1
    # if alpha>1 adjust the scale factor to match theoretical results
    if P/N > 1:
        scale = N/P
        pcs = pcs[pcs != 0]
        pcs1 = pcs1[pcs1 != 0]

    # create figure
    if ax is None:
        fig, ax = plt.subplots()
    # first plot
    counts, bins = np.histogram(pcs, bins=edges, density=True)
    ax.plot(bins[1:], scale*counts, color='#3182bd', linewidth=2, label='original data')
    ax.fill_between(bins[1:], scale*counts, 0, color='#9ecae1', alpha=.4)
    # second plot
    counts, bins = np.histogram(pcs1, bins=edges, density=True)
    ax.plot(bins[1:], scale*counts, color='#de2d26', linewidth=2, label='scrambled data')
    ax.fill_between(bins[1:], scale*counts, 0, color='#fc9272', alpha=.4)
    # plot analytical Marchenko-Pastur distribution
    x = np.linspace(-0.1, x_max, 100)
    y = [mp_distribution(val, P / N) for val in x]
    ax.plot(x, y, color='#756bb1', linestyle='dashed', label='MP')
    # labels and limits
    if x_label:
        ax.set_xlabel(r"$\lambda$")
    if y_label:
        ax.set_ylabel(r"$\rho(\lambda)$")
    ax.set_ylim(0, y_max)
    ax.set_xlim(0, x_max)
    # set x_ticks with difference of 2
    ax.set_xticks(np.arange(0, (x_max // 2) * 2 + 2, 2))
    # set y_ticks with difference of 0.1
    ax.set_yticks(np.arange(0, (y_max // 0.1) * 0.1 + 0.1, 0.1))
    ax.legend(facecolor='white', framealpha=1)
    # print max eigenvalue value in the plot
    #ax.text(x_max * 0.5, y_max * 0.5, r"$\lambda_{max}$: "+f" {round(max(pcs), 2)}")
    # change font size of labels and axes
    return ax

def calculate_entropy(DE_timepoints, reg_param):
    model = sklearn.covariance.GraphicalLasso(alpha=reg_param,max_iter=1000, tol=0.001)
    model.fit(DE_timepoints)
    sign, entropy = np.linalg.slogdet(model.covariance_)
    # Plot the clustered heatmap
    # We use a diverging colormap ('vlag' or 'coolwarm') centered at 0 for covariance
    g = sns.clustermap(
        model.covariance_,
        cmap="viridis",
        annot=False,  # Set to False if your matrix is very large
        figsize=(6, 6),
        xticklabels=False,
        yticklabels=False,
        linewidths=0,
        cbar_kws={'label': 'Covariance'}
    )
    # Hide the hierarchical clustering trees (dendrograms)
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    g.cax.set_position([1.02, 0.3, 0.03, 0.4])
    # Adjust layout and show the plot
    g.fig.suptitle(fr'GLASSO regularized Covariance Matrix $\rho={reg_param}$', y=0.85)
    plt.show()
    return entropy


def run_de_analysis(data, metadata, condition_col='condition', test_time=1):
    """
    Performs DE analysis on GMM rate data.

    Parameters:
    data: pd.DataFrame (rows=reactions/genes, cols=samples)
    metadata: pd.DataFrame (rows=samples, must contain condition_col)
    """
    results = []
    genes = data.index

    # 1. Log2 Transformation (Handling zeros if necessary)
    # Since rates are continuous and > 0, log-transforming stabilizes variance
    log_data = np.log2(data + 1e-9)

    # Identify groups
    groups = metadata[condition_col].unique()
    cond1 = metadata[metadata[condition_col] == groups[0]].index
    cond2 = metadata[metadata[condition_col] == groups[test_time]].index

    for gene in genes:
        sample1 = log_data.loc[gene, cond1]
        sample2 = log_data.loc[gene, cond2]

        # 2. Calculate Log2 Fold Change
        l2fc = sample2.mean() - sample1.mean()

        # 3. Statistical Test
        # We use Welch's t-test (doesn't assume equal variance)
        # because GMM parameters can cause wild variance differences
        t_stat, p_val = stats.ttest_ind(sample2, sample1, equal_var=False)

        results.append({
            'gene': gene,
            'log2FoldChange': l2fc,
            'pvalue': p_val,
            'stat': t_stat
        })

    res_df = pd.DataFrame(results).set_index('gene')

    # 4. Multiple Testing Correction (Benjamini-Hochberg)
    # Crucial for associating reliable significance to the L2FC
    res_df['padj'] = multipletests(res_df['pvalue'], method='fdr_bh')[1]

    return res_df.sort_values('padj')



def plot_volcano(res_df, lfc_thresh=1.0, p_thresh=0.05):
    """
    Plots a volcano plot from the DE analysis results.
    """
    # 1. Prepare data
    plot_df = res_df.copy()
    plot_df['-log10p'] = -np.log10(plot_df['padj'])

    # 2. Define Significance Categories
    plot_df['group'] = 'Not Significant'
    plot_df.loc[(plot_df['log2FoldChange'] > lfc_thresh) & (plot_df['padj'] < p_thresh), 'group'] = 'Up-regulated'
    plot_df.loc[(plot_df['log2FoldChange'] < -lfc_thresh) & (plot_df['padj'] < p_thresh), 'group'] = 'Down-regulated'

    # 3. Create Plot
    plt.figure(figsize=(6, 6))
    sns.set_style("whitegrid")

    colors = {'Not Significant': 'grey', 'Up-regulated': '#d62728', 'Down-regulated': '#1f77b4'}

    sns.scatterplot(
        data=plot_df, x='log2FoldChange', y='-log10p',
        hue='group', palette=colors, alpha=0.7, edgecolor=None
    )

    # 4. Add Threshold Lines
    plt.axhline(-np.log10(p_thresh), color='black', linestyle='--', alpha=0.5)
    plt.axvline(lfc_thresh, color='black', linestyle='--', alpha=0.5)
    plt.axvline(-lfc_thresh, color='black', linestyle='--', alpha=0.5)

    # 5. Aesthetics
    plt.title('Differential Expression Volcano Plot', fontsize=15)
    plt.xlabel(r'$\log_2$ Fold Change', fontsize=12)
    plt.ylabel(r'$-\log_{10}$ Adjusted P-value', fontsize=12)
    plt.legend(title='Status', loc='upper right')
    plt.xlim([-8,4])
    # Annotate top 5 most significant hits
    plt.show()
