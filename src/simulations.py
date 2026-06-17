import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm, nbinom, kstest, pearsonr
from scipy import stats as scipy_stats
from sklearn.covariance import GraphicalLasso
from statsmodels.stats.multitest import multipletests


# ── Covariance matrix generation ────────────────────────────────────────────

def generate_gram_hub_matrix(n, alpha, shape, hub_probability, seed=42):
    """
    Build an n×n correlation matrix with cluster + hub factor structure.

    Parameters
    ----------
    n               : number of variables (genes)
    alpha           : shared-variance fraction; controls off-diagonal strength
    shape           : Pareto shape for cluster-size distribution (heavier tail = larger clusters)
    hub_probability : probability that a cluster hub connects to the global hub factor
    seed            : RNG seed for reproducibility
    """
    np.random.seed(seed)
    # Loading matrix: identity part gives each gene unique (noise) variance
    A = np.eye(n) * np.sqrt(1 - alpha)

    remaining_vars = n
    current_idx = 0
    hubs = []

    # Cluster factors — sizes drawn from a Pareto distribution
    while remaining_vars > 0:
        size = int(np.random.pareto(shape) + 1)
        size = max(1, min(size, remaining_vars))
        idx_set = np.arange(current_idx, current_idx + size)

        cluster_col = np.zeros((n, 1))
        strength = np.random.uniform(0.1 * alpha, alpha)
        cluster_col[idx_set] = np.sqrt(strength)

        A = np.hstack([A, cluster_col])
        hubs.append(np.random.choice(idx_set))

        current_idx += size
        remaining_vars -= size

    # Global hub factor connecting a random subset of cluster hubs
    global_hub_col = np.zeros((n, 1))
    for h_idx in hubs:
        if np.random.rand() < hub_probability:
            weight = np.random.uniform(-alpha, alpha)
            global_hub_col[h_idx] = weight
    A = np.hstack([A, global_hub_col])

    C = A @ A.T
    d = np.sqrt(np.diag(C))
    R = C / np.outer(d, d)
    return R


# ── scRNA-seq count simulation ───────────────────────────────────────────────

def simulate_scRNA_data(n_cells=1000, n_genes=2000, sigma=None, rho=0.9, dropout_rate=2, inv_gamma_shape=1.5, inv_gamma_scale=0.01, seed=None):
    """
    Simulate single-cell RNA-seq count data with realistic noise.

    Parameters
    ----------
    n_cells      : number of cells
    n_genes      : number of genes
    sigma        : pre-computed correlation matrix; generated with rho if None
    rho          : shared-variance fraction passed to generate_gram_hub_matrix
    dropout_rate : controls dropout probability: P(dropout) = exp(-dropout_rate * count)
    seed         : RNG seed (int or None)

    Returns
    -------
    true_counts     : counts before dropout (n_cells × n_genes)
    observed_counts : counts after dropout  (n_cells × n_genes)
    """
    if seed is not None:
        np.random.seed(seed)

    if sigma is None:
        sigma = generate_gram_hub_matrix(n_genes, rho, 1.5, 0.2, seed=31)

    # Latent correlated normal → uniform via CDF → negative-binomial counts
    mu_vec = np.zeros(n_genes)
    latent_data = np.random.multivariate_normal(mu_vec, sigma, size=n_cells)
    uniform_data = norm.cdf(latent_data)

    counts = np.zeros_like(uniform_data)
    rng = np.random.default_rng(seed)
    for i in range(n_genes):
        # Gene mean drawn from inverse-Gamma (heavy-tailed expression distribution)
        gene_mu = 1.0 / rng.gamma(inv_gamma_shape, 1.0 / inv_gamma_scale, 1)
        gene_r = 0.5  # NB dispersion typical of scRNA-seq
        p_param = gene_r / (gene_r + gene_mu)
        counts[:, i] = nbinom.ppf(uniform_data[:, i], gene_r, p_param)

    # Cell-level amplitude variability (library-size differences)
    cell_amp = np.random.lognormal(0.5, 1, n_cells)
    cell_amp = np.maximum(cell_amp, 0)
    counts = np.round(cell_amp[:, np.newaxis] * counts)

    # Expression-dependent dropout
    dropout_prob = np.exp(-dropout_rate * counts)
    mask = np.random.rand(*counts.shape) > dropout_prob
    observed_counts = counts * mask

    return counts, observed_counts


# ── Distribution validation ──────────────────────────────────────────────────

def plot_ccdf(data, label=None, ax=None):
    """Log-log complementary CDF plot for a 1-D array."""
    sorted_data = np.sort(data)
    n = len(sorted_data)
    ccdf = 1 - np.arange(1, n + 1) / n + (1 / n)

    if ax is None:
        ax = plt.gca()
    ax.loglog(sorted_data, ccdf, marker='.', linestyle='none', label=label)
    ax.set_xlabel('Value (log)')
    ax.set_ylabel('P(X > x) (log)')
    ax.set_title('Complementary CDF (Log-Log Scale)')
    ax.grid(True, which='both', ls='-', alpha=0.5)
    return ax


def validate_cell_expression_distribution(counts, ax=None):
    """
    Compare per-cell total expression CCDF against data.

    Plots the empirical CCDF of gene totals and overlays a reference power-law
    slope, then returns the KS statistic against an exponential null (a rough
    check that the tail is heavy).
    """
    cell_totals = counts.sum(axis=1)
    cell_totals = cell_totals[cell_totals > 0]
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    ax = plot_ccdf(cell_totals, label='cell totals', ax=ax)
    ax.set_title('Cell expression CCDF')
    ax.set_xlim(10,2000)
    # KS test: heavy tail means it deviates significantly from exponential
    ks_stat, ks_p = kstest(cell_totals, 'expon',
                           args=(cell_totals.min(), cell_totals.mean()))
    return {'ks_stat': ks_stat, 'ks_p': ks_p, 'n_cells': len(cell_totals)}


def validate_gene_expression_distribution(counts, n_sample=200, ax=None):
    """
    Compare per-gene total expression CCDF against an inverse-Gamma expectation.

    Plots the empirical CCDF of gene totals and overlays a reference power-law
    slope, then returns the KS statistic against an exponential null (a rough
    check that the tail is heavy).
    """
    gene_totals = counts.sum(axis=0)
    gene_totals = gene_totals[gene_totals > 0]
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    ax = plot_ccdf(gene_totals, label='gene totals', ax=ax)
    ax.set_title('Gene expression CCDF')

    # KS test: heavy tail means it deviates significantly from exponential
    ks_stat, ks_p = kstest(gene_totals, 'expon',
                           args=(gene_totals.min(), gene_totals.mean()))
    return {'ks_stat': ks_stat, 'ks_p': ks_p, 'n_genes': len(gene_totals)}


def validate_dropout_rate(true_counts, observed_counts):
    """
    Compare observed zero fraction against the expected Poisson-dropout model.

    Returns a dict with:
      - observed_zero_fraction
      - expected_zero_fraction  (mean exp(-lambda) across all entries)
      - correlation between count magnitude and expression (should be positive)
    """
    obs_zero_frac = np.mean(observed_counts == 0)
    # Expected zeros under the dropout model: E[exp(-rate * count)]
    exp_zero_frac = np.mean(np.exp(-observed_counts))

    # Per-gene: mean expression vs fraction of non-zero cells
    gene_means = true_counts.mean(axis=0)
    detection_rate = np.mean(observed_counts > 0, axis=0)
    corr, _ = pearsonr(gene_means, detection_rate)

    return {
        'observed_zero_fraction': obs_zero_frac,
        'expected_zero_fraction': exp_zero_frac,
        'mean_detection_rate': detection_rate.mean(),
        'expression_detection_correlation': corr,
    }


def validate_count_distribution(counts, n_genes_sample=10, seed=0, ax=None):
    """
    Spot-check that per-gene count distributions look negative-binomial.

    Plots histograms of `n_genes_sample` randomly chosen genes and reports
    the fraction with variance > mean (overdispersion, expected for NB counts).
    """
    rng = np.random.default_rng(seed)
    gene_means = counts.mean(axis=0)
    gene_vars = counts.var(axis=0)
    overdispersed_frac = np.mean(gene_vars > gene_means)

    sample_idx = rng.choice(counts.shape[1], size=n_genes_sample, replace=False)
    if ax is None:
        fig, axes = plt.subplots(2, n_genes_sample // 2, figsize=(n_genes_sample * 1.5, 4))
        axes = axes.flatten()
    else:
        axes = [ax] * n_genes_sample

    for i, gi in enumerate(sample_idx):
        axes[i].hist(counts[:, gi], bins=20, color='steelblue', edgecolor='k', alpha=0.7)
        axes[i].set_title(f'gene {gi}', fontsize=8)
        axes[i].set_xlabel('counts', fontsize=7)
        axes[i].set_ylabel('freq', fontsize=7)

    plt.suptitle('Per-gene count distributions (sample)', fontsize=10)
    plt.tight_layout()

    return {'overdispersed_fraction': overdispersed_frac, 'n_genes': counts.shape[1]}


def validate_correlation_structure(sigma, observed_counts, n_genes_sample=500, seed=0):
    """
    Check that the simulated data preserves the designed correlation structure.

    Computes the Pearson correlation between the upper-triangle entries of the
    true sigma and the empirical correlation matrix from observed_counts (on a
    random gene subsample for speed).

    Returns
    -------
    dict with:
      - sigma_empirical_corr  : Pearson r between true and empirical correlations
      - mean_abs_sigma        : mean |sigma| in the ground-truth matrix
      - mean_abs_empirical    : mean |r| in the empirical matrix
    """
    rng = np.random.default_rng(seed)
    n = min(n_genes_sample, sigma.shape[0], observed_counts.shape[1])
    idx = rng.choice(sigma.shape[0], size=n, replace=False)

    sigma_sub = sigma[np.ix_(idx, idx)]
    counts_sub = observed_counts[:, idx].astype(float)

    # Normalize columns for empirical correlation
    col_std = counts_sub.std(axis=0)
    col_std[col_std == 0] = 1
    counts_z = (counts_sub - counts_sub.mean(axis=0)) / col_std
    emp_corr = (counts_z.T @ counts_z) / counts_sub.shape[0]

    tri_idx = np.triu_indices(n, k=1)
    true_vals = sigma_sub[tri_idx]
    emp_vals = emp_corr[tri_idx]

    corr, _ = pearsonr(true_vals, emp_vals)

    return {
        'sigma_empirical_corr': corr,
        'mean_abs_sigma': np.mean(np.abs(true_vals)),
        'mean_abs_empirical': np.mean(np.abs(emp_vals)),
    }


def run_validation(true_counts, observed_counts, sigma):
    """Run all validation checks and print a summary report."""
    print('=== Simulation Validation ===\n')

    res_exp = validate_gene_expression_distribution(true_counts)
    print(f'[Expression distribution]')
    print(f'  genes with non-zero expression : {res_exp["n_genes"]}')
    print(f'  KS stat vs exponential         : {res_exp["ks_stat"]:.3f}  (p={res_exp["ks_p"]:.2e})\n')

    res_umi = validate_cell_expression_distribution(true_counts)
    print(f'[Expression distribution]')
    print(f'  genes with non-zero expression : {res_umi["n_cells"]}')
    print(f'  KS stat vs exponential         : {res_umi["ks_stat"]:.3f}  (p={res_umi["ks_p"]:.2e})\n')

    res_do = validate_dropout_rate(true_counts, observed_counts)
    print(f'[Dropout]')
    print(f'  observed zero fraction         : {res_do["observed_zero_fraction"]:.3f}')
    print(f'  mean detection rate            : {res_do["mean_detection_rate"]:.3f}')
    print(f'  expression-detection corr      : {res_do["expression_detection_correlation"]:.3f}\n')

    res_nb = validate_count_distribution(true_counts)
    print(f'[Count distributions]')
    print(f'  fraction overdispersed genes   : {res_nb["overdispersed_fraction"]:.3f}\n')

    res_corr = validate_correlation_structure(sigma, observed_counts)
    print(f'[Correlation structure]')
    print(f'  sigma vs empirical corr (r)    : {res_corr["sigma_empirical_corr"]:.3f}')
    print(f'  mean |sigma|                   : {res_corr["mean_abs_sigma"]:.4f}')
    print(f'  mean |empirical|               : {res_corr["mean_abs_empirical"]:.4f}')

    return {**res_exp, **res_do, **res_nb, **res_corr}


# ── GMP-Cor calibration sweep ───────────────────────────────────────────────

def rho_sweep(rho_values=None, n_cells=500, n_genes=500, dropout_rate=1,
              shape=1.5, hub_probability=0.2,
              n_repeats=3, seed=42, log_file=None):
    """
    Sweep correlation strength (rho) and compute GMP-Cor at each level.

    GMP-Cor = Σ max(λᵢ − λ*_scrambled, 0) for all eigenvalues λᵢ,
    where λ*_scrambled is the maximum eigenvalue of the scrambled matrix.

    Parameters
    ----------
    rho_values      : list of rho values to sweep (default [0.1, 0.3, 0.5, 0.7, 0.9])
    n_cells         : number of cells per simulation
    n_genes         : number of genes per simulation
    dropout_rate    : dropout rate passed to simulate_scRNA_data
    shape           : Pareto shape for cluster-size distribution in generate_gram_hub_matrix
    hub_probability : hub connectivity probability in generate_gram_hub_matrix
    n_repeats       : independent repeats per rho value (for mean ± SD)
    seed            : base RNG seed; repeat k uses seed + k
    log_file        : path to write a JSON parameter log; prints to stdout if None

    Returns
    -------
    pd.DataFrame with columns: rho, repeat, gmp_cor
    """
    import json
    import datetime
    from .analysis_functions import get_eig_dist

    if rho_values is None:
        rho_values = [0.1, 0.3, 0.5, 0.7, 0.9]

    params = dict(
        rho_values=rho_values,
        n_cells=n_cells,
        n_genes=n_genes,
        dropout_rate=dropout_rate,
        shape=shape,
        hub_probability=hub_probability,
        n_repeats=n_repeats,
        sigma_seed=seed,
        count_seeds=list(range(seed, seed + n_repeats)),
        note='sigma matrix fixed per rho (sigma_seed); only count sampling varies across repeats (count_seeds)',
        timestamp=datetime.datetime.now().isoformat(),
    )
    log_str = json.dumps(params, indent=2)
    if log_file:
        with open(log_file, 'w') as f:
            f.write(log_str + '\n')
    else:
        print('rho_sweep parameters:\n' + log_str)

    records = []
    for rho in rho_values:
        # Generate sigma once per rho so cluster structure is identical across repeats;
        # only the count-sampling RNG varies, isolating rho as the sole variable.
        sigma = generate_gram_hub_matrix(n_genes, rho, shape, hub_probability, seed=seed)
        for rep in range(n_repeats):
            _, observed = simulate_scRNA_data(
                n_cells=n_cells, n_genes=n_genes, sigma=sigma,
                dropout_rate=dropout_rate, seed=seed + rep,
            )
            pcs, pcs1, _ = get_eig_dist(observed, norm=True, log=False, norm_sum=100)
            threshold = pcs1.max()
            gmp_cor = float(np.sum(np.maximum(pcs - threshold, 0)))
            records.append({'rho': rho, 'repeat': rep, 'gmp_cor': gmp_cor})
            print(f'  rho={rho:.2f}  rep={rep}  GMP-Cor={gmp_cor:.3f}')

    return pd.DataFrame(records)


# ── Subpopulation mixing scenario ───────────────────────────────────────────

def subpopulation_mixing(
    n_cells=1000,
    n_genes=2000,
    mixing_ratio=0.5,
    rho_low=0.1,
    rho_high=0.8,
    seed_a=20,
    seed_b=21,
    count_seed_a=0,
    count_seed_b=1,
    dropout_rate=1.0,
    shape=1.5,
    hub_probability=0.2,
    log_file=None,
):
    """
    Simulate two contrasting subpopulation-mixing conditions and compute GMP-Cor.

    Two conditions are run:
      - **Dysregulated**: sub-pop A (rho_low, seed_a) + sub-pop B (rho_low, seed_b)
      - **Regulated**:    sub-pop A (rho_high, seed_a) + sub-pop B (rho_high, seed_b)

    Using different sigma seeds (seed_a vs seed_b) gives each sub-population a
    distinct hub-network topology even when rho is identical, matching the
    reviewer scenario of two internally-regulated but transcriptomically distinct
    cell populations.

    GMP-Cor = Σ max(λᵢ − λ*_scrambled, 0) for all eigenvalues λᵢ.

    Parameters
    ----------
    n_cells         : total cells (split mixing_ratio / (1-mixing_ratio) between subs)
    n_genes         : number of genes
    mixing_ratio    : fraction of cells assigned to sub-pop A (0.5 = 50/50)
    rho_low         : shared-variance fraction for the dysregulated condition
    rho_high        : shared-variance fraction for the regulated condition
    seed_a          : sigma-matrix RNG seed for sub-population A
    seed_b          : sigma-matrix RNG seed for sub-population B (distinct network)
    count_seed_a    : count-generation RNG seed for sub-population A
    count_seed_b    : count-generation RNG seed for sub-population B
    dropout_rate    : dropout-rate parameter passed to simulate_scRNA_data
    shape           : Pareto shape for cluster-size distribution
    hub_probability : probability of cluster hub connecting to global hub
    log_file        : path for JSON results log; prints to stdout if None

    Returns
    -------
    dict with keys 'params', 'dysregulated', 'regulated'.
    Each condition dict contains:
      rho, n_cells_a, n_cells_b, n_cells_total,
      gmp_cor_subpop_a, gmp_cor_subpop_b, gmp_cor_combined
    """
    import json
    import datetime
    from .analysis_functions import get_eig_dist

    n_cells_a = int(round(n_cells * mixing_ratio))
    n_cells_b = n_cells - n_cells_a

    params = dict(
        n_cells=n_cells,
        n_cells_a=n_cells_a,
        n_cells_b=n_cells_b,
        n_genes=n_genes,
        mixing_ratio=mixing_ratio,
        rho_low=rho_low,
        rho_high=rho_high,
        seed_a=seed_a,
        seed_b=seed_b,
        count_seed_a=count_seed_a,
        count_seed_b=count_seed_b,
        dropout_rate=dropout_rate,
        shape=shape,
        hub_probability=hub_probability,
        gmp_cor_definition='sum(max(lambda_i - max_scrambled_lambda, 0))',
        note=(
            'seed_a / seed_b determine hub-network topology (sigma matrix); '
            'count_seed_a / count_seed_b determine count-sampling noise. '
            'Two distinct networks per rho level — different seeds → different hub structures.'
        ),
        timestamp=datetime.datetime.now().isoformat(),
    )

    def _gmp_cor(observed):
        pcs, pcs1, _ = get_eig_dist(observed, norm=True, log=False, norm_sum=100)
        return float(np.sum(np.maximum(pcs - pcs1.max(), 0)))

    results = {'params': params}

    for condition, rho in [('dysregulated', rho_low), ('regulated', rho_high)]:
        # Different hub-network topologies via seed_a vs seed_b; same rho within condition
        sigma_a = generate_gram_hub_matrix(n_genes, rho, shape, hub_probability, seed=seed_a)
        sigma_b = generate_gram_hub_matrix(n_genes, rho, shape, hub_probability, seed=seed_b)

        _, obs_a = simulate_scRNA_data(
            n_cells=n_cells_a, n_genes=n_genes, sigma=sigma_a,
            dropout_rate=dropout_rate, seed=count_seed_a,
        )
        _, obs_b = simulate_scRNA_data(
            n_cells=n_cells_b, n_genes=n_genes, sigma=sigma_b,
            dropout_rate=dropout_rate, seed=count_seed_b,
        )

        obs_combined = np.vstack([obs_a, obs_b])

        print(f'\n[{condition.upper()}]  rho={rho}')
        print(f'  Computing GMP-Cor for sub-pop A  ({n_cells_a} cells) ...')
        gmp_a = _gmp_cor(obs_a)
        print(f'  Computing GMP-Cor for sub-pop B  ({n_cells_b} cells) ...')
        gmp_b = _gmp_cor(obs_b)
        print(f'  Computing GMP-Cor for combined   ({n_cells} cells) ...')
        gmp_combined = _gmp_cor(obs_combined)

        results[condition] = {
            'rho': rho,
            'n_cells_a': n_cells_a,
            'n_cells_b': n_cells_b,
            'n_cells_total': n_cells,
            'gmp_cor_subpop_a': gmp_a,
            'gmp_cor_subpop_b': gmp_b,
            'gmp_cor_combined': gmp_combined,
        }

        print(f'  => sub-pop A: {gmp_a:.4f} | sub-pop B: {gmp_b:.4f} | combined: {gmp_combined:.4f}')

    log_str = json.dumps(results, indent=2)
    if log_file:
        with open(log_file, 'w') as fh:
            fh.write(log_str + '\n')
        print(f'\nJSON log written to: {log_file}')
    else:
        print('\nsubpopulation_mixing results:\n' + log_str)

    return results


# ── Graphical-lasso covariance entropy ───────────────────────────────────────

def calculate_entropy(DE_timepoints, reg_param):
    """Fit a GraphicalLasso model and return the log-determinant of its covariance."""
    model = GraphicalLasso(alpha=reg_param, max_iter=1000, tol=0.001)
    model.fit(DE_timepoints)
    sign, entropy = np.linalg.slogdet(model.covariance_)

    g = sns.clustermap(
        model.covariance_,
        cmap='viridis',
        annot=False,
        figsize=(6, 6),
        xticklabels=False,
        yticklabels=False,
        linewidths=0,
        cbar_kws={'label': 'Covariance'},
    )
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    g.cax.set_position([1.02, 0.3, 0.03, 0.4])
    g.figure.suptitle(fr'GLASSO regularized Covariance Matrix $\rho={reg_param}$', y=0.85)
    plt.show()
    return entropy


# ── Differential-expression helpers ─────────────────────────────────────────

def run_de_analysis(data, metadata, condition_col='condition', test_time=1):
    """
    Differential expression via Welch t-test with BH correction.

    Parameters
    ----------
    data         : pd.DataFrame  rows=genes, cols=samples
    metadata     : pd.DataFrame  rows=samples, must contain condition_col
    condition_col: column in metadata identifying the experimental condition
    test_time    : index into unique conditions identifying the test group
    """
    results = []
    genes = data.index
    log_data = np.log2(data + 1e-9)

    groups = metadata[condition_col].unique()
    cond1 = metadata[metadata[condition_col] == groups[0]].index
    cond2 = metadata[metadata[condition_col] == groups[test_time]].index

    for gene in genes:
        sample1 = log_data.loc[gene, cond1]
        sample2 = log_data.loc[gene, cond2]
        l2fc = sample2.mean() - sample1.mean()
        t_stat, p_val = scipy_stats.ttest_ind(sample2, sample1, equal_var=False)
        results.append({'gene': gene, 'log2FoldChange': l2fc, 'pvalue': p_val, 'stat': t_stat})

    res_df = pd.DataFrame(results).set_index('gene')
    res_df['padj'] = multipletests(res_df['pvalue'], method='fdr_bh')[1]
    return res_df.sort_values('padj')


def plot_volcano(res_df, lfc_thresh=1.0, p_thresh=0.05):
    """Volcano plot from a DE results DataFrame produced by run_de_analysis."""
    plot_df = res_df.copy()
    plot_df['-log10p'] = -np.log10(plot_df['padj'])
    plot_df['group'] = 'Not Significant'
    plot_df.loc[(plot_df['log2FoldChange'] > lfc_thresh) & (plot_df['padj'] < p_thresh), 'group'] = 'Up-regulated'
    plot_df.loc[(plot_df['log2FoldChange'] < -lfc_thresh) & (plot_df['padj'] < p_thresh), 'group'] = 'Down-regulated'

    colors = {'Not Significant': 'grey', 'Up-regulated': '#d62728', 'Down-regulated': '#1f77b4'}
    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=plot_df, x='log2FoldChange', y='-log10p',
                    hue='group', palette=colors, alpha=0.7, edgecolor=None)
    plt.axhline(-np.log10(p_thresh), color='black', linestyle='--', alpha=0.5)
    plt.axvline(lfc_thresh, color='black', linestyle='--', alpha=0.5)
    plt.axvline(-lfc_thresh, color='black', linestyle='--', alpha=0.5)
    plt.title('Differential Expression Volcano Plot', fontsize=15)
    plt.xlabel(r'$\log_2$ Fold Change', fontsize=12)
    plt.ylabel(r'$-\log_{10}$ Adjusted P-value', fontsize=12)
    plt.legend(title='Status', loc='upper right')
    plt.xlim([-8, 4])
    plt.show()