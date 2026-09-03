"""
Synthetic scRNA-seq data generation and the metrics computed on it, used to
validate the GMP-Cor eigenvalue-spectrum statistic (see CLAUDE.md / README.md)
against ground truth that the real experimental data cannot provide.

Two things live here:

  - A generative model of bacterial scRNA-seq counts (`generate_gram_hub_matrix`
    + `simulate_scRNA_data`) with a known gene-gene correlation structure
    (`sigma`), realistic negative-binomial marginals, library-size variation and
    expression-dependent dropout. `validate_*` / `run_validation` check that the
    simulated data reproduce the qualitative statistical properties (heavy-tailed
    expression, overdispersion, dropout-expression coupling) seen in real data.

  - Experiment runners built on that model: `rho_sweep` (does GMP-Cor track the
    known correlation strength rho?), `subpopulation_mixing` and
    `inverted_subpopulation_mixing` (can mixing two internally-regulated but
    distinct sub-populations spuriously inflate GMP-Cor, mimicking loss of
    regulation?). These were written for the reviewer-response simulations under
    `simulations/` and write their parameters/results to JSON logs for
    reproducibility (see `results/simulation_results/` in CLAUDE.md).

Not part of the eigenvalue-spectrum work: `calculate_entropy` (graphical-lasso
covariance entropy) and `run_de_analysis` / `plot_volcano` (a simple Welch
t-test differential-expression pipeline) are exploratory helpers used from
`simulations/notebook.ipynb`, independent of the rest of the module.

This module is imported as `src.simulations` (or `af`/`sf`-style aliases) from
notebooks and the standalone scripts under `simulations/`; it is not run
directly.
"""
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

def simulate_scRNA_data(n_cells=1000, n_genes=2000, sigma=None, rho=0.9, dropout_rate=2, inv_gamma_shape=1.5, inv_gamma_scale=0.01, seed=None, gene_mu=None):
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
    gene_mu      : optional length-n_genes array of per-gene NB means. If None (default)
                   the means are drawn from the inverse-Gamma prior exactly as before.
                   Supplying them explicitly lets two simulated populations share, or
                   deliberately differ in, their marginal expression profile — see
                   `draw_gene_means` / `invert_gene_means`.

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
    if gene_mu is not None:
        gene_mu = np.asarray(gene_mu, dtype=float).ravel()
        if gene_mu.size != n_genes:
            raise ValueError(f'gene_mu has {gene_mu.size} entries, expected {n_genes}')
    for i in range(n_genes):
        # Gene mean drawn from inverse-Gamma (heavy-tailed expression distribution),
        # unless an explicit expression profile was supplied.
        mu_i = gene_mu[i] if gene_mu is not None else 1.0 / rng.gamma(inv_gamma_shape, 1.0 / inv_gamma_scale, 1)
        gene_r = 0.5  # NB dispersion typical of scRNA-seq
        p_param = gene_r / (gene_r + mu_i)
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
    Compare per-cell total expression CCDF against an exponential-tail null.

    Plots the empirical CCDF of cell totals (sum over genes, i.e. per-cell
    library size), then returns the KS statistic against an exponential null
    (a rough check that the tail is heavier than exponential).
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
    Compare per-gene total expression CCDF against an exponential-tail null.

    Plots the empirical CCDF of gene totals (sum over cells), then returns the
    KS statistic against an exponential null (a rough check that the tail is
    heavier than exponential, consistent with the inverse-Gamma mean prior used
    by simulate_scRNA_data). `n_sample` is accepted for API symmetry with other
    validate_* helpers but is currently unused - the full gene set is plotted.
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


# ── Distinct (rank-inverted) sub-populations ─────────────────────────────────

def draw_gene_means(n_genes, seed=0, inv_gamma_shape=1.5, inv_gamma_scale=0.01):
    """
    Draw a per-gene negative-binomial mean profile from the inverse-Gamma prior
    used by `simulate_scRNA_data`.

    Returned as an explicit array so that two simulated sub-populations can be
    given deliberately related (e.g. rank-inverted) expression profiles instead
    of two independent draws.
    """
    rng = np.random.default_rng(seed)
    return 1.0 / rng.gamma(inv_gamma_shape, 1.0 / inv_gamma_scale, n_genes)


def invert_gene_means(gene_mu):
    """
    Reverse the expression ranking of a gene-mean profile.

    The most lowly-expressed gene receives the highest mean, the second-lowest
    receives the second-highest, and so on. The multiset of means — and hence the
    marginal expression distribution of the population — is exactly preserved;
    only its assignment to genes is reversed.

    This produces two populations that are transcriptomically opposite while
    remaining statistically identical in every global property (dynamic range,
    sparsity, library-size distribution), so any separation between them is a
    genuine difference in *which* genes are expressed, not a depth artefact.
    """
    gene_mu = np.asarray(gene_mu, dtype=float).ravel()
    order = np.argsort(gene_mu)              # gene indices, lowest mean first
    inverted = np.empty_like(gene_mu)
    inverted[order] = gene_mu[order[::-1]]   # lowest ↔ highest
    return inverted


def gmp_cor(matrix, norm=True, norm_sum=50):
    """
    GMP-Cor = Σ max(λᵢ − λ*_scrambled, 0), computed exactly as for the
    experimental data via `analysis_functions.get_eig_dist`.

    Returns a dict. `p_kept` is the number of genes surviving the all-zero filter;
    GMP-Cor is extensive in it, so `gmp_cor_per_gene` is the form that may be
    compared across matrices of different gene dimension.
    """
    from .analysis_functions import get_eig_dist
    pcs, pcs1, frac_nz = get_eig_dist(matrix, norm=norm, log=False,
                                      norm_method='sum', norm_sum=norm_sum)
    thr = pcs1.max()
    excess = pcs[pcs > thr] - thr
    p_kept = int(pcs.size)
    return dict(gmp_cor=float(excess.sum()),
                gmp_cor_per_gene=float(excess.sum() / p_kept) if p_kept else np.nan,
                lambda_max=float(pcs.max()), lambda_max_scrambled=float(thr),
                n_modes_above_threshold=int(excess.size), p_kept=p_kept,
                fraction_non_zero=float(frac_nz))


def _preprocess(matrix, norm_sum=50):
    """Row-normalize then z-transform columns — the preprocessing inside get_eig_dist."""
    from .analysis_functions import normalize, z_transform
    return z_transform(normalize(np.asarray(matrix, dtype=float),
                                 method='sum', target_sum=norm_sum))


def cell_scores(matrix, n_components=2, norm_sum=50):
    """Cell coordinates on the leading principal components of the gene-gene structure."""
    z = _preprocess(matrix, norm_sum=norm_sum)
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    return u[:, :n_components] * s[:n_components]


def group_centered(matrix, labels, norm_sum=50):
    """
    Row-normalize, then subtract each group's own gene means.

    This removes all *between-group* structure while leaving within-cell
    coordination intact, and is the basis of the dGMP diagnostic: a GMP-Cor
    driven by a mixture collapses under this operation, genuine coordination
    does not.
    """
    from .analysis_functions import normalize
    m = normalize(np.asarray(matrix, dtype=float), method='sum', target_sum=norm_sum)
    labels = np.asarray(labels)
    out = m.copy()
    for lab in np.unique(labels):
        sel = labels == lab
        out[sel] = m[sel] - m[sel].mean(axis=0)
    return out


def separation_metrics(matrix, labels, n_components=5, norm_sum=50):
    """
    How strongly do two labelled cell groups separate, and does the separating
    direction show up as a mode of the correlation spectrum?

    The group axis is the difference between the two groups' mean gene profiles in
    the same z-transformed space in which GMP-Cor is computed. Reporting the AUC on
    that axis — rather than on PC1 — matters: when a population carries strong
    internal structure the separating direction is often PC2 or PC3, and a PC1-only
    statistic then reads as "no separation" for two perfectly separated groups.

    Returns
      auc_group_axis    : AUC separating the groups along the group axis (0.5–1)
      cohens_d          : standardized group difference on that axis
      silhouette        : silhouette of the true labels in PC1–PC2 space
      best_pc_auc/index : best-separating principal component and its AUC
      group_mode_index  : mode most aligned with the group axis (0-based)
      group_mode_alignment  : |cos| between that mode and the group axis
      group_mode_eigenvalue : its eigenvalue
      bimodality_coef   : Sarle's coefficient on the group axis (> 5/9 ⇒ bimodal)
      gmm_bic_gain      : BIC(1 component) − BIC(2) on the group axis; > 0 favours two modes
      n_genes_dz_gt_0.5 : genes whose standardized group difference exceeds 0.5 —
                          the count that determines how large the separating mode is
    """
    from scipy.stats import mannwhitneyu, skew, kurtosis
    from sklearn.metrics import silhouette_score
    from sklearn.mixture import GaussianMixture

    labels = np.asarray(labels)
    z = _preprocess(matrix, norm_sum=norm_sum)
    u, sv, vt = np.linalg.svd(z, full_matrices=False)
    eigvals = sv ** 2 / z.shape[0]
    scores = u * sv

    def _auc(a, b):
        if not a.size or not b.size:
            return np.nan
        stat = mannwhitneyu(a, b, alternative='two-sided').statistic
        return max(stat, a.size * b.size - stat) / (a.size * b.size)

    dz = z[labels == 0].mean(axis=0) - z[labels == 1].mean(axis=0)
    axis = dz / np.linalg.norm(dz) if np.linalg.norm(dz) > 0 else dz
    proj = z @ axis
    pa, pb = proj[labels == 0], proj[labels == 1]

    pooled = np.sqrt((pa.var(ddof=1) + pb.var(ddof=1)) / 2) if pa.size > 1 and pb.size > 1 else 0.0
    k = min(n_components, scores.shape[1])
    pc_aucs = [_auc(scores[labels == 0, j], scores[labels == 1, j]) for j in range(k)]
    align = np.abs(vt[:k] @ axis)

    n = proj.size
    g, kt = skew(proj), kurtosis(proj, fisher=True)
    denom = kt + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)) if n > 3 else np.nan
    x = proj.reshape(-1, 1)

    return {
        'auc_group_axis': float(_auc(pa, pb)),
        'cohens_d': float(abs(pa.mean() - pb.mean()) / pooled) if pooled > 0 else np.nan,
        'silhouette': float(silhouette_score(scores[:, :2], labels)) if np.unique(labels).size > 1 else np.nan,
        'best_pc_auc': float(np.nanmax(pc_aucs)),
        'best_pc_index': int(np.nanargmax(pc_aucs)),
        'group_mode_index': int(np.argmax(align)),
        'group_mode_alignment': float(np.max(align)),
        'group_mode_eigenvalue': float(eigvals[int(np.argmax(align))]),
        'bimodality_coef': float((g ** 2 + 1) / denom) if denom and denom > 0 else np.nan,
        'gmm_bic_gain': float(GaussianMixture(1, random_state=0).fit(x).bic(x)
                              - GaussianMixture(2, random_state=0).fit(x).bic(x)),
        'n_genes_dz_gt_0.5': int((np.abs(dz) > 0.5).sum()),
    }


def inverted_subpopulation_mixing(
    n_cells=1000,
    n_genes=2000,
    ratios=(0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0),
    rho_high=0.9,
    rho_low=0.1,
    seed_a=20,
    seed_b=21,
    mu_seed=7,
    count_seed_a=0,
    count_seed_b=1,
    dropout_rate=1.0,
    shape=1.5,
    hub_probability=0.2,
    inv_gamma_scale=0.04,
    norm_sum=50,
    sigma_mode='shared',
    example_ratio=None,
    umap_neighbors=15,
    umap_min_dist=0.1,
    umap_n_pcs=50,
    umap_seed=0,
    repeat=0,
):
    """
    Mix two *genuinely distinct* internally-regulated sub-populations and track
    what the separating mode does to GMP-Cor.

    Sub-population B is built from sub-population A's expression profile with the
    gene ranking **inverted** (`invert_gene_means`): A's most lowly-expressed
    genes are B's most highly-expressed. Each population additionally carries its
    own hub-network topology (`seed_a` vs `seed_b`) and is internally regulated at
    `rho_high`. The two are therefore transcriptomically opposite — as the two
    clusters of an experimental mixture are — while sharing an identical marginal
    expression distribution, so the separation cannot be attributed to depth.

    At every mixing ratio the function reports GMP-Cor, GMP-Cor after removing the
    between-population means (`group_centered`), and how strongly the two
    populations separate along PC1. A single dysregulated population (`rho_low`)
    of the same size is computed as the reference for genuine loss of coordination.

    Parameters
    ----------
    ratios     : fractions of the total cells taken from sub-population A. 0 and 1 are
                 the pure populations; the interior values are the mixtures.
    sigma_mode : 'shared'   — both populations use the same hub network, so they differ
                              only in which genes they express (the experimental case:
                              two states of the same organism share regulatory
                              architecture). This is the configuration that reproduces
                              the behaviour of a real mixture.
                 'distinct' — each population additionally gets its own network
                              topology. The two networks dilute one another, which
                              works against the mixture signal; reported as the
                              conservative case.
    inv_gamma_scale : scale of the inverse-Gamma expression prior. The default 0.04
                 (rather than the 0.01 used elsewhere in this module) is calibrated so
                 that the simulated cells detect ~85 genes each, matching the
                 experimental matrices; at 0.01 the data are so sparse that only a
                 handful of genes carry any between-population difference and no
                 separating mode can form.
    repeat     : offset added to every count seed, for independent realizations.

    Every ratio is evaluated at the same total cell number, so the pure populations
    and the mixtures share n — without which the scrambled threshold (a pure
    function of matrix shape) differs between them and the comparison is invalid.

    example_ratio : if given, the record for that ratio additionally returns a UMAP
                 embedding, principal-component coordinates, group-axis projection and
                 leading eigenvalues of that one mixture, for plotting. The UMAP is run
                 on the top `umap_n_pcs` principal components of the same preprocessed
                 matrix GMP-Cor is computed from, matching `data_functions.get_umap`.

    Returns
    -------
    dict with 'params', 'gene_means', 'reference' (the dysregulated population),
    'ratios' (one record per mixing ratio) and 'example'.
    """
    ratios = list(ratios)

    if sigma_mode not in ('shared', 'distinct'):
        raise ValueError("sigma_mode must be 'shared' or 'distinct'")
    sigma_a = generate_gram_hub_matrix(n_genes, rho_high, shape, hub_probability, seed=seed_a)
    sigma_b = (sigma_a if sigma_mode == 'shared' else
               generate_gram_hub_matrix(n_genes, rho_high, shape, hub_probability, seed=seed_b))
    sigma_dys = generate_gram_hub_matrix(n_genes, rho_low, shape, hub_probability, seed=seed_a)

    mu_a = draw_gene_means(n_genes, seed=mu_seed, inv_gamma_scale=inv_gamma_scale)
    mu_b = invert_gene_means(mu_a)

    # Cell pools, generated once and subsampled at every ratio.
    _, pool_a = simulate_scRNA_data(n_cells=n_cells, n_genes=n_genes, sigma=sigma_a,
                                    dropout_rate=dropout_rate, seed=count_seed_a + 100 * repeat,
                                    gene_mu=mu_a)
    _, pool_b = simulate_scRNA_data(n_cells=n_cells, n_genes=n_genes, sigma=sigma_b,
                                    dropout_rate=dropout_rate, seed=count_seed_b + 100 * repeat,
                                    gene_mu=mu_b)
    _, pool_dys = simulate_scRNA_data(n_cells=n_cells, n_genes=n_genes, sigma=sigma_dys,
                                      dropout_rate=dropout_rate, seed=count_seed_a + 100 * repeat,
                                      gene_mu=mu_a)

    rng = np.random.default_rng(1000 + repeat)
    idx_a = rng.permutation(n_cells)
    idx_b = rng.permutation(n_cells)

    records = []
    example = None
    for r in ratios:
        n_a = int(round(n_cells * r))
        n_b = n_cells - n_a
        parts, labels = [], []
        if n_a:
            parts.append(pool_a[idx_a[:n_a]])
            labels.append(np.zeros(n_a, dtype=int))
        if n_b:
            parts.append(pool_b[idx_b[:n_b]])
            labels.append(np.ones(n_b, dtype=int))
        mixture = np.vstack(parts)
        lab = np.concatenate(labels)

        raw = gmp_cor(mixture, norm=True, norm_sum=norm_sum)
        cen = gmp_cor(group_centered(mixture, lab, norm_sum=norm_sum), norm=False)
        rec = dict(ratio_a=r, n_cells_a=n_a, n_cells_b=n_b, **raw)
        rec['gmp_cor_group_centered'] = cen['gmp_cor']
        rec['gmp_cor_group_centered_per_gene'] = cen['gmp_cor_per_gene']
        rec['d_gmp'] = float(1 - cen['gmp_cor'] / raw['gmp_cor']) if raw['gmp_cor'] > 0 else np.nan
        if n_a and n_b:
            rec.update(separation_metrics(mixture, lab, norm_sum=norm_sum))
        records.append(rec)

        if example_ratio is not None and n_a and n_b and abs(r - example_ratio) < 1e-9:
            import umap as _umap
            z = _preprocess(mixture, norm_sum=norm_sum)
            u, sv, vt = np.linalg.svd(z, full_matrices=False)
            dz = z[lab == 0].mean(axis=0) - z[lab == 1].mean(axis=0)
            axis = dz / np.linalg.norm(dz)
            # UMAP on the top principal components, as in data_functions.get_umap
            n_pcs = min(umap_n_pcs, z.shape[1])
            embedding = _umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist,
                                   n_components=2, random_state=umap_seed
                                   ).fit_transform(z @ vt[:n_pcs].T)
            example = dict(ratio_a=r, labels=lab.tolist(),
                           umap=np.asarray(embedding).tolist(),
                           umap_params=dict(n_neighbors=umap_neighbors,
                                            min_dist=umap_min_dist, n_pcs=n_pcs,
                                            random_state=umap_seed),
                           pc_scores=(u[:, :2] * sv[:2]).tolist(),
                           group_axis_projection=(z @ axis).tolist(),
                           eigenvalues=(sv[:20] ** 2 / z.shape[0]).tolist())

    ref = gmp_cor(pool_dys, norm=True, norm_sum=norm_sum)

    return dict(
        params=dict(n_cells=n_cells, n_genes=n_genes, ratios=ratios,
                    rho_high=rho_high, rho_low=rho_low, seed_a=seed_a, seed_b=seed_b,
                    mu_seed=mu_seed, count_seed_a=count_seed_a, count_seed_b=count_seed_b,
                    dropout_rate=dropout_rate, shape=shape,
                    hub_probability=hub_probability, inv_gamma_scale=inv_gamma_scale,
                    norm_sum=norm_sum,
                    sigma_mode=sigma_mode, repeat=repeat,
                    gmp_cor_definition='sum(max(lambda_i - max_scrambled_lambda, 0))'),
        gene_means=dict(mu_a=mu_a.tolist(), mu_b=mu_b.tolist(),
                        spearman_a_vs_b=float(scipy_stats.spearmanr(mu_a, mu_b).statistic)),
        reference=dict(rho=rho_low, **ref),
        ratios=records,
        example=example,
    )

# ── Graphical-lasso covariance entropy ───────────────────────────────────────

def calculate_entropy(DE_timepoints, reg_param):
    """Fit a GraphicalLasso model and return the log-determinant of its covariance."""
    model = GraphicalLasso(alpha=reg_param, max_iter=1000, tol=0.001)
    model.fit(DE_timepoints)
    # log-determinant of the regularized covariance as an entropy proxy (up to an
    # additive constant); `sign` is discarded since a valid covariance estimate is
    # always positive semi-definite (sign should be +1).
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
    # tiny pseudocount avoids log2(0) = -inf for genes with zero counts in a sample
    log_data = np.log2(data + 1e-9)

    groups = metadata[condition_col].unique()
    # groups[0] is always the reference; test_time selects which other condition to
    # contrast it against. Assumes at least `test_time + 1` distinct condition values
    # are present in metadata[condition_col] - raises IndexError otherwise.
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