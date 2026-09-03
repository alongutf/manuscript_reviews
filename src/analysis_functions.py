"""
Pure numerical routines for the gene-gene correlation / eigenvalue-spectrum
analysis described in CLAUDE.md: matrix preprocessing (normalization,
log/z-transforms, column scrambling), principal-component / eigenvalue
computation, the empirical-vs-scrambled eigenvalue spectrum comparison that
GMP-Cor is built on (get_eig_dist), the analytic Marchenko-Pastur null density
(mp_distribution) it is compared against, and a handful of per-mode diagnostics
(coordination_score, participation_ratio, eigenvector_entropy) used to
interpret which genes drive a given eigenmode. Everything here operates on
plain numpy arrays (no I/O) and is imported as `src.analysis_functions`
(conventionally aliased `af`) by data_functions.py, simulations.py and the
scripts/ notebooks and figure scripts - it is not run directly.
"""
import numpy as np
import pandas as pd
from numba.np.arraymath import np_average
from numpy import linalg as la
from sklearn.decomposition import SparsePCA
from scipy.stats import rankdata

def scramble(m):
    # Independently permute the *cell* (row) values within each gene (column),
    # i.e. destroy any gene-gene correlation while preserving each gene's own
    # marginal distribution across cells. This is the "scrambled" null model
    # get_eig_dist compares the real eigenvalue spectrum against: a matrix with
    # the same per-gene marginals but no cross-gene structure.
    #m = np.array([np.random.permutation(row) for row in m])
    # Scramble the row indices in each column of a matrix m
    m = np.array([np.random.permutation(row) for row in m.T]).T
    return m


def normalize(m, method='norm', target_sum=1):
    # Rescale each row (cell) of m to a common total: method='sum' divides by the
    # row sum (library-size / depth normalization, the usual scRNA-seq choice),
    # any other method value divides by the row's Euclidean norm instead. Both
    # are then scaled to target_sum. Rows that summed/normed to zero produce a
    # 0/0 division -> NaN, which is reset to 0 (an all-zero cell stays all-zero
    # rather than propagating NaN into downstream computations).
    if method == 'sum':
        m = target_sum*m / m.sum(axis=1)[:, None]
    else:
        m = target_sum*m / la.norm(m, axis=1)[:, None]
    m[np.isnan(m)] = 0
    return m

def spearman_ranking(m):
    # Replace each gene's (column's) values across cells with their ranks
    # (ties handled by scipy.stats.rankdata's default averaging), i.e. convert m
    # to a Spearman-rank representation before correlation. Not called by
    # get_eig_dist itself - available for callers who want a rank-based rather
    # than Pearson-style (z-transform-based) eigenvalue spectrum.
    m = np.array([rankdata(row) for row in m.T]).T
    return m


def z_transform(m):
    # Standardize each gene (column) to zero mean, unit variance across cells -
    # this is what makes the subsequent SVD/eigenvalues equivalent to an
    # eigendecomposition of the gene-gene Pearson correlation matrix rather than
    # the raw covariance. Constant (zero-variance) columns give 0/0 -> NaN,
    # reset to 0 so they contribute no variance rather than corrupting the SVD.
    m = (m - m.mean(axis=0)) / m.std(axis=0)
    m[np.isnan(m)] = 0
    return m


def log_transform(m):
    # log1p-style transform (log(m + 1), so log(0) = 0 rather than -inf) to
    # reduce the influence of a small number of very highly expressed genes.
    return np.log(m + 1)


def get_pcs(m):
    # Eigenvalues of the gene-gene correlation/covariance structure of m
    # (n_cells x p_genes), via the singular values of m itself rather than
    # forming the p x p covariance matrix explicitly. Returned as a length-p
    # array (eigenvalue_i = singular_value_i**2 / n_cells) padded with zeros
    # for the min(n, p) .. p-1 structurally-zero directions when p > n (more
    # genes than cells, the usual regime for this data) - la.svd(m) with the
    # default full_matrices=True still returns only min(n, p) singular values,
    # so this padding keeps the returned array's length equal to p regardless.
    n = m.shape[0]
    p = m.shape[1]
    pcs = np.zeros(p)
    pcs[:min(n, p)] = la.svd(m)[1]**2/n
    return pcs


def get_sparse_pcs(m, n_components, alpha=0.5):
    # Sparse-PCA gene loadings for m (n_cells x p_genes): n_components components,
    # each a p-length loading vector with an L1 penalty (alpha) encouraging most
    # entries to be exactly zero, i.e. each component driven by relatively few
    # genes. Returned as (n_components, p_genes); not used by get_eig_dist, which
    # uses the dense SVD (get_pcs) instead.
    spca = SparsePCA(n_components=n_components, alpha=alpha)
    spca.fit(m)
    return spca.components_


def get_eig_dist(m, norm=True, log=False, norm_method='sum', norm_sum=1):
    # Core GMP-Cor computation: return the empirical eigenvalue spectrum of m's
    # gene-gene correlation structure alongside the null spectrum obtained by
    # scrambling away gene-gene correlations, so the two can be compared (the
    # excess of empirical eigenvalues above the scrambled maximum is GMP-Cor -
    # see src/simulations.py's gmp_cor() and analysis notebooks).
    #
    # m : (n_cells x p_genes) raw count matrix.
    # norm / norm_method / norm_sum : passed to normalize() - row (per-cell)
    #     normalization applied before z-transforming, if norm is True.
    # log : whether to log_transform after normalizing and before z-transforming.
    #
    # Returns
    #   pcs             : (p_kept,) empirical eigenvalues of the (filtered,
    #                      preprocessed) matrix, see get_pcs.
    #   pcs1            : (p_kept,) eigenvalues of the scrambled matrix, averaged
    #                      over `rep` independent scrambles to reduce noise in the
    #                      null spectrum (and in particular in its maximum, the
    #                      GMP-Cor threshold).
    #   fraction_non_zero: scalar, fraction of non-zero entries in m after the
    #                      zero-gene/zero-cell filter below (a sparsity diagnostic).
    #m = log_transform(m)  # z-transform the matrix m
    rep=10
    # remove zero genes:
    # drop genes with no detected counts in any cell, and cells with no detected
    # genes - both would otherwise contribute rows/columns of all-zero variance
    gene_sums = (m>0).sum(axis=0)
    min_cells = 1
    m = m[:,gene_sums >= min_cells]
    min_genes = 1
    cell_sums = (m>0).sum(axis=1)
    m = m[cell_sums >= min_genes,:]
    fraction_non_zero = np.sum(m.flatten()>0)/len(m.flatten())
    if norm:
        m = normalize(m, method=norm_method, target_sum=norm_sum)  # normalize the rows of the matrix m
    #   m1 = normalize(m1, method=norm_method, target_sum=norm_sum)
    if log:
        m = log_transform(m)  # log-transform the matrix m
    #    m1 = log_transform(m1)
    m = z_transform(m)
    pcs1 = np.zeros(m.shape[1])
    for _ in range(rep):
        m1 = m.copy()  # copy the matrix m for scrambling

        m1 = scramble(m1)  # scramble the matrix m
#        m1 = z_transform(m1)  # z-transform the matrix m1
        pcs1 += get_pcs(m1)  # get the principal components of the matrix m1
    pcs1 = pcs1 / rep
    pcs = get_pcs(m)  # get the principal components of the matrix m
    print(m.shape)


    return pcs, pcs1, fraction_non_zero


def get_eig_vectors(m, n_top=5, norm=True, log=False, norm_method='sum', norm_sum=1):
    # Compute the leading eigenvectors (gene loadings) of the gene-gene correlation
    # structure of m (cells x genes), together with their eigenvalues and the
    # scrambled noise threshold (lambda_max^scr) used by GMP-Cor.
    #
    # Preprocessing mirrors get_eig_dist exactly so the spectrum is identical:
    # drop all-zero genes/cells, normalize rows, optional log, z-transform columns.
    #
    # Returns
    #   eigvals   : (n_top,)            leading eigenvalues (singular_value**2 / n_cells)
    #   eigvecs   : (n_top, n_kept)     leading right singular vectors = gene loadings
    #   threshold : float              mean over reps of the max scrambled eigenvalue
    #   kept_cols : (n_genes,) bool    mask of original columns of m that survived filtering
    #
    # If n_top is None, all modes are returned.
    rep = 10
    # remove zero genes / cells (same rule as get_eig_dist)
    gene_sums = (m > 0).sum(axis=0)
    kept_cols = gene_sums >= 1
    m = m[:, kept_cols]
    cell_sums = (m > 0).sum(axis=1)
    m = m[cell_sums >= 1, :]
    if norm:
        m = normalize(m, method=norm_method, target_sum=norm_sum)
    if log:
        m = log_transform(m)
    m = z_transform(m)
    n = m.shape[0]
    # scrambled noise threshold: mean of the max eigenvalue over rep scrambles
    max_scr = 0.0
    for _ in range(rep):
        m1 = scramble(m.copy())
        max_scr += get_pcs(m1).max()
    threshold = max_scr / rep
    # SVD: rows of Vt are the gene-space eigenvectors of the correlation matrix
    _, s, vt = la.svd(m, full_matrices=False)
    eigvals = s ** 2 / n
    if n_top is None:
        n_top = vt.shape[0]
    else:
        n_top = min(n_top, vt.shape[0])
    return eigvals[:n_top], vt[:n_top], threshold, kept_cols


def coordination_score(eigvals, eigvecs, threshold):
    # Per-gene coordination (leverage) score: the diagonal of the de-noised signal
    # covariance reconstructed from the modes above the scrambled threshold,
    #   score_i = sum_{k: eig_k > threshold} (eig_k - threshold) * eigvecs[k, i]**2 .
    # This is rotation-invariant within the (possibly near-degenerate) signal subspace,
    # unlike any individual eigenvector, so it gives a stable per-gene readout of how
    # much each gene participates in the coordinated structure. Returns zeros when no
    # mode exceeds the threshold (fully dysregulated sample).
    eigvals = np.asarray(eigvals)
    signal = eigvals > threshold
    score = np.zeros(eigvecs.shape[1])
    for k in np.where(signal)[0]:
        score += (eigvals[k] - threshold) * eigvecs[k] ** 2
    return score


def participation_ratio(v):
    # Participation ratio of a (unit-norm) eigenvector v: ranges from 1 (localized on a
    # single gene) to len(v) (fully delocalized). PR = (sum v^2)^2 / sum v^4.
    v = np.asarray(v, dtype=float)
    s2 = np.sum(v ** 2)
    s4 = np.sum(v ** 4)
    return (s2 ** 2) / s4 if s4 > 0 else 0.0


def eigenvector_entropy(v):
    # Shannon entropy of an eigenvector treated as a probability distribution over genes.
    # Each gene's squared loading is its "participation probability" p_i = v_i^2 / sum_j v_j^2
    # (sum_i p_i = 1; for a unit-norm eigenvector the denominator is 1). The Shannon entropy
    #   H = -sum_i p_i * ln(p_i)
    # measures how the mode's variance is spread across genes:
    #   H -> 0          a single gene carries the whole mode (fully localized),
    #   H -> ln(n)      every gene participates equally (fully delocalized).
    # We report three quantities:
    #   entropy        : H in nats.
    #   norm_entropy   : H / ln(n) in [0, 1] -- 0 = one gene, 1 = all genes participate.
    #   effective_n    : exp(H), the effective number of participating ("dominant") genes,
    #                    ranging from 1 to n (compare with participation_ratio, which is the
    #                    exp of the order-2 Renyi entropy and is dominated by the few largest
    #                    loadings; the Shannon effective_n weights moderate loadings more).
    v = np.asarray(v, dtype=float)
    p = v ** 2
    s = p.sum()
    n = len(v)
    if s <= 0 or n == 0:
        return 0.0, 0.0, 0.0
    p = p / s
    nz = p[p > 0]
    entropy = float(-np.sum(nz * np.log(nz)))
    norm_entropy = entropy / np.log(n) if n > 1 else 0.0
    effective_n = float(np.exp(entropy))
    return entropy, norm_entropy, effective_n


def mp_distribution(x, a):
    # Marchenko-Pastur distribution with ratio a
    l_min = (1-np.sqrt(a))**2
    l_max = (1+np.sqrt(a))**2
    if l_min < x < l_max:
        f = (1/(2*np.pi*x*a))*np.sqrt((x-l_min)*(l_max-x))
    else:
        f = 0
    return f

def get_entropy(pcs):
    p = len(pcs)
    P = pcs / p
    P = P[P > 0]
    return np.exp(-np.sum(P * np.log(P[P > 0])))

