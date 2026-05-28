import numpy as np
import pandas as pd
from numba.np.arraymath import np_average
from numpy import linalg as la
from sklearn.decomposition import SparsePCA
from scipy.stats import rankdata

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

def spearman_ranking(m):
    m = np.array([rankdata(row) for row in m.T]).T
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


def get_sparse_pcs(m, n_components, alpha=0.5):
    # Get the principal components of a matrix m using sparse PCA
    spca = SparsePCA(n_components=n_components, alpha=alpha)
    spca.fit(m)
    return spca.components_


def get_eig_dist(m, norm=True, log=False, norm_method='sum', norm_sum=1):
    # get the eigenvalue distribution of the normalized matrix m
    # scramble the matrix m, and get the eigenvalue distribution of the normalized matrix
    #m = log_transform(m)  # z-transform the matrix m
    rep=10
    # remove zero genes:
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

