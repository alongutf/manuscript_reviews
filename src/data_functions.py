"""
Defines AnnMat, the annotated cell x gene count-matrix class used throughout the
pipeline (see CLAUDE.md), plus the functions that build, filter, transform and
visualize it: get_annotated_data() converts a raw probe-count CSV into an
AnnMat; the AnnMat filtering methods (filter_by_umi_count, filter_by_gene_
dispersion, filter_by_mean_expression, remove_highly_expressed_genes,
remove_non_protein_coding_genes) mark cells/genes for exclusion without
copying data until get_filtered_matrix() materializes the result;
plot_eig_dist / pcs_to_csv / get_histogram_data turn eigenvalue spectra
(computed in analysis_functions.py) into figures and CSV histograms;
project_on_pcs / get_umap produce PCA/UMAP embeddings; concat_datasets and
equate_dims combine or size-match multiple AnnMat datasets (e.g. for the UMAP
and cross-dataset comparison notebooks). This module is imported (not run
directly) as `src.data_functions`, typically aliased `df`, from the scripts/
notebooks.
"""
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
import src.analysis_functions as af
import umap
import os


class AnnMat:
    """
    Annotated cell x gene count matrix: `m` is the raw (n_objects x n_variables)
    numeric array, `obj_names`/`var_names` are its row/column labels (cell
    barcodes and gene names for the usual scRNA-seq use, but also used
    transposed - see transpose_matrix), and `batch` is an optional per-row
    dataset/batch label (set by concat_datasets).

    Filtering is lazy: `filtered_obj` / `filtered_var` are boolean masks over
    obj_names / var_names that the filter_by_* methods narrow (AND together)
    in place, without ever copying or slicing `m`. Call get_filtered_matrix()
    to materialize a new AnnMat containing only the currently-unmasked rows
    and columns; until then `m` itself is unchanged.
    """
    def __init__(self, m, obj_names=None, var_names=None, batch=None):
        self.m = m
        self.obj_names = np.array(obj_names)
        self.var_names = np.array(var_names)
        self.filtered_obj = np.full(len(self.obj_names), True)
        self.filtered_var = np.full(len(self.var_names), True)
        self.batch = batch

    def __str__(self):
        n, p = self.m.shape
        return f"annotated matrix with {n} objects and {p} variables"

    def to_csv(self, path):
        # save the (unfiltered) matrix to csv: rows = obj_names, columns = var_names
        df = pd.DataFrame(self.m, index=self.obj_names, columns=self.var_names)
        df.to_csv(path)
        print('Data saved to', path)

    def filter_by_umi_count(self, umi_min, umi_max, target_cells=None, plot=True):
        # Mark cells (rows) whose total UMI count falls outside (umi_min, umi_max)
        # for exclusion. If target_cells is given, umi_min is instead derived so
        # that (after already-excluded top cells above umi_max) exactly
        # target_cells cells remain - i.e. keep the top target_cells cells by UMI
        # count. Updates self.filtered_obj in place; does not modify self.m.
        data_mat = self.m
        umi_count = np.sum(data_mat, axis=1)
        num_cells_above_max = np.sum(umi_count > umi_max)
        if target_cells is not None:
            umi_count_sorted = np.flip(np.sort(umi_count))
            target_cells = min(len(umi_count_sorted), num_cells_above_max+target_cells)
            umi_min = umi_count_sorted[num_cells_above_max+target_cells]
        filter_ind = np.logical_and(umi_count > umi_min, umi_count < umi_max)
        self.filtered_obj = np.logical_and(self.filtered_obj, filter_ind)
        # plot barcode rank
        if plot:
            umi_count = np.flip(np.sort(umi_count))
            ind = np.logical_and(umi_count > umi_min, umi_count < umi_max)
            plt.scatter(np.array(range(len(umi_count)))[np.invert(ind)], umi_count[np.invert(ind)], color='k')
            plt.scatter(np.array(range(len(umi_count)))[ind], umi_count[ind])
            plt.xscale('log')
            plt.yscale('log')
            plt.title('barcode rank plot')
            plt.xlabel('barcode')
            plt.ylabel('UMI count')
            plt.grid(visible=True, which='major')
            plt.show()

    def filter_by_gene_dispersion(self, min_dispersion=1, target_number=None, plot=True):
        # Mark genes (columns) below min_dispersion for exclusion, using only the
        # currently-filtered cells (self.filtered_obj) to compute mean/variance.
        # If target_number is given, min_dispersion is instead set so exactly the
        # target_number genes with the highest dispersion survive.
        data_mat = self.m[self.filtered_obj, :]
        # get mean along matrix columns
        mean = np.mean(data_mat, axis=0)
        # get variance along matrix columns
        var = np.var(data_mat, axis=0)
        # get dispersion, if mean is 0, set dispersion to 0
        dispersion = var / mean
        dispersion[np.isnan(dispersion)] = 0
        # sort dispersion
        sorted_dispersion = np.flip(np.sort(dispersion))
        if target_number is not None:
            min_dispersion = sorted_dispersion[target_number]

        filter_ind = dispersion > min_dispersion
        self.filtered_var = np.logical_and(filter_ind, self.filtered_var)
        # plot gene cv, color genes that are filtered out in red
        # note: excluded (red) points use sqrt(dispersion) but kept (blue) points use
        # sqrt(sqrt(dispersion)) - a 4th root, not a square root - so the two point
        # sets are not on the same y-scale; see FINDINGS in logs/src/data_functions.txt.
        if plot:
            plt.scatter(mean[np.invert(filter_ind)], np.sqrt(dispersion[np.invert(filter_ind)]), color='r')
            plt.scatter(mean[filter_ind], np.sqrt(np.sqrt(dispersion[filter_ind])), color='b')
            plt.xscale('log')
            plt.yscale('log')
            plt.title('gene dispersion plot')
            plt.xlabel('mean')
            plt.ylabel('dispersion')
            plt.grid(visible=True, which='major')
            plt.show()


    def filter_by_mean_expression(self, min_mean=1e-2, target_number=None, plot=True):
        # Mark genes (columns) with mean expression below min_mean for exclusion
        # (over currently-filtered cells). If target_number is given, min_mean is
        # instead set so exactly the target_number highest-mean genes survive.
        data_mat = self.m[self.filtered_obj, :]
        # get mean along matrix columns
        mean = np.mean(data_mat, axis=0)
        # get variance along matrix columns
        std = np.std(data_mat, axis=0)
        # sort mean
        sorted_mean = np.flip(np.sort(mean))
        if target_number is not None:
            min_mean = sorted_mean[target_number]

        filter_ind = mean > min_mean
        self.filtered_var = np.logical_and(filter_ind, self.filtered_var)
        # plot gene cv, color genes that are filtered out in red
        if plot:
            plt.scatter(mean[np.invert(filter_ind)], std[np.invert(filter_ind)], color='r')
            plt.scatter(mean[filter_ind], std[filter_ind], color='b')
            plt.xscale('log')
            plt.yscale('log')
            plt.title('gene dispersion plot')
            plt.xlabel('mean')
            plt.ylabel('dispersion')
            plt.grid(visible=True, which='major')
            plt.show()

    def get_gene_stats(self):
        # mean/median number of distinct genes detected (non-zero) per cell, over
        # currently-filtered cells - a standard scRNA-seq QC summary statistic
        data_mat = self.m[self.filtered_obj, :]
        # binarize the matrix
        bin_mat = np.zeros(data_mat.shape)
        bin_mat[data_mat > 0] = 1
        # get the mean and median number of genes per cell
        mean_genes = np.mean(np.sum(bin_mat, axis=1))
        median_genes = np.median(np.sum(bin_mat, axis=1))
        return mean_genes, median_genes

    def remove_highly_expressed_genes(self, percentile=90):
        # Mark genes whose mean expression is at or above the given percentile for
        # exclusion (e.g. percentile=90 drops the top 10% most highly expressed
        # genes) - used to down-weight the influence of a few dominant genes.
        data_mat = self.m[self.filtered_obj, :]
        mean_expression = np.mean(data_mat, axis=0)
        max_val = np.percentile(mean_expression, percentile)
        filter_ind = mean_expression < max_val
        self.filtered_var = np.logical_and(filter_ind, self.filtered_var)


    def get_filtered_matrix(self):
        # Materialize a new AnnMat containing only the rows/columns currently
        # marked True in filtered_obj/filtered_var (the lazy filters are applied
        # here, once). Note: batch is not carried over to the returned AnnMat
        # (it defaults to None) even though self.batch may hold per-row labels -
        # see FINDINGS in logs/src/data_functions.txt.
        m = self.m[self.filtered_obj, :][:, self.filtered_var]
        obj_names = self.obj_names[self.filtered_obj]
        var_names = self.var_names[self.filtered_var]
        return AnnMat(m, obj_names, var_names)

    def reset_filters(self):
        # clear both filters back to "everything included"
        self.filtered_obj = np.full(len(self.obj_names), True)
        self.filtered_var = np.full(len(self.var_names), True)


    def sort_names_by_expression(self):
        # gene names among currently-filtered cells, ordered from highest to lowest
        # mean expression
        data_mat = self.m[self.filtered_obj, :]
        mean_expression = np.mean(data_mat, axis=0)
        # return the sorted gene names
        sorted_indices = np.argsort(mean_expression)
        sorted_genes = self.var_names[np.flip(sorted_indices)]
        # return the sorted genes
        return sorted_genes


    def random_knockout(self):
        # For each currently-filtered cell with more than 10 detected genes,
        # remove exactly one read, chosen at random with probability proportional
        # to each gene's share of that cell's counts (i.e. higher-expressed genes
        # are more likely to lose a read) - used to test robustness of downstream
        # statistics to a small, expression-weighted perturbation of the counts.
        # Cells with <=10 detected genes are left untouched.
        data_mat = self.m[self.filtered_obj, :]
        for i in range(data_mat.shape[0]):
            # get the indices of non-zero elements
            non_zero_indices = np.nonzero(data_mat[i, :])[0]
            if len(non_zero_indices) > 10:
                # randomly select one index to knock out
                probabilities = data_mat[i, non_zero_indices] / np.sum(data_mat[i, non_zero_indices])
                index_to_knockout = np.random.choice(non_zero_indices, p=probabilities)
                data_mat[i, index_to_knockout] -= 1
        # return the annotated matrix with the knocked out reads
        return AnnMat(data_mat, self.obj_names[self.filtered_obj], self.var_names[self.filtered_var], self.batch)

def read_from_csv(path):
    # Load a previously-saved AnnMat (e.g. written by AnnMat.to_csv) back from a
    # delimited text file: first column = object (cell) names, first row =
    # variable (gene) names. Any NaN cell (e.g. from a ragged/missing entry) is
    # treated as zero counts. Delimiter is guessed from the path: '.csv' -> comma,
    # anything else -> tab.
    if 'csv' in path:
        df = pd.read_csv(path, index_col=0)
    else:
        df = pd.read_csv(path, sep='\t', index_col=0)
    obj_names = df.index.to_numpy()
    var_names = df.columns.to_numpy()
    m = df.to_numpy()
    m[np.isnan(m)] = 0
    return AnnMat(m, obj_names, var_names)


def transpose_matrix(amat):
    # swap rows/columns of an AnnMat (e.g. to treat genes as objects and cells as
    # variables); the transposed matrix's filters are reset to "all included"
    return AnnMat(amat.m.T, amat.var_names, amat.obj_names)


def gene_id(val):
    # Collapse a probe/feature ID down to its shared gene identifier. Probe IDs
    # from this platform are of the form '<gene>_<probe-number>' when the probed
    # gene has multiple probes (trailing '_xx' token is purely numeric) - that
    # numeric suffix is stripped so all probes for one gene are counted together.
    # If the last '_'-separated token is not numeric, val is assumed to already be
    # a plain gene id and is returned unchanged.
    x = re.split('_', val)
    if bool(re.findall(r'\d+', x[-1:][0])):
        return '_'.join(x[:-1])
    else:
        return val


def get_annotated_data(path, method='max_probe'):
    # Build a dense cell x gene AnnMat from a raw per-probe count CSV (one row per
    # cell-barcode/probe/UMI-count observation). method controls how multiple
    # probes mapping to the same gene (see gene_id) are combined:
    #   'max_probe'   - take the probe with the highest UMI count for that gene
    #   'sum_probes'  - sum UMI counts across all probes for that gene
    df = pd.read_csv(path)
    df.columns = ['cell_barcode', 'Feature_ID', 'Feature_name', 'Feature_type', 'UMI_count']
    # create list of cell barcodes
    cell_bar = df.cell_barcode.to_numpy()
    cell_bar = np.unique(cell_bar)
    cell_bar = np.ndarray.tolist(cell_bar)
    # combine all gene identifiers that end with '_xx' and create a unique list of features
    features = []
    for val in df.Feature_ID:
        temp = gene_id(val)
        if temp not in features:
            features.append(temp)
    features.sort()
    # create a cell-gene matrix
    M = np.zeros((len(cell_bar), len(features)), dtype=int)
    for i in range(len(df.cell_barcode)):
        cell_index = cell_bar.index(df.cell_barcode[i])
        gene_index = features.index(gene_id(df.Feature_ID[i]))
        if method == 'max_probe':
            # keep the highest-count probe seen so far for this (cell, gene) pair
            if M[cell_index, gene_index] < df.at[i, 'UMI_count']:
                M[cell_index, gene_index] = df.at[i, 'UMI_count']
        elif method == 'sum_probes':
            M[cell_index, gene_index] = M[cell_index, gene_index] + df.at[i, 'UMI_count']
        # print progress
        if i % 1000 == 0:
            print(str(round(100 * i / len(df.cell_barcode))) + '% of run', end='\r', flush=True)
    # create AnnMat object from matrix and annotations
    return AnnMat(M, cell_bar, features)


def plot_eig_dist(pcs, pcs1, N, x_max, y_max, n_bins, ax=None, x_label=True, y_label=True):
    # Plot the empirical eigenvalue distribution (pcs, from the real data) against
    # the scrambled/null distribution (pcs1) and the analytic Marchenko-Pastur
    # density (analysis_functions.mp_distribution), as stacked-area histograms.
    # pcs / pcs1 are the eigenvalue arrays returned by analysis_functions.get_eig_dist;
    # N is the number of cells (rows) of the matrix they were computed from.
    P = len(pcs)
    scale = 1 # scale factor for the Marchenko-Pastur distribution
    edges = np.linspace(-0.1, x_max, num=n_bins)

    # remove zeros in pcs and pcs1
    # aspect ratio P/N > 1 (more genes than cells): MP density is defined for the
    # p/n <= 1 convention, so genes/cells are swapped (scale = N/P) and the
    # structural zero eigenvalues (rank-deficient directions) are dropped before
    # histogramming, matching the analytic MP formula's normalization.
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
    y = [af.mp_distribution(val, P / N) for val in x]
    ax.plot(x, y, color='#756bb1', linestyle='dashed', label='MP')
    # labels and limits
    if x_label:
        ax.set_xlabel("$\lambda$")
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


def pcs_to_csv(amat, output_dir, file_name, norm='sum', log=False, norm_method='sum', norm_sum=1):
    # Compute the empirical eigenvalue spectrum of amat and write its histogram
    # (as an X/Y density table) to <output_dir>/<file_name>_N<N>_P<P>.csv.
    # NOTE: get_eig_dist returns a 3-tuple (pcs, pcs1, fraction_non_zero) - see
    # FINDINGS in logs/src/data_functions.txt, this unpacking is broken as written.
    m = amat.m
    N, P = m.shape
    pcs, pcs1 = af.get_eig_dist(m, norm=norm, log=log, norm_method=norm_method, norm_sum=norm_sum)
    pcs = pcs[pcs > 0]
    # get histogram values of pcs
    scale = 1  # scale factor for the Marchenko-Pastur distribution
    x_max = 50
    n_bins = 501
    edges = np.linspace(0, x_max, num=n_bins)
    bin_diff = edges[1] - edges[0]
    # remove zeros in pcs and pcs1
    # if alpha>1 adjust the scale factor to match theoretical results
    if P / N > 1:
        scale = N / P
        pcs = pcs[pcs != 0]

    # calculate histogram
    counts, bins = np.histogram(pcs, bins=edges, density=True)
    df = pd.DataFrame({'X': bins[:-1] + bin_diff / 2, 'Y': scale * counts})
    # save to csv
    file_name = file_name + f'_N{N}_P{P}.csv'
    df.to_csv(os.path.join(output_dir, file_name), index=False)


def equate_dims(amat1, amat2, target_cells):
    # Make two AnnMat datasets directly comparable by restricting both to their
    # shared genes (via set intersection, not a precondition of equality) and
    # subsampling both to the same number of cells, target_cells (capped by
    # whichever input has fewer cells/genes to draw from). Filters are set on the
    # input objects in place; the *filtered* matrices are returned, unmaterialized
    # copies of amat1/amat2 remain in their original, now-mutated filter state.
    # np.random.seed(0) is called once, globally, before both cell draws - this
    # makes the subsample reproducible run-to-run, but also resets the global
    # numpy RNG state as a side effect for any unseeded random call made
    # afterwards elsewhere in the same process.
    target_cells = min(len(amat1.obj_names), len(amat2.obj_names), target_cells)
    # get the intersection of the gene names
    intersect, indices1, indices2 = np.intersect1d(amat1.var_names, amat2.var_names, return_indices=True)
    amat1.filtered_var = [val in indices1 for val in range(len(amat1.var_names))]
    amat2.filtered_var = [val in indices2 for val in range(len(amat2.var_names))]
    # random subsample of cells to target_cells
    np.random.seed(0)
    indices1 = np.random.choice(range(len(amat1.obj_names)), target_cells, replace=False)
    indices2 = np.random.choice(range(len(amat2.obj_names)), target_cells, replace=False)
    # filter the matrices
    amat1.filtered_obj = [val in indices1 for val in range(len(amat1.obj_names))]
    amat2.filtered_obj = [val in indices2 for val in range(len(amat2.obj_names))]
    # return filtered matrices
    return amat1.get_filtered_matrix(), amat2.get_filtered_matrix()


def remove_non_protein_coding_genes(amat):
    # Drop genes that are not protein-coding (or pseudogenes) according to a
    # biotype/synonym lookup table, matched case-insensitively and after
    # stripping a known probe-name prefix ('lelobekk_'). Assumes the process's
    # current working directory's *parent* contains a filtered_data/ folder with
    # biotype_map_syn.csv (see bulk_functions.py's similar metadata/ assumption in
    # CLAUDE.md - run notebooks from scripts/ so this path resolves correctly).
    dir = os.path.join(os.path.dirname(os.getcwd()), 'filtered_data')
    file_name = 'biotype_map_syn.csv'
    path = os.path.join(dir, file_name)
    # read from csv
    biotype_map = pd.read_csv(path, index_col=0, dtype=str)
    biotype_map.fillna('', inplace=True)
    syn_array = get_synonym_array(biotype_map)
    # get casefold list of gene names:
    var_names = [val.casefold() for val in amat.var_names]
    # remove pattern 'lelobekk_' from var_names
    var_names = [val.replace('lelobekk_', '') for val in var_names]
    # get gene type
    gene_types = [get_gene_type(val, biotype_map, syn_array) for val in var_names]
    filter_protein_coding = [np.logical_or(g_type == 'protein_coding', g_type == 'pseudogene') for g_type in gene_types]
    not_in_list = [g_type == 'not found' for g_type in gene_types]
    # filter_protein_coding = np.logical_or(filter_protein_coding, not_in_list)
    amat.filtered_var = np.logical_and(filter_protein_coding, amat.filtered_var)
    print(amat.var_names[np.invert(filter_protein_coding)])
    return amat.get_filtered_matrix()


def get_gene_type(gene, biotype_map, synonym_array):
    # Look up gene's biotype (e.g. 'protein_coding') by finding which row of
    # synonym_array contains it (any synonym or the canonical name matches), and
    # reading biotype_map at that row index. Returns 'not found' if no row
    # contains the gene, or if any other lookup error occurs (broad except).
    gene = gene.casefold()
    try:
        index = np.where(synonym_array == gene)[0][0]
        g_type = biotype_map.loc[index, 'biotype']
    except:
        g_type = 'not found'

    return g_type


def get_synonym_array(biotype_map):
    # Turn biotype_map (one row per gene, with a comma-separated 'synonym' column
    # and a 'gene' column) into a dense 2-D object array, one row per gene, where
    # each row lists that gene's synonyms plus its canonical name (all casefolded)
    # padded with zeros to the longest row - so get_gene_type can search it with a
    # single vectorized np.where(synonym_array == gene) instead of per-row parsing.
    # convert to list of lists:
    list_of_lists = []
    for i in range(len(biotype_map)):
        l = str(biotype_map.loc[i, 'synonym']).split(',')
        l = [val.casefold() for val in l]
        l.append(str(biotype_map.loc[i, 'gene']).casefold())
        list_of_lists.append(l)
    max_len = max([len(val) for val in list_of_lists])
    np_array = np.zeros((len(list_of_lists), max_len), dtype=object)
    for i in range(len(list_of_lists)):
        for j in range(len(list_of_lists[i])):
            np_array[i, j] = list_of_lists[i][j]
    return np_array


def project_on_pcs(amat, n_pcs=50):
    # Project amat.m onto its own top n_pcs principal components: normalize rows
    # to a fixed total (target_sum=1e4, a standard scRNA-seq CPM-style depth
    # normalization), log-transform, z-transform columns, then take the SVD and
    # keep the leading n_pcs right-singular directions. Returns an
    # (n_cells x n_pcs) score matrix used as the input to UMAP in get_umap.
    data = af.normalize(amat.m, method='sum', target_sum=1e4)
    data = af.log_transform(data)
    data = af.z_transform(data)
    u, s, vh = np.linalg.svd(data, full_matrices=False)
    return np.dot(data, vh[:n_pcs, :].T)


def get_umap(amat, n_neighbors=15, min_dist=0.1, n_components=2, n_pcs=50):
    # 2-D (by default) UMAP embedding of amat's cells, computed on their top n_pcs
    # PCA coordinates (see project_on_pcs) rather than on the raw gene matrix -
    # the usual scRNA-seq practice of denoising via PCA before UMAP.
    data = project_on_pcs(amat, n_pcs=n_pcs)
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components)
    embedding = reducer.fit_transform(data)
    return embedding


def concat_datasets(datasets):
    # Stack multiple AnnMat datasets into one, on the union of their gene sets
    # (a gene absent from a given dataset is filled with zeros for that
    # dataset's cells). Returns a single AnnMat whose `batch` array records
    # which input dataset (0-indexed, in `datasets` order) each cell came from.
    #
    # IMPORTANT: gene_names = np.union1d(...) is always alphabetically sorted, and
    # each dataset's rows are written into gene_names-sorted column positions via
    # np.isin/np.where - this only lines up correctly with dataset.m's own column
    # order if dataset.var_names is *also* already sorted (true for AnnMat objects
    # built by get_annotated_data, which sorts `features` before building M, and
    # preserved by get_filtered_matrix's boolean masking). Passing a dataset whose
    # var_names are not sorted will silently misassign gene values into the wrong
    # columns - see FINDINGS in logs/src/data_functions.txt.
    gene_names = np.array([])
    for dataset in datasets:
        gene_names = np.union1d(gene_names, dataset.var_names)
    # create an empty matrix
    n_genes = len(gene_names)
    n_cells = 0
    for dataset in datasets:
        n_cells += len(dataset.obj_names)
    M = np.zeros((n_cells, n_genes))
    # fill the matrix
    # store info on original dataset in batch
    batch = np.zeros(n_cells)
    batch_index = 0
    cell_index = 0
    cell_names = np.array([])
    for dataset in datasets:
        for i in range(len(dataset.obj_names)):
            # cell_index can never reach n_cells here since M/n_cells is sized to
            # exactly the sum of all datasets' cell counts - this check is inert
            if cell_index == n_cells:
                break
            cell_names = np.append(cell_names, dataset.obj_names[i])
            gene_index = np.where(np.isin(gene_names, dataset.var_names))[0]
            M[cell_index, gene_index] = dataset.m[i, :]
            # update batch info
            batch[cell_index] = batch_index
            cell_index += 1
        batch_index += 1
        print(f"batch {batch_index} out of {len(datasets)} completed")
    return AnnMat(M, cell_names, gene_names, batch)


def get_histogram_data(pcs,N,P):
    # Turn an eigenvalue array pcs (from analysis_functions.get_eig_dist, for an
    # N-cell x P-gene matrix) into a density-histogram DataFrame (X = bin center,
    # Y = density), applying the same P/N>1 rescaling as plot_eig_dist / pcs_to_csv
    # so it is directly comparable to the analytic Marchenko-Pastur curve.
    pcs = pcs[pcs > 0]
    # get histogram values of pcs
    scale = 1  # scale factor for the Marchenko-Pastur distribution
    x_max = 50
    n_bins = 501
    edges = np.linspace(0, x_max, num=n_bins)
    bin_diff = edges[1] - edges[0]
    # remove zeros in pcs and pcs1
    # if alpha>1 adjust the scale factor to match theoretical results
    if P / N > 1:
        scale = N / P
    # calculate histogram values
    counts, bins = np.histogram(pcs, bins=edges, density=True)
    pcs_df = pd.DataFrame({'X': bins[:-1] + bin_diff / 2, 'Y': scale * counts})
    return pcs_df