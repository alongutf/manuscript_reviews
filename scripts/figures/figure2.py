"""Figure 2: the eigenvalue-spectrum method behind GMP-Cor, from theory to data.

Assembles a seven-panel figure that walks the reader through the correlation-
spectrum argument used throughout the paper:
  A - Marchenko-Pastur (MP) law for an uncorrelated random matrix, as a sanity
      check that the theory matches a simulated null.
  B - the Generalized MP (GMP) family: how the spectrum's shape changes as the
      underlying gene-gene correlation strength (chi) increases.
  C - a small synthetic sparse count matrix and its per-gene permutation-
      scrambled counterpart, to make the empirical-vs-null comparison concrete.
  D - eigenvalue density (PDF) of that synthetic data, colour-coded into MP-null,
      other-spurious, and true-signal eigenvalues, with a CCDF inset.
  E-G - CCDF of the real correlation spectrum (signal vs. scrambled null) for the
      paper's own exponential-phase E. coli sample and two published reference
      datasets (Ma et al. E. coli and K. pneumoniae), each annotated with its
      permutation-test p-value.

Run interactively (e.g. from scripts/figures/, as a Jupyter cell or via
`run figure2.py`) with the working directory two levels below the repo root, so
that `root_dir` below resolves to the repo root.

Input:
  ev_data/simulated_pcs.npy                 - synthetic-data eigenvalues (panel D)
  ev_data/Expira_biorep_t0A_filtered.npy     - exponential E. coli spectrum (panel E)
  ev_data/deb_Ec_CDS_untreated.npy           - Ma et al. E. coli spectrum (panel F)
  ev_data/deb_KP_CDS_untreated.npy           - Ma et al. K. pneumoniae spectrum (panel G)
  each .npy is a (2, n_genes) array: row 0 = empirical eigenvalues, row 1 =
  scrambled-null eigenvalues (see src.analysis_functions.get_eig_dist)
  "model fit"/model_alpha2_sigma*.txt        - precomputed GMP PDF curves (panel B)

Output:
  figure2.svg   (300 dpi, publication)
  figure2_preview.png (200 dpi, quick look)
"""

import src.analysis_functions as af
import src.data_functions as df
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch
import seaborn as sns
from figure_functions import PanelFigure
import permutation_pvalues as pv
import numpy as np
import os
import importlib

importlib.reload(af)
importlib.reload(df)
# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------
fsize = 10
plt.close("all")
# repo root is two levels above this script's usual working directory
# (scripts/figures/); see the module docstring for the assumed cwd
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
ev_data_dir = os.path.join(root_dir, 'ev_data')


def panel_A(ax):
    # Random N×P Gaussian matrix correlation spectrum vs. MP distribution
    N = 500
    P = 1000
    matrix = np.random.randn(N, P)
    # sample (Pearson-style) correlation matrix of P independent Gaussian genes
    # over N samples; no true correlation exists, so its spectrum is pure MP noise
    corr_matrix = matrix.T @ matrix / N
    eigvals, _ = np.linalg.eig(corr_matrix)
    # drop numerically-zero eigenvalues (P - N of them are exactly rank-deficient
    # since corr_matrix has rank <= N < P) before histogramming
    eigvals = np.real(eigvals[eigvals > 1e-6])
    bins = np.linspace(0, 6, 40)
    # weight each count by 1/(P * bin width) so the histogram is a density over
    # eigenvalue index normalised by P, directly comparable to mp_distribution
    ax.hist(eigvals,
            weights=np.ones_like(eigvals) / (P * (bins[1] - bins[0])),
            bins=bins, color='r', alpha=0.5, density=False,
            label='simulated')
    x = np.linspace(0, 6, 100)
    # analytic Marchenko-Pastur PDF at aspect ratio a = P/N, overlaid on the
    # simulated histogram as the theoretical prediction for uncorrelated data
    y = np.array([af.mp_distribution(val, P / N) for val in x])
    ax.plot(x, y, color='r', linewidth=1, linestyle='--', label='MP')
    ax.set_title('Random matrix\nMP distribution', fontsize=fsize-2, pad=0)
    ax.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize - 2, labelpad=0)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax.grid(False)
    ax.legend(fontsize=fsize - 2)


def panel_B(ax):
    # GMP distribution: random Gaussian matrix with underlying correlations
    # each file is a precomputed GMP PDF (columns: lambda, rho(lambda)) at aspect
    # ratio alpha=2 and correlation strength chi=sigma, from model_fit.nb (Mathematica)
    files = ['model_alpha2_sigma0.txt', 'model_alpha2_sigma07.txt',
             'model_alpha2_sigma08.txt', 'model_alpha2_sigma09.txt']
    labels = [r'$\chi=0$ (MP)', r'$\chi=0.7$', r'$\chi=0.8$', r'$\chi=0.9$']
    cmap = plt.cm.RdBu
    colors = [cmap(i) for i in [0.1, 0.7, 0.85, 0.95]]
    for i, file in enumerate(files):
        data = np.loadtxt(os.path.join(root_dir, 'model fit', file))
        ax.plot(data[:, 0], data[:, 1], label=labels[i], color=colors[i], linewidth=1.5)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize - 2, labelpad=0)
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 0.35)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax.set_title('Random matrix\nwith correlations\nGeneralized MP', fontsize=fsize-2, pad=0)
    ax.legend(fontsize=fsize - 2)


def panel_C(axes):
    # Illustrative sparse simulation matrix + scrambled version
    np.random.seed(42)   # fixed seed so the illustrative matrix is reproducible
    n_cells, n_genes = 22, 14

    # small, purely illustrative dataset: each gene is nonzero in a random subset
    # of cells (15-38% detection rate) with exponentially distributed counts,
    # mimicking the sparsity of real scRNA-seq probe counts
    matrix = np.zeros((n_cells, n_genes))
    for j in range(n_genes):
        n_nz = max(1, int(n_cells * np.random.uniform(0.15, 0.38)))
        rows = np.random.choice(n_cells, size=n_nz, replace=False)
        matrix[rows, j] = np.random.exponential(2.5, size=n_nz)

    # the permutation null: shuffle each gene's values independently across cells,
    # which destroys gene-gene correlations while keeping each gene's own
    # marginal distribution (and overall sparsity) unchanged
    scrambled = matrix.copy()
    for j in range(n_genes):
        perm = np.random.permutation(n_cells)
        scrambled[:, j] = matrix[perm, j]

    vmax = matrix.max()

    axes[0, 0].imshow(matrix, aspect='auto', cmap='Greys', vmin=0, vmax=vmax,
                      interpolation='nearest')
    axes[0, 0].set_xticks([])
    axes[0, 0].set_yticks([])
    axes[0, 0].set_xlabel('genes', fontsize=fsize - 2, labelpad=2)
    axes[0, 0].set_ylabel('cells', fontsize=fsize - 2, labelpad=2)
    axes[0, 0].set_title('Synthetic data\nwith correlations', fontsize=fsize-2, pad=0)

    axes[1, 0].imshow(scrambled, aspect='auto', cmap='Greys', vmin=0, vmax=vmax,
                      interpolation='nearest')
    axes[1, 0].set_xticks([])
    axes[1, 0].set_yticks([])
    axes[1, 0].set_xlabel('genes', fontsize=fsize - 2, labelpad=2)
    axes[1, 0].set_ylabel('cells', fontsize=fsize - 2, labelpad=2)
    #axes[1, 0].set_title('Scrambled data', fontsize=fsize)


def panel_D(axes):
    # Simulation eigenvalue distributions: original (color-coded) and scrambled
    # row 0 = empirical (correlated) spectrum, row 1 = per-gene-permuted null
    # spectrum, both from the same synthetic dataset used to validate GMP-Cor
    pcs_data = np.load(os.path.join(ev_data_dir, 'simulated_pcs.npy'))
    pcs, pcs1 = pcs_data[0], pcs_data[1]

    data1 = pcs[pcs > 0]
    data2 = pcs1[pcs1 > 0]
    bin_width = 0.15
    x1 = (1 + np.sqrt(2)) ** 2 + bin_width   # MP upper edge (γ = P/N = 2)
    # empirical scrambled-null maximum eigenvalue -- this IS the GMP-Cor threshold
    # (lambda_max^scr): empirical eigenvalues above it are counted as true signal
    x2 = float(np.max(pcs1))+bin_width     # scrambled maximum = GMP-Cor threshold

    all_data = np.concatenate([data1, data2])
    bin_edges = np.arange(min(all_data), max(all_data) + bin_width, bin_width)

    ax_top = axes[0, 0]
    ax_bot = axes[1, 0]

    # Top: original pcs with three-color coding
    _, _, patches1 = ax_top.hist(
        data1, bins=bin_edges, width=bin_width * 0.8, align='right',
        edgecolor='black', color='#d9d9d9', alpha=0.7, density=True)
    # recolor each bar by which of the three eigenvalue regimes its left edge
    # falls in: below the analytic MP edge (pure noise, expected even with no
    # scrambling), between the MP edge and the scrambled maximum (extra spurious
    # correlation from sparsity/finite-size effects, caught by the null), and
    # above the scrambled maximum (true correlation signal, counted in GMP-Cor)
    for patch in patches1:
        bx = patch.get_x()
        if bx < x1:
            patch.set_facecolor('darkgray')
        elif bx < x2:
            patch.set_facecolor('salmon')
        else:
            patch.set_facecolor('skyblue')

    # Threshold lines on BOTH subplots
    ax_top.axvline(x1, color='k', linestyle='--', alpha=0.6, label=r'$\lambda_{max}^{MP}$')
    ax_top.axvline(x2, color='dimgray', linestyle=':', alpha=0.8, label=r'$\lambda_{max}^{scr}$')
    ax_top.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize - 2)
    ax_top.set_title('Correlation eigenvalue density', fontsize=fsize-2, pad=0)
    ax_top.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax_top.set_xlim([2, 12])
    ax_top.set_ylim(0, 0.2)
    ax_top.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax_top.grid(False)

    # CCDF inset on the top subplot
    inset = ax_top.inset_axes([0.37, 0.4, 0.54, 0.57])
    _draw_ccdf(inset, data1, data2, show_legend=True, markersize=2)
    inset.set_xlabel(r'$\lambda$', fontsize=fsize - 3, labelpad=0)
    inset.set_ylabel('CCDF', fontsize=fsize - 3, labelpad=0)
    inset.tick_params(labelsize=fsize - 3)
    inset.set_xlim([0.2,14])

    # Bottom: scrambled pcs1 (same x-axis), color intermediate (sparsity) bars too.
    # scrambled data has no eigenvalues above x2 by construction (x2 is its own
    # max), so only the below/above-MP-edge split applies here.
    _, _, patches2 = ax_bot.hist(
        data2, bins=bin_edges, width=bin_width * 0.8, align='right',
        edgecolor='black', color='darkgray', alpha=0.7, density=True)
    for patch in patches2:
        bx = patch.get_x()
        if bx < x1:
            patch.set_facecolor('darkgray')
        else:
            patch.set_facecolor('salmon')
    ax_bot.text(x1,0.1,r'$\lambda_{max}^{MP}$', fontsize=fsize - 3, horizontalalignment='right')
    ax_bot.text(x2,0.1,r'$\lambda_{max}^{scr}$', fontsize=fsize - 3, horizontalalignment='right')
    ax_bot.axvline(x1, color='k', linestyle='--', alpha=0.6)
    ax_bot.axvline(x2, color='dimgray', linestyle=':', alpha=0.8)
    ax_bot.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax_bot.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize - 2)
    ax_bot.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax_bot.set_xlim([2, 12])
    ax_bot.set_ylim(0, 0.2)
    ax_bot.tick_params(axis='both', which='major', labelsize=fsize - 2)
    ax_bot.grid(False)

    # Shared legend for both subplots
    legend_handles = [
        Patch(facecolor='darkgray', edgecolor='black', label='MP spurious correlations'),
        Patch(facecolor='salmon', edgecolor='black', label='spurious correlations'),
        Patch(facecolor='skyblue', edgecolor='black', label='true correlation signal'),
    ]
    ax_bot.legend(handles=legend_handles, fontsize=fsize - 3,
                  loc='upper right', framealpha=1)


def _draw_ccdf(ax, data1, data2, dataset=None, show_legend=True, markersize=3):
    """Draw the loglog CCDF of original (data1) vs scrambled (data2) eigenvalues.

    Shared color scheme across all CCDF plots:
      grey  -> spurious correlations (eigenvalues below scrambled threshold)
      blue  -> true correlations     (eigenvalues at/above scrambled threshold)
      black -> scrambled data
    """
    x2 = float(np.max(data2))
    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    # empirical CCDF P(X >= x_i) at each sorted value; the "+ 1/p1" shifts the
    # largest point up to 1/p1 instead of 0, which would be unplottable on a
    # log axis
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    spurious = d1s < x2
    ax.loglog(d1s[spurious], ccdf1[spurious], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='spurious',
              markersize=markersize)
    ax.loglog(d1s[~spurious], ccdf1[~spurious], '.', linestyle='-',
              color='skyblue', label='signal', markersize=markersize)
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled',
              markersize=markersize)
    ax.set_xlim([0.1, np.max(d1s) * 1.5])
    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
    if show_legend:
        # p-value rides in the legend box, under the three series keys
        pv.legend_with_p(ax, dataset, fontsize=fsize - 2)


def panel_E(ax):
    # Regulated dataset: CCDF only (no PDF, no inset)
    _plot_ccdf(ax, 'Expira_biorep_t0A_filtered.npy', r'Exponential $\it{E. coli}$')


def _plot_ccdf(ax, npy_file, title):
    """Load a (2, n_genes) empirical/scrambled eigenvalue array and draw its CCDF."""
    arr = np.load(os.path.join(ev_data_dir, npy_file))
    data1 = arr[0, :];  data1 = data1[data1 > 0]
    data2 = arr[1, :];  data2 = data2[data2 > 0]
    _draw_ccdf(ax, data1, data2, dataset=npy_file)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)
    ax.set_title(title, fontsize=fsize - 2, pad=0)
    ax.tick_params(labelsize=fsize - 2)


def panel_F(ax):
    _plot_ccdf(ax, 'deb_Ec_CDS_untreated.npy', 'Exponential $\it{E. coli}$, Ma et. al')


def panel_G(ax):
    _plot_ccdf(ax, 'deb_KP_CDS_untreated.npy', 'Exponential $\it{K. pneumoniae}$,\nMa et. al')


# ------------------------------------------------------------------
# ASSEMBLE FIGURE
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 6.5), label_offset=(-0.03, 0.04))

panel_pos = [
    [0.08, 0.78, 0.19, 0.16],   # A – MP histogram only (single)
    [0.08, 0.44, 0.19, 0.21],   # B – GMP curves (single)
    [0.33, 0.42, 0.12, 0.52],   # C – sparse simulation (2×1 grid; matches A+B height)
    [0.55, 0.42, 0.44, 0.52],   # D – simulation eigenvalues (2×1 grid; matches A+B height)
    [0.08, 0.08, 0.24, 0.24],   # E – CCDF (our Exponential E. coli)   ← single row
    [0.40, 0.08, 0.24, 0.24],   # F – CCDF (Ma et al. E. coli)         ← single row
    [0.72, 0.08, 0.24, 0.24],   # G – CCDF (Ma et al. K. pneumoniae)   ← single row
]

# Panel A
pf.add_panel(panel_pos[0], draw_func=panel_A)

# Panel B
pf.add_panel(panel_pos[1], draw_func=panel_B)

# Panel C
axes_panel_C = pf.add_grid_panel(panel_pos[2], 2, 1, hspace=0.30)
panel_C(axes_panel_C)

# Panel D
axes_panel_D = pf.add_grid_panel(panel_pos[3], 2, 1, hspace=0.3)
panel_D(axes_panel_D)

# Panel E
pf.add_panel(panel_pos[4], draw_func=panel_E)

# Panel F
pf.add_panel(panel_pos[5], draw_func=panel_F)

# Panel G
pf.add_panel(panel_pos[6], draw_func=panel_G)

pf.save("figure2.svg", dpi=300)
pf.save("figure2_preview.png", dpi=200)
plt.show()
