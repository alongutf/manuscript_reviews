import src.analysis_functions as af
import src.data_functions as df
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch
import seaborn as sns
from figure_functions import PanelFigure
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
root_dir = os.path.dirname(os.path.dirname(os.getcwd()))
ev_data_dir = os.path.join(root_dir, 'ev_data')


def panel_A(axes):
    # Random N×P Gaussian matrix and its correlation spectrum vs. MP distribution
    N = 500
    P = 1000
    matrix = np.random.randn(N, P)
    corr_matrix = matrix.T @ matrix / N
    sns.heatmap(corr_matrix[:20, :20], ax=axes[0, 0], cmap='bwr', center=0,
                cbar=True, vmin=-.5, vmax=.5)
    colorbar = axes[0, 0].collections[0].colorbar
    colorbar.set_ticks([-0.5, 0, 0.5])
    colorbar.set_label('correlation', fontsize=fsize - 2, labelpad=0)
    colorbar.ax.tick_params(labelsize=fsize - 2)

    eigvals, _ = np.linalg.eig(corr_matrix)
    eigvals = np.real(eigvals[eigvals > 1e-6])
    bins = np.linspace(0, 6, 40)
    axes[1, 0].hist(eigvals,
                    weights=np.ones_like(eigvals) / (P * (bins[1] - bins[0])),
                    bins=bins, color='red', alpha=0.5, density=False,
                    label='simulated\neigenvalues')
    x = np.linspace(0, 6, 100)
    y = np.array([af.mp_distribution(val, P / N) for val in x])
    axes[1, 0].plot(x, y, color='red', linewidth=1, label='MP')
    axes[0, 0].set_title('Random matrix', fontsize=fsize)
    axes[0, 0].set_yticks([])
    axes[0, 0].set_xticks([])
    axes[1, 0].set_title('Eigenvalue density', fontsize=fsize)
    axes[1, 0].set_ylabel(r'$\rho(\lambda)$', fontsize=fsize - 2, labelpad=0)
    axes[1, 0].set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    axes[1, 0].tick_params(axis='both', which='major', labelsize=fsize - 2)
    axes[1, 0].grid(False)
    axes[1, 0].legend(fontsize=fsize - 3)


def panel_B(ax):
    # GMP distribution: random Gaussian matrix with underlying correlations
    files = ['model_alpha2_sigma0.txt', 'model_alpha2_sigma07.txt',
             'model_alpha2_sigma08.txt', 'model_alpha2_sigma09.txt']
    labels = [r'$\chi=0$ (MP)', r'$\chi=0.7$', r'$\chi=0.8$', r'$\chi=0.9$']
    cmap = plt.cm.RdBu
    colors = [cmap(i) for i in [0.1, 0.7, 0.85, 0.95]]
    for i, file in enumerate(files):
        data = np.loadtxt(os.path.join(root_dir, 'model fit', file))
        ax.plot(data[:, 0], data[:, 1], label=labels[i], color=colors[i], linewidth=1.5)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize, labelpad=0)
    ax.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize, labelpad=0)
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 0.35)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.tick_params(axis='both', which='major', labelsize=fsize)
    ax.set_title('Generalized MP\n(with correlations)', fontsize=fsize)
    ax.legend(fontsize=fsize - 2)


def panel_C(axes):
    # Illustrative sparse simulation matrix + scrambled version with arrows
    np.random.seed(42)
    n_cells, n_genes = 22, 14

    matrix = np.zeros((n_cells, n_genes))
    for j in range(n_genes):
        n_nz = max(1, int(n_cells * np.random.uniform(0.15, 0.38)))
        rows = np.random.choice(n_cells, size=n_nz, replace=False)
        matrix[rows, j] = np.random.exponential(2.5, size=n_nz)

    scrambled = matrix.copy()
    col_perm = {}
    for j in range(n_genes):
        perm = np.random.permutation(n_cells)
        scrambled[:, j] = matrix[perm, j]
        col_perm[j] = perm

    vmax = matrix.max()

    axes[0, 0].imshow(matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=vmax,
                      interpolation='nearest')
    axes[0, 0].set_xticks([])
    axes[0, 0].set_yticks([])
    axes[0, 0].set_xlabel('genes', fontsize=fsize - 2, labelpad=2)
    axes[0, 0].set_ylabel('cells', fontsize=fsize - 2, labelpad=2)
    axes[0, 0].set_title('Simulation data', fontsize=fsize)

    axes[0, 1].imshow(scrambled, aspect='auto', cmap='YlOrRd', vmin=0, vmax=vmax,
                      interpolation='nearest')
    axes[0, 1].set_xticks([])
    axes[0, 1].set_yticks([])
    axes[0, 1].set_xlabel('genes', fontsize=fsize - 2, labelpad=2)
    axes[0, 1].set_title('Scrambled', fontsize=fsize)

    # Arc arrows in the scrambled panel showing cells that moved within columns
    for j in [2, 6, 11]:
        orig_rows = np.where(matrix[:, j] > 0.5)[0]
        if len(orig_rows) == 0:
            continue
        orig_row = orig_rows[0]
        new_row_arr = np.where(col_perm[j] == orig_row)[0]
        if len(new_row_arr) == 0 or abs(new_row_arr[0] - orig_row) < 4:
            continue
        axes[0, 1].annotate('', xy=(j, new_row_arr[0]), xytext=(j, orig_row),
                            arrowprops=dict(arrowstyle='->', color='#2171b5',
                                           lw=1.3, connectionstyle='arc3,rad=0.4'))

    # Schematic arrow between the two panels in figure coordinates
    fig = axes[0, 0].figure
    p0 = axes[0, 0].get_position()
    p1 = axes[0, 1].get_position()
    xmid_start = p0.x1 + 0.002
    xmid_end = p1.x0 - 0.002
    ymid = (p0.y0 + p0.y1) / 2
    fig.add_artist(FancyArrowPatch(
        (xmid_start, ymid), (xmid_end, ymid),
        transform=fig.transFigure,
        arrowstyle='->', mutation_scale=10, color='#636363', lw=1.5))
    fig.text((xmid_start + xmid_end) / 2, ymid + 0.007, 'shuffle',
             transform=fig.transFigure,
             fontsize=fsize - 3, ha='center', va='bottom', color='#636363')


def panel_D(axes):
    # Simulation eigenvalue distributions: original (color-coded) and scrambled
    pcs_data = np.load(os.path.join(ev_data_dir, 'simulated_pcs.npy'))
    pcs, pcs1 = pcs_data[0], pcs_data[1]

    data1 = pcs[pcs > 0]
    data2 = pcs1[pcs1 > 0]
    bin_width = 0.2
    x1 = (1 + np.sqrt(2)) ** 2   # MP upper edge (γ = P/N = 2)
    x2 = float(np.max(pcs1))     # scrambled maximum = GMP-Cor threshold

    all_data = np.concatenate([data1, data2])
    bin_edges = np.arange(min(all_data), max(all_data) + bin_width, bin_width)
    xlim = (bin_edges[0], bin_edges[-1] + bin_width)

    ax_top = axes[0, 0]
    ax_bot = axes[1, 0]

    # Top: original pcs with three-color coding
    _, _, patches1 = ax_top.hist(
        data1, bins=bin_edges, width=bin_width * 0.8, align='right',
        edgecolor='black', color='#d9d9d9', alpha=0.7, density=True)
    for patch in patches1:
        bx = patch.get_x()
        if bx < x1:
            patch.set_facecolor('darkgray')
        elif bx < x2:
            patch.set_facecolor('salmon')
        else:
            patch.set_facecolor('skyblue')

    ax_top.axvline(x1, color='k', linestyle='--', alpha=0.6,
                   label=r'$\lambda^{max}_{MP}$')
    ax_top.axvline(x2, color='dimgray', linestyle=':', alpha=0.8,
                   label=r'$\lambda^{max}_{scr}$')
    legend_handles = [
        Patch(facecolor='darkgray', edgecolor='black', label='MP noise'),
        Patch(facecolor='salmon', edgecolor='black', label='above MP'),
        Patch(facecolor='skyblue', edgecolor='black', label='GMP-Cor'),
    ]
    ax_top.legend(handles=legend_handles, fontsize=fsize - 3,
                  loc='upper right', framealpha=1)
    ax_top.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize)
    ax_top.set_title('Simulation', fontsize=fsize)
    ax_top.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax_top.set_xlim(xlim)
    ax_top.set_xticks([])
    ax_top.grid(False)

    # Bottom: scrambled pcs1 (same x-axis)
    ax_bot.hist(
        data2, bins=bin_edges, width=bin_width * 0.8, align='right',
        edgecolor='black', color='darkgray', alpha=0.7, density=True,
        label='scrambled')
    ax_bot.axvline(x2, color='dimgray', linestyle=':', alpha=0.8,
                   label=r'$\lambda^{max}_{scr}$')
    ax_bot.legend(fontsize=fsize - 3, loc='upper right', framealpha=1)
    ax_bot.set_xlabel(r'$\lambda$', fontsize=fsize, labelpad=0)
    ax_bot.set_ylabel(r'$\rho(\lambda)$', fontsize=fsize)
    ax_bot.set_title('Scrambled', fontsize=fsize - 1)
    ax_bot.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax_bot.set_xlim(xlim)
    ax_bot.grid(False)


def panel_E(ax):
    # Regulated dataset: side-by-side PDF histogram + inset CCDF
    arr = np.load(os.path.join(ev_data_dir, 'sample_13b_filtered.npy'))
    data1 = arr[0, :];  data1 = data1[data1 > 0]
    data2 = arr[1, :];  data2 = data2[data2 > 0]
    x2 = float(np.max(data2))
    bin_width = 0.2
    bin_edges = np.arange(min(np.concatenate([data1, data2])), 10 + bin_width, bin_width)

    _, _, patches1 = ax.hist(
        data1, bins=bin_edges, width=bin_width * 0.5, align='left',
        edgecolor='black', color='#d9d9d9', alpha=0.7, density=True)
    ax.hist(data2, bins=bin_edges + bin_width * 0.5, width=bin_width * 0.5,
            edgecolor='black', color='black', alpha=0.7, density=True,
            align='right', label='scrambled')
    for patch in patches1:
        patch.set_facecolor('skyblue' if patch.get_x() >= x2 else 'darkgray')

    ax.set_xlabel(r'Eigenvalue - $\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel(r'Density - $\rho(\lambda)$', fontsize=fsize - 2, labelpad=0)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xlim([bin_edges[0], bin_edges[-1] + bin_width])
    ax.grid(False)
    ax.set_title('Reg-Arrest (rep. 1)', fontsize=fsize)
    ax.tick_params(axis='both', labelsize=fsize - 2)

    # Inset CCDF
    inset = ax.inset_axes([0.40, 0.38, 0.56, 0.56])
    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    noise = d1s < x2
    inset.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
                 color='darkgray', alpha=0.7, label='noise')
    inset.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
                 color='skyblue', label='signal')
    inset.loglog(d2s, ccdf2, '.', linestyle='-',
                 color='black', alpha=0.5, label='scrambled')
    inset.axvline(x2, color='k', linestyle='--', alpha=0.6)
    inset.set_xlim([0.1, np.max(d1s)])
    inset.set_xlabel(r'$\lambda$', fontsize=fsize - 3)
    inset.set_ylabel('CCDF', fontsize=fsize - 3)
    inset.legend(fontsize=fsize - 4)
    inset.tick_params(labelsize=fsize - 4)


def _plot_ccdf(ax, npy_file, title):
    arr = np.load(os.path.join(ev_data_dir, npy_file))
    data1 = arr[0, :];  data1 = data1[data1 > 0]
    data2 = arr[1, :];  data2 = data2[data2 > 0]
    x2 = float(np.max(data2))
    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    noise = d1s < x2
    ax.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='noise')
    ax.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
              color='skyblue', label='signal')
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled')
    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
    ax.set_xlabel(r'$\lambda$', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)
    ax.set_title(title, fontsize=fsize)
    ax.legend(fontsize=fsize - 3)
    ax.tick_params(labelsize=fsize - 2)


def panel_F(ax):
    _plot_ccdf(ax, 'sample_2b_filtered.npy', 'Exponential')


def panel_G(ax):
    _plot_ccdf(ax, 'sample_15b_filtered.npy', 'Reg-Arrest (rep. 2)')


# ------------------------------------------------------------------
# ASSEMBLE FIGURE
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 9.5), label_offset=(-0.04, 0.04))

panel_pos = [
    [0.04, 0.61, 0.18, 0.35],   # A – random matrix (2×1 grid)
    [0.27, 0.71, 0.21, 0.23],   # B – GMP curves (single)
    [0.52, 0.61, 0.44, 0.35],   # C – sparse simulation (1×2 grid)
    [0.04, 0.05, 0.26, 0.50],   # D – simulation eigenvalues (2×1 grid, tall)
    [0.36, 0.30, 0.59, 0.24],   # E – regulated dataset PDF + CCDF inset (wide)
    [0.36, 0.05, 0.28, 0.21],   # F – CCDF (Exponential)
    [0.69, 0.05, 0.27, 0.21],   # G – CCDF (Reg-Arrest rep. 2)
]

# Panel A
axes_panel_A = pf.add_grid_panel(panel_pos[0], 2, 1, hspace=0.4)
panel_A(axes_panel_A)

# Panel B
pf.add_panel(panel_pos[1], draw_func=panel_B)

# Panel C
axes_panel_C = pf.add_grid_panel(panel_pos[2], 1, 2, wspace=0.12)
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

# pf.save("figure2.svg", dpi=300)
plt.show()
