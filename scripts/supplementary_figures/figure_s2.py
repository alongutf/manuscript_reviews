"""
Supplementary Figure S2 — Synthetic scRNA-seq data generator.

Explains, end to end, how the simulated single-cell data used in the paper is
produced, with schematics + representative example renderings:

  1. Design the gene-gene correlation structure (cluster + hub factor model)
     and representative correlation matrices (covariance heatmaps, high vs low chi)
  2. Generate counts (MVN -> CDF -> NB -> library -> dropout)
  3. Representative simulated output (cell / gene rank plots) compared against
     real data (sample_2b_filtered.csv)

Run from this directory:
    cd scripts/supplementary_figures
    python figure_s2.py
The figure is written next to this script as figure_s2.pdf.
"""
import os
import sys

# --- import bootstrap: this script lives in scripts/supplementary_figures/ ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                       # repo root
sys.path.insert(0, _REPO)                                            # import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))        # figure_functions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm, nbinom

from figure_functions import PanelFigure
from src.simulations import generate_gram_hub_matrix, simulate_scRNA_data

# ------------------------------------------------------------------
# Global style
# ------------------------------------------------------------------
fsize = 10
plt.close("all")

GENE_C = '#9ecae1'    # ordinary gene node
HUB_C = '#E07B54'    # cluster hub
GHUB_C = '#c0392b'    # global hub
NB_C = 'steelblue'   # count distributions (flow-chart insets)
RANK_C = '#404040'   # dark grey — all rank plots

# Heatmap generation parameters (match the chi_sweep run behind the Fig. 3 CCDFs)
_HEATMAP_PARAMS = dict(n=2000, shape=1.5, hub_probability=0.2, seed=31)

# Representative real dataset used for the simulated-vs-data rank comparison
_REAL_DATA_PATH = os.path.join(_REPO, 'data_for_paper', 'sample_2b_filtered.csv')

# ------------------------------------------------------------------
# Pre-computed data (generated once, reused across panels)
# ------------------------------------------------------------------
# Representative correlation matrices for the heatmaps (row 1): same
# cluster/hub topology (fixed seed), only the shared-variance fraction alpha
# differs, so the two heatmaps isolate the effect of correlation strength.
R_high = generate_gram_hub_matrix(alpha=0.9, **_HEATMAP_PARAMS)
R_low = generate_gram_hub_matrix(alpha=0.5, **_HEATMAP_PARAMS)

# Small correlation matrix used to drive the count-generation schematic (row 2 flow inset).
# Kept deliberately small (400x400, not the paper-scale 2000x2000) since only a
# handful of cells/genes are actually rendered in the schematic's mini-panels.
_N_SMALL = 400
R_small = generate_gram_hub_matrix(n=_N_SMALL, alpha=0.9, shape=1.5,
                                   hub_probability=0.2, seed=31)
# inv_gamma_shape/scale here are chosen for a visually clear NB histogram in the
# schematic, not to match real data (contrast with the defaults used below for
# rank_obs_counts, which are tuned to reproduce the real marginals).
true_counts, obs_counts = simulate_scRNA_data(
    n_cells=_N_SMALL, n_genes=_N_SMALL, sigma=R_small, dropout_rate=1, inv_gamma_shape=2, inv_gamma_scale=0.1, seed=0)

# Larger simulation, shaped to match the real data (1000 cells x 2000 genes), for the
# simulated-vs-data rank plots. Parameters mirror the rho_sweep GMP-Cor calibration
# (scripts/simulated_data.ipynb): dropout_rate=1, seed=31, and the default count model
# (inv_gamma_shape=1.5, inv_gamma_scale=0.01) so the marginals track the real data.
# R_high is exactly the sweep's sigma at rho=0.9 (generate_gram_hub_matrix(2000, 0.9, 1.5, 0.2, seed=31)).
_N_CELLS_RANK, _N_GENES_RANK = 1000, 2000
_, rank_obs_counts = simulate_scRNA_data(
    n_cells=_N_CELLS_RANK, n_genes=_N_GENES_RANK, sigma=R_high,
    dropout_rate=1, seed=31)
sim_cell_totals = rank_obs_counts.sum(axis=1)
sim_gene_totals = rank_obs_counts.sum(axis=0)

# Real data (cells x genes); totals drive the comparison rank plots
_real = pd.read_csv(_REAL_DATA_PATH, index_col=0).fillna(0).to_numpy(dtype=float)
real_cell_totals = _real.sum(axis=1)
real_gene_totals = _real.sum(axis=0)


def _rank(values, drop_top=False):
    """Return (ranks, sorted_values) for a rank plot: values sorted descending.
    If drop_top, omit the single highest value (an outlier that stretches the scale),
    keeping the original rank index (so ranks start at 2)."""
    v = np.sort(values[values > 0])[::-1]
    ranks = np.arange(1, len(v) + 1)
    if drop_top and len(v) > 1:
        v, ranks = v[1:], ranks[1:]
    return ranks, v


def _pair_ylim(a, b, drop_top=False, pad=1.6):
    """Shared y-limits spanning both series of a simulated-vs-data pair."""
    all_v = np.concatenate([_rank(a, drop_top)[1], _rank(b, drop_top)[1]])
    return all_v.min() / pad, all_v.max() * pad


_CELL_YLIM = _pair_ylim(sim_cell_totals, real_cell_totals)
_GENE_YLIM = _pair_ylim(sim_gene_totals, real_gene_totals, drop_top=True)


# ==================================================================
# Row 1 — Design the correlation structure (schematic + heatmaps)
# ==================================================================
def panel_A(ax):
    """Cluster + hub factor-model network schematic."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    rng = np.random.default_rng(1)

    gx, gy = 0.5, 0.5
    clusters = [(0.17, 0.80), (0.17, 0.20), (0.83, 0.78), (0.83, 0.22)]
    connect = [True, False, True, True]   # subset linked to the global hub
    radius = 0.10   # plotting radius of each cluster's gene scatter, in axes units

    cluster_hubs = []
    for (cx, cy) in clusters:
        n_pts = 7   # genes drawn per cluster, purely for the schematic's visual density
        ang = rng.uniform(0, 2 * np.pi, n_pts)
        # sqrt of a uniform draw gives points uniform over the disc area, not just the radius
        rr = radius * np.sqrt(rng.uniform(0.15, 1, n_pts))
        xs, ys = cx + rr * np.cos(ang), cy + rr * np.sin(ang)
        for i in range(n_pts):                              # star edges hub->members
            ax.plot([cx, xs[i]], [cy, ys[i]], color='0.82', lw=0.5, zorder=1)
        ax.scatter(xs, ys, s=16, color=GENE_C, edgecolors='white', lw=0.3, zorder=2)
        ax.scatter([cx], [cy], s=64, color=HUB_C, edgecolors='white', lw=0.6, zorder=3)
        cluster_hubs.append((cx, cy))

    for (cx, cy), linked in zip(cluster_hubs, connect):     # cluster hub -> global hub
        if linked:
            ax.plot([gx, cx], [gy, cy], color=GHUB_C, lw=1.3, alpha=0.7, zorder=1)
    ax.scatter([gx], [gy], s=180, marker='*', color=GHUB_C,
               edgecolors='white', lw=0.6, zorder=4)

    # annotations
    ax.annotate('gene', xy=(clusters[1][0] + 0.07, clusters[1][1] + 0.05),
                xytext=(0.30, 0.06), fontsize=fsize - 3, ha='left',
                arrowprops=dict(arrowstyle='->', color='0.4', lw=0.7))
    ax.annotate('cluster hub', xy=clusters[2], xytext=(0.55, 0.88),
                fontsize=fsize - 3, ha='left', color=HUB_C,
                arrowprops=dict(arrowstyle='->', color=HUB_C, lw=0.7))
    ax.annotate('global hub', xy=(gx, gy), xytext=(0.50, 0.30),
                fontsize=fsize - 3, ha='center', color=GHUB_C,
                arrowprops=dict(arrowstyle='->', color=GHUB_C, lw=0.7))
    ax.set_title('Cluster + hub factor model', fontsize=fsize)


def _heatmap(ax, R, title):
    """Draw the top-left 100x100 sub-block of a genes x genes correlation matrix R."""
    # only a corner of the full 2000x2000 matrix is shown; large enough to see the
    # block/hub structure without the individual cells becoming illegibly small
    im = ax.imshow(R[:100, :100], cmap='RdBu_r', vmin=-1, vmax=1,
                   aspect='auto', interpolation='nearest')
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_ticks([-1, 0, 1])
    cb.ax.tick_params(labelsize=fsize - 2)
    ax.set_title(title, fontsize=fsize)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('genes', fontsize=fsize - 2, labelpad=1)
    ax.set_ylabel('genes', fontsize=fsize - 2, labelpad=1)


def panel_C(ax):
    _heatmap(ax, R_high, 'High correlation\n($\\chi=0.9$)')


def panel_D(ax):
    _heatmap(ax, R_low, 'Low correlation\n($\\chi=0.5$)')


# ==================================================================
# Row 2 — Generate counts (schematic + mini examples)
# ==================================================================
def panel_E(ax):
    """Flow-chart schematic of the count-generation step: correlated MVN -> NB counts
    -> observed (scaled + dropout) counts. Uses its own small 2-gene toy example,
    independent of R_small/obs_counts above, purely to illustrate the mechanism."""
    rng = np.random.default_rng(7)
    cov = [[1, 0.8], [0.8, 1]]   # toy 2-gene correlation, just for the scatter inset
    g = rng.multivariate_normal([0, 0], cov, size=300)
    u = norm.cdf(g)              # Gaussian copula: correlated normal -> correlated uniform
    mu, r = 3.0, 0.5              # illustrative NB mean/dispersion for the histogram inset
    counts = nbinom.ppf(u[:, 0], r, r / (r + mu))

    w, h, y0 = 0.22, 0.60, 0.08
    xs = [0.04, 0.39, 0.74]

    # (i) correlated latent MVN
    a1 = ax.inset_axes([xs[0], y0, w, h])
    a1.scatter(g[:, 0], g[:, 1], s=3, color=NB_C, alpha=0.5, edgecolors='none')
    a1.set_title('correlated\nlatent (MVN)', fontsize=fsize - 3, pad=2)
    a1.set_xticks([]); a1.set_yticks([])
    a1.set_xlabel('gene $i$', fontsize=fsize - 3, labelpad=1)
    a1.set_ylabel('gene $j$', fontsize=fsize - 3, labelpad=1)

    # (ii) negative-binomial counts
    a2 = ax.inset_axes([xs[1], y0, w, h])
    a2.hist(counts, bins=np.arange(0, counts.max() + 2) - 0.5,
            color=NB_C, edgecolor='white', lw=0.3)
    a2.set_title('NB counts', fontsize=fsize - 3, pad=2)
    a2.set_xlabel('count', fontsize=fsize - 3, labelpad=1)
    a2.set_yticks([])
    a2.tick_params(labelsize=fsize - 4)

    # (iii) observed counts (scaled + dropout) — the real generator output
    a3 = ax.inset_axes([xs[2], y0, w, h])
    sub = obs_counts[:40, :40].astype(float)     # small corner, just for a legible thumbnail
    # clip the color scale below the max so a few very high counts don't wash out
    # the rest of the (mostly small/zero) matrix
    vmax = np.percentile(sub[sub > 0], 95) if np.any(sub > 0) else 1
    a3.imshow(sub, cmap='Greys', vmin=0, vmax=vmax, aspect='auto',
              interpolation='nearest')
    a3.set_title('observed counts', fontsize=fsize - 3, pad=2)
    a3.set_xticks([]); a3.set_yticks([])
    a3.set_xlabel('genes', fontsize=fsize - 3, labelpad=1)
    a3.set_ylabel('cells', fontsize=fsize - 3, labelpad=1)
    zero_frac = np.mean(obs_counts == 0)
    a3.text(0.97, 0.04, f'{zero_frac * 100:.0f}% zeros', transform=a3.transAxes,
            ha='right', va='bottom', fontsize=fsize - 4, color='black',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='0.7', lw=0.5))

    # arrows + operator labels between the mini-axes
    labels = [r'$\Phi,\ \mathrm{NB}^{-1}$', 'scale,\ndrop']
    for i in range(2):
        x1 = xs[i] + w + 0.02
        x2 = xs[i + 1] - 0.02
        ax.annotate('', xy=(x2, y0 + h / 2), xytext=(x1, y0 + h / 2),
                    xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='-|>', color='0.3', lw=2))
        ax.text((x1 + x2) / 2, y0 + h / 2 + 0.06, labels[i], transform=ax.transAxes,
                ha='center', va='bottom', fontsize=fsize - 3)


# ==================================================================
# Row 3 — Representative simulated output vs real data (rank plots)
# ==================================================================
def _rank_plot(ax, totals, title, xlabel, show_ylabel=True, drop_top=False, ylim=None):
    """Semilog-y plot of totals sorted descending against their rank (1 = largest)."""
    ranks, v = _rank(totals, drop_top=drop_top)
    ax.semilogy(ranks, v, '-', color=RANK_C, lw=1)
    ax.set_title(title, fontsize=fsize)
    ax.set_xlabel(xlabel, fontsize=fsize - 2, labelpad=0)
    if show_ylabel:
        ax.set_ylabel('total expression', fontsize=fsize - 2, labelpad=0)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.tick_params(labelsize=fsize - 2)
    if not show_ylabel:                 # matched y-axis within the pair -> drop duplicate ticks
        ax.tick_params(labelleft=False)


def panel_F(ax):      # simulated — cells
    _rank_plot(ax, sim_cell_totals, 'Simulated', 'cell rank', ylim=_CELL_YLIM)


def panel_F_data(ax):  # real data — cells
    _rank_plot(ax, real_cell_totals, 'Data', 'cell rank', show_ylabel=False, ylim=_CELL_YLIM)


def panel_G(ax):      # simulated — genes
    _rank_plot(ax, sim_gene_totals, 'Simulated', 'gene rank', drop_top=True, ylim=_GENE_YLIM)


def panel_G_data(ax):  # real data — genes
    _rank_plot(ax, real_gene_totals, 'Data', 'gene rank', show_ylabel=False,
               drop_top=True, ylim=_GENE_YLIM)


# ------------------------------------------------------------------
# Assemble — portrait, 3-row narrative
# ------------------------------------------------------------------
# panel_pos entries are [left, bottom, width, height] in normalized figure
# coordinates, one per PanelFigure.add_panel call below, in the same order.
pf = PanelFigure(figsize=(7, 6.5), label_offset=(-0.02, 0.02))

panel_pos = [
    # Row 1 — design: network schematic + two correlation heatmaps
    [0.02, 0.64, 0.30, 0.24],   # A — factor-model network
    [0.41, 0.66, 0.23, 0.20],   # C — heatmap chi=0.9
    [0.72, 0.66, 0.23, 0.20],   # D — heatmap chi=0.5
    # Row 2 — count-generation flow
    [0.05, 0.36, 0.90, 0.17],   # E — count-generation flow
    # Row 3 — rank plots: simulated vs data, paired (cells | genes)
    [0.09, 0.07, 0.16, 0.15],   # F  — simulated cell rank
    [0.28, 0.07, 0.16, 0.15],   # F' — data cell rank
    [0.57, 0.07, 0.16, 0.15],   # G  — simulated gene rank
    [0.76, 0.07, 0.16, 0.15],   # G' — data gene rank
]

pf.add_panel(panel_pos[0], draw_func=panel_A, hide_axis=True, label='A')
pf.add_panel(panel_pos[1], draw_func=panel_C, label='B')
pf.add_panel(panel_pos[2], draw_func=panel_D, label='C')
pf.add_panel(panel_pos[3], draw_func=panel_E, hide_axis=True, label='D')
pf.add_panel(panel_pos[4], draw_func=panel_F, label='E')
pf.add_panel(panel_pos[5], draw_func=panel_F_data, label=' ')
pf.add_panel(panel_pos[6], draw_func=panel_G, label='F')
pf.add_panel(panel_pos[7], draw_func=panel_G_data, label=' ')

# Step headers down the left margin
for y, txt in [(0.925, '1 · Design correlation structure'),
               (0.575, '2 · Generate counts'),
               (0.255, '3 · Representative output vs. data')]:
    pf.fig.text(0.02, y, txt, fontsize=fsize - 1, fontweight='bold', color='0.25',
                ha='left', va='bottom')

pf.save("figure_s2_preview.png", dpi=200)
pf.save("figure_s2.pdf", dpi=300, transparent=True)
plt.show()
