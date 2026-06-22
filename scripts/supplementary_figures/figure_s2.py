"""
Supplementary Figure S2 — Synthetic scRNA-seq data generator.

Explains, end to end, how the simulated single-cell data used in the paper is
produced, with schematics + representative example renderings:

  1. Design the gene-gene correlation structure (cluster + hub factor model)
  2. Representative correlation matrices (covariance heatmaps, high vs low rho)
  3. Generate counts via a Gaussian copula (MVN -> CDF -> NB -> library -> dropout)
  4. Representative simulated output (count marginal, heavy tail, sparse matrix)

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
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
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
NB_C = 'steelblue'   # count distributions

# Heatmap generation parameters (match the rho_sweep run behind the Fig. 3 CCDFs)
_HEATMAP_PARAMS = dict(n=2000, shape=1.5, hub_probability=0.2, seed=31)

# ------------------------------------------------------------------
# Pre-computed data (generated once, reused across panels)
# ------------------------------------------------------------------
# Representative correlation matrices for the heatmaps (row 2)
R_high = generate_gram_hub_matrix(alpha=0.9, **_HEATMAP_PARAMS)
R_low = generate_gram_hub_matrix(alpha=0.5, **_HEATMAP_PARAMS)

# Small correlation matrix used to drive the representative count simulation (row 4)
_N_SMALL = 400
R_small = generate_gram_hub_matrix(n=_N_SMALL, alpha=0.9, shape=1.5,
                                   hub_probability=0.2, seed=31)
true_counts, obs_counts = simulate_scRNA_data(
    n_cells=_N_SMALL, n_genes=_N_SMALL, sigma=R_small, dropout_rate=1, inv_gamma_shape=2, inv_gamma_scale=0.1, seed=0)


def _illustrative_loading_matrix(n=40, alpha=0.8, seed=3):
    """Build a small, clean loading matrix A in the same form the generator uses:
    A = sqrt(1-alpha) * I  |  cluster columns  |  one global-hub column.
    Returns A (n x k) and the normalized correlation matrix R (n x n)."""
    rng = np.random.default_rng(seed)
    A = np.eye(n) * np.sqrt(1 - alpha)
    cluster_sizes = [10, 8, 12, 6, 4]
    cols, hubs, idx = [], [], 0
    for size in cluster_sizes:
        if idx >= n:
            break
        size = min(size, n - idx)
        col = np.zeros((n, 1))
        col[idx:idx + size] = np.sqrt(rng.uniform(0.4 * alpha, alpha))
        cols.append(col)
        hubs.append(idx + size // 2)
        idx += size
    global_hub = np.zeros((n, 1))
    for h in hubs:
        if rng.random() < 0.6:
            global_hub[h] = rng.uniform(0.4 * alpha, alpha)
    cols.append(global_hub)
    A = np.hstack([A] + cols)
    C = A @ A.T
    d = np.sqrt(np.diag(C))
    R = C / np.outer(d, d)
    return A, R


def _ccdf(values):
    v = np.sort(values[values > 0])
    p = len(v)
    return v, 1 - np.arange(1, p + 1) / p + 1 / p


# ==================================================================
# Row 1 — Design the correlation structure (schematics)
# ==================================================================
def panel_A(ax):
    """Cluster + hub factor-model network schematic."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    rng = np.random.default_rng(1)

    gx, gy = 0.5, 0.5
    clusters = [(0.17, 0.80), (0.17, 0.20), (0.83, 0.78), (0.83, 0.22)]
    connect = [True, False, True, True]   # subset linked to the global hub
    radius = 0.10

    cluster_hubs = []
    for (cx, cy) in clusters:
        n_pts = 7
        ang = rng.uniform(0, 2 * np.pi, n_pts)
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


def panel_B(ax):
    """A -> AAᵀ -> normalize -> R construction schematic."""
    A, R = _illustrative_loading_matrix()

    ax_A = ax.inset_axes([0.02, 0.10, 0.30, 0.82])
    ax_A.imshow(A, cmap='Reds', aspect='auto', interpolation='nearest')
    ax_A.set_xticks([]); ax_A.set_yticks([])
    ax_A.set_title('loading matrix $A$', fontsize=fsize - 2, pad=4)
    ax_A.set_xlabel('factors', fontsize=fsize - 3, labelpad=1)
    ax_A.set_ylabel('genes', fontsize=fsize - 3, labelpad=1)

    ax_R = ax.inset_axes([0.60, 0.10, 0.36, 0.82])
    im = ax_R.imshow(R, cmap='RdBu_r', vmin=-1, vmax=1, interpolation='nearest')
    ax_R.set_xticks([]); ax_R.set_yticks([])
    ax_R.set_title('positive definite covariance', fontsize=fsize - 2, pad=4)
    cb = ax.figure.colorbar(im, ax=ax_R, fraction=0.046, pad=0.04)
    cb.set_ticks([-1, 0, 1])
    cb.ax.tick_params(labelsize=fsize - 3)

    ax.annotate('', xy=(0.585, 0.5), xytext=(0.345, 0.5), xycoords='axes fraction',
                arrowprops=dict(arrowstyle='-|>', color='0.3', lw=1.6))
    ax.text(0.465, 0.60, r'$AA^{\top}$' + '\nnormalize', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=fsize - 3)
    ax.set_title('Build the correlation matrix', fontsize=fsize, pad=10)


# ==================================================================
# Row 2 — Representative correlation matrices (covariance heatmaps)
# ==================================================================
def _heatmap(ax, R, title):
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
    _heatmap(ax, R_high, 'High correlation\n($\\rho=0.9$)')


def panel_D(ax):
    _heatmap(ax, R_low, 'Low correlation\n($\\rho=0.5$)')


# ==================================================================
# Row 3 — Generate counts via a Gaussian copula (schematic + mini examples)
# ==================================================================
def panel_E(ax):
    rng = np.random.default_rng(7)
    cov = [[1, 0.8], [0.8, 1]]
    g = rng.multivariate_normal([0, 0], cov, size=300)
    u = norm.cdf(g)
    mu, r = 3.0, 0.5
    counts = nbinom.ppf(u[:, 0], r, r / (r + mu))

    w, h, y0 = 0.165, 0.58, 0.10
    xs = [0.025, 0.285, 0.545, 0.805]

    # (i) correlated latent MVN
    a1 = ax.inset_axes([xs[0], y0, w, h])
    a1.scatter(g[:, 0], g[:, 1], s=3, color=NB_C, alpha=0.5, edgecolors='none')
    a1.set_title('correlated\nlatent (MVN)', fontsize=fsize - 3, pad=2)
    a1.set_xticks([]); a1.set_yticks([])
    a1.set_xlabel('gene $i$', fontsize=fsize - 3, labelpad=1)
    a1.set_ylabel('gene $j$', fontsize=fsize - 3, labelpad=1)

    # (ii) uniform copula
    a2 = ax.inset_axes([xs[1], y0, w, h])
    a2.scatter(u[:, 0], u[:, 1], s=3, color=NB_C, alpha=0.5, edgecolors='none')
    a2.set_title('uniform', fontsize=fsize - 3, pad=2)
    a2.set_xticks([0, 1]); a2.set_yticks([0, 1])
    a2.tick_params(labelsize=fsize - 4)

    # (iii) negative-binomial counts
    a3 = ax.inset_axes([xs[2], y0, w, h])
    a3.hist(counts, bins=np.arange(0, counts.max() + 2) - 0.5,
            color=NB_C, edgecolor='white', lw=0.3)
    a3.set_title('NB counts', fontsize=fsize - 3, pad=2)
    a3.set_xlabel('count', fontsize=fsize - 3, labelpad=1)
    a3.set_yticks([])
    a3.tick_params(labelsize=fsize - 4)

    # (iv) sparse observed matrix (reuse the row-4 simulation)
    a4 = ax.inset_axes([xs[3], y0, w, h])
    a4.imshow(obs_counts[:40, :40] > 0, cmap='Greys', aspect='auto',
              interpolation='nearest')
    a4.set_title('+ scale cell total\n+ model dropouts', fontsize=fsize - 3, pad=2)
    a4.set_xticks([]); a4.set_yticks([])
    a4.set_xlabel('genes', fontsize=fsize - 3, labelpad=1)
    a4.set_ylabel('cells', fontsize=fsize - 3, labelpad=1)

    # arrows + operator labels between the mini-axes
    labels = [r'$\Phi$', r'NB$^{-1}$', 'scale,\ndrop']
    for i in range(3):
        x1 = xs[i] + w + 0.02
        x2 = xs[i + 1] - 0.02
        ax.annotate('', xy=(x2, y0 + h / 2), xytext=(x1, y0 + h / 2),
                    xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='-|>', color='0.3', lw=2))
        ax.text((x1 + x2) / 2, y0 + h / 2 + 0.06, labels[i], transform=ax.transAxes,
                ha='center', va='bottom', fontsize=fsize - 3)


# ==================================================================
# Row 4 — Representative simulated output (examples)
# ==================================================================
def panel_F(ax):
    # cell total
    v, ccdf = _ccdf(true_counts.sum(axis=1))
    ax.loglog(v, ccdf, '.', linestyle='-', color=NB_C, markersize=3)
    ax.set_title('Cell distribution', fontsize=fsize)
    ax.set_xlabel('cell total', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)
    ax.tick_params(labelsize=fsize - 2)


def panel_G(ax):
    v, ccdf = _ccdf(true_counts.sum(axis=0))
    ax.loglog(v, ccdf, '.', linestyle='-', color=NB_C, markersize=3)
    ax.set_title('Gene distribution', fontsize=fsize)
    ax.set_xlabel('gene total', fontsize=fsize - 2, labelpad=0)
    ax.set_ylabel('CCDF', fontsize=fsize - 2, labelpad=0)
    ax.tick_params(labelsize=fsize - 2)


def panel_H(ax):
    sub = obs_counts[:120, :120].astype(float)
    vmax = np.percentile(sub[sub > 0], 95) if np.any(sub > 0) else 1
    ax.imshow(sub, cmap='Greys', vmin=0, vmax=vmax, aspect='auto',
              interpolation='nearest')
    zero_frac = np.mean(obs_counts == 0)
    ax.set_title('Observed counts (sparse)', fontsize=fsize)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('genes', fontsize=fsize - 2, labelpad=1)
    ax.set_ylabel('cells', fontsize=fsize - 2, labelpad=1)
    ax.text(0.97, 0.04, f'{zero_frac * 100:.0f}% zeros', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=fsize - 3, color='black',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='0.7', lw=0.5))


# ------------------------------------------------------------------
# Assemble — portrait, 4-row narrative
# ------------------------------------------------------------------
pf = PanelFigure(figsize=(7, 7.5), label_offset=(-0.02, 0.02))

panel_pos = [
    # Row 1 — schematics
    [0.05, 0.78, 0.42, 0.14],   # A — factor-model network
    [0.5, 0.78, 0.42, 0.14],   # B — A -> R construction
    # Row 2 — correlation heatmaps
    [0.10, 0.54, 0.2, 0.14],   # C — heatmap rho=0.9
    [0.5, 0.54, 0.2, 0.14],   # D — heatmap rho=0.5
    # Row 3 — copula flow
    [0.04, 0.30, 0.93, 0.15],   # E — count-generation flow
    # Row 4 — outputs
    [0.09, 0.05, 0.21, 0.16],   # F — NB count histogram
    [0.41, 0.05, 0.21, 0.16],   # G — gene-total CCDF
    [0.71, 0.05, 0.25, 0.16],   # H — sparse observed matrix
]

pf.add_panel(panel_pos[0], draw_func=panel_A, hide_axis=True, label='A')
pf.add_panel(panel_pos[1], draw_func=panel_B, hide_axis=True, label='B')
pf.add_panel(panel_pos[2], draw_func=panel_C, label='C')
pf.add_panel(panel_pos[3], draw_func=panel_D, label='D')
pf.add_panel(panel_pos[4], draw_func=panel_E, hide_axis=True, label='E')
pf.add_panel(panel_pos[5], draw_func=panel_F, label='F')
pf.add_panel(panel_pos[6], draw_func=panel_G, label='G')
pf.add_panel(panel_pos[7], draw_func=panel_H, label='H')

# Step headers down the left margin
for y, txt in [(0.965, '1 · Design correlation structure'),
               (0.735, '2 · Representative correlation matrices'),
               (0.48, '3 · Generate counts (Gaussian copula)'),
               (0.24, '4 · Representative simulated output')]:
    pf.fig.text(0.02, y, txt, fontsize=fsize - 1, fontweight='bold', color='0.25',
                ha='left', va='bottom')

pf.save("figure_s2.pdf", dpi=300, transparent=True)
plt.show()
