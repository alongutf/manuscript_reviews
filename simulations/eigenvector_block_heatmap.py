"""
Visual check: do the leading eigenvector modes line up with the true correlation blocks?

Plots the ground-truth correlation matrix R next to |v_k| for the 5 leading modes of the
post-dropout counts (via analysis_functions.get_eig_vectors), on a SHARED gene axis, so a
block in R can be read straight across into the mode panel.

Genes are in generator order, and the generator's clusters are contiguous, so R is
block-diagonal as drawn -- no reordering is applied anywhere.

Parameters match simulations/eigenvector_block_recovery_run.py (alpha=0.75,
inv_gamma_scale=0.03, the depth-matched condition), so the panels correspond to the
numbers in that run's log.

Two variants of the mode panel are written:
  abs     |v_k| on a sequential scale -- shows only WHERE a mode puts its weight
  signed  v_k on the same diverging blue-red scale as R -- also shows the sign
          structure, i.e. genes that move together vs. in opposition within a mode.
          Note the overall sign of an eigenvector is arbitrary (v and -v are the same
          mode), so only sign differences WITHIN a panel column are meaningful.

Outputs -> results/simulation_results/figures/
  eigenvector_block_heatmap_<abs|signed>_<ts>.svg / .png       full 2000-gene range
  eigenvector_block_heatmap_zoom_<abs|signed>_<ts>.svg / .png  one row per strongest block
"""

import sys
import os
import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.simulations import generate_gram_hub_matrix, simulate_scRNA_data  # noqa: E402
from src.analysis_functions import get_eig_vectors  # noqa: E402

# ── Parameters (mirror the depth-matched recovery run) ───────────────────────

PARAMS = dict(
    n_genes=2000, n_cells=1000, alpha=0.75, shape=1.5, hub_probability=0.2,
    sigma_seed=31, count_seed=0, dropout_rate=1.0,
    inv_gamma_shape=1.5, inv_gamma_scale=0.03, n_top=5,
)
N_ZOOM = 5        # strongest blocks to zoom into
ZOOM_PAD = 1.0    # zoom window padding, as a multiple of the block's size

_FIG_DIR = os.path.join(_REPO_ROOT, 'results', 'simulation_results', 'figures')
os.makedirs(_FIG_DIR, exist_ok=True)
_ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

# ── Simulate ─────────────────────────────────────────────────────────────────

p = PARAMS
R, structure = generate_gram_hub_matrix(
    p['n_genes'], p['alpha'], p['shape'], p['hub_probability'],
    seed=p['sigma_seed'], return_structure=True)
labels, sizes, strengths = structure['labels'], structure['sizes'], structure['strengths']

print('Simulating counts...')
_, observed = simulate_scRNA_data(
    n_cells=p['n_cells'], n_genes=p['n_genes'], sigma=R,
    dropout_rate=p['dropout_rate'], seed=p['count_seed'],
    inv_gamma_shape=p['inv_gamma_shape'], inv_gamma_scale=p['inv_gamma_scale'])

_, vecs_obs, threshold, kept = get_eig_vectors(observed, n_top=p['n_top'],
                                               norm_method='sum', norm_sum=50)

# Re-index the modes onto the FULL gene axis: get_eig_vectors drops all-zero genes, so
# without this the mode panel would be off by the dropped columns against R and the
# alignment the figure is meant to show would be false.
V_signed = np.zeros((p['n_top'], p['n_genes']))
V_signed[:, np.asarray(kept, dtype=bool)] = vecs_obs
print(f'{int(np.sum(kept))} of {p["n_genes"]} genes kept by the pipeline; '
      f'dropped genes are drawn as zero loading')

# blocks worth annotating: strongest by size * strength
rank = np.argsort(sizes * strengths)[::-1][:N_ZOOM]
starts = {b: int(np.flatnonzero(labels == b)[0]) for b in rank}
ends = {b: int(np.flatnonzero(labels == b)[-1]) + 1 for b in rank}

# ── Shared colour scales ─────────────────────────────────────────────────────

off = R[~np.eye(p['n_genes'], dtype=bool)]
r_vmax = float(np.percentile(np.abs(off), 99.9))
v_vmax = float(np.percentile(np.abs(V_signed), 99.9))
print(f'colour scales: |R| off-diagonal vmax={r_vmax:.3f}, |v| vmax={v_vmax:.4f}')

# the two mode-panel variants: (name, matrix, cmap, vmin, panel title)
VARIANTS = [
    ('abs', np.abs(V_signed), 'magma', 0.0, 'leading modes, |v|'),
    ('signed', V_signed, 'RdBu_r', -v_vmax, 'leading modes, v (signed)'),
]


def draw(ax_R, ax_b, lo, hi, V, cmap, vmin, panel_title, title=None, annotate=True):
    """Draw R and the mode panel over gene range [lo, hi) on a shared axis."""
    ax_R.imshow(R[lo:hi, lo:hi], cmap='RdBu_r', vmin=-r_vmax, vmax=r_vmax,
                aspect='auto', interpolation='nearest',
                extent=[lo, hi, hi, lo])
    ax_R.set_ylabel('gene index')
    if title:
        ax_R.set_title(title, loc='left', fontsize=9, fontweight='bold')
    ax_b.imshow(V[:, lo:hi].T, cmap=cmap, vmin=vmin, vmax=v_vmax,
                aspect='auto', interpolation='nearest',
                extent=[0.5, p['n_top'] + 0.5, hi, lo])
    ax_b.set_xticks(np.arange(1, p['n_top'] + 1))
    ax_b.set_xlabel('mode')
    ax_b.set_title(panel_title, fontsize=8)
    ax_b.set_yticklabels([])
    # block boundaries, drawn across both panels so the eye can carry them over
    if annotate:
        for b in rank:
            s, e = starts[b], ends[b]
            if e < lo or s > hi:
                continue
            for ax in (ax_R, ax_b):
                for y in (s, e):
                    ax.axhline(y, color='0.15', lw=0.6, ls=':', alpha=0.8)
            ax_b.annotate(f'#{b} ({sizes[b]}g)', xy=(p['n_top'] + 0.7, (s + e) / 2),
                          xycoords=('data', 'data'), fontsize=6.5,
                          va='center', ha='left', annotation_clip=False)


# ── Figures: one pair (full range + zoom) per mode-panel variant ─────────────

for vname, V, cmap, vmin, panel_title in VARIANTS:
    fig_main = os.path.join(_FIG_DIR, f'eigenvector_block_heatmap_{vname}_{_ts}')
    fig_zoom = os.path.join(_FIG_DIR, f'eigenvector_block_heatmap_zoom_{vname}_{_ts}')

    # full range
    fig = plt.figure(figsize=(11, 7.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[5, 1.3], wspace=0.06)
    axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
    draw(axes[0], axes[1], 0, p['n_genes'], V, cmap, vmin, panel_title,
         title=f'ground-truth correlation matrix R  (alpha = {p["alpha"]}, '
               f'{len(sizes)} blocks)')
    axes[0].set_xlabel('gene index')
    fig.suptitle('Do the leading modes sit on the true correlation blocks?  '
                 'Shared gene axis, generator order (no reordering)',
                 fontsize=11, fontweight='bold')
    note = (f'R on a +/-{r_vmax:.2f} diverging scale (off-diagonal 99.9th pct); '
            f'modes on a {"+/-" if vmin < 0 else "0-"}{v_vmax:.3f} scale. '
            f'Dotted lines mark the {N_ZOOM} strongest blocks.')
    if vmin < 0:
        note += ('  The overall sign of an eigenvector is arbitrary, so only sign '
                 'differences within one mode are meaningful.')
    fig.text(0.5, 0.015, note, ha='center', fontsize=7.5, color='0.3')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(fig_main + '.svg', format='svg', bbox_inches='tight')
    fig.savefig(fig_main + '.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # zoom: one row per strongest block
    fig, axarr = plt.subplots(N_ZOOM, 2, figsize=(8.0, 3.0 * N_ZOOM),
                              gridspec_kw=dict(width_ratios=[5, 1.3], wspace=0.08,
                                               hspace=0.35))
    for row, b in enumerate(rank):
        st, en = starts[b], ends[b]
        pad = int(max(10, ZOOM_PAD * (en - st)))
        lo, hi = max(0, st - pad), min(p['n_genes'], en + pad)
        draw(axarr[row, 0], axarr[row, 1], lo, hi, V, cmap, vmin, panel_title,
             title=f'block #{b}: {sizes[b]} genes, w = {strengths[b]:.2f}  '
                   f'(genes {st}-{en - 1})')
        axarr[row, 0].set_xlabel('gene index')
    fig.suptitle('Zoom on the strongest blocks — each block and the modes '
                 'over the same genes', fontsize=11, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(fig_zoom + '.svg', format='svg', bbox_inches='tight')
    fig.savefig(fig_zoom + '.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    print(f'{vname:7s} full range : {fig_main}.svg / .png')
    print(f'{vname:7s} zoom       : {fig_zoom}.svg / .png')
