import os
import sys
import json

# --- import bootstrap: this script lives in scripts/supplementary_figures/ ----
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))                       # repo root
sys.path.insert(0, _REPO)                                             # import src.*
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))         # figure_functions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from figure_functions import PanelFigure

# ------------------------------------------------------------------
# Supplementary Figure S12
#   Subpopulation mixing and subsampling behaviour of GMP-Cor
#   (Reviewer #1, comments 1 and 2.3).
#
#   Row 1 -- mixing two transcriptomically distinct sub-populations RAISES
#   GMP-Cor; it never lowers it, so a low observed value cannot be an artefact
#   of population heterogeneity.
#     A. simulated regulated sub-populations    (rho = 0.7)
#     B. the experimental counterpart: the VapC UMAP of figure5.py panel C,
#        annotated with the matched-n GMP-Cor of Exponential, VapC-2h and their
#        50/50 mixture, read from the dataset_mixing_ratio run log.
#
#   Row 2 -- how GMP-Cor scales with the size of the matrix.
#     C. cells subsampled, gene panel complete: experiment and simulation follow
#        the same curve.
#     D. cells and genes subsampled together at a fixed cell:gene ratio: GMP-Cor
#        is extensive in the gene count while the per-gene index is invariant.
#        Simulation only -- the experimental matrices have no spare gene
#        dimension to subsample at a fixed ratio.
#
#   Every GMP-Cor in row 1 is computed at a matched cell number (n = 1000), and
#   in the simulated panels on a matched gene panel, because the scrambled
#   threshold is a pure function of matrix shape.
# ------------------------------------------------------------------

fsize = 10

# ── Source runs ──────────────────────────────────────────────────────────────

SIM_LOG_DIR = os.path.join(_REPO, 'results', 'simulation_results', 'logs')
REG_LOG = os.path.join(SIM_LOG_DIR, 'inverted_subpopulation_mixing_20260830_105044.json')
GENE_LOG = os.path.join(SIM_LOG_DIR, 'subsampling_robustness_rho09_20260615_114715.json')
CELL_LOG = os.path.join(_REPO, 'results', 'subsampling_experimental', 'logs',
                        'subsampling_experimental_2b_20260830_113604.json')
# Mixture computed on the data_for_paper matrices with the FULL union gene panel and
# ALL cells (no subsampling, no Fano cut), under the paper's own exact-case reporter
# drop list -- so the pure populations reproduce their published per-sample values and
# the mixture sits on the same matrices.
MIX_LOG = os.path.join(SIM_LOG_DIR, 'dataset_mixing_ratio_20260830_131921.json')
MIX_LOG_SUB = os.path.join(SIM_LOG_DIR, 'dataset_mixing_ratio_20260830_122400.json')
# Published per-sample GMP-Cor (sum_denoised_ev), the values behind Fig. 2/3.
DATA_METRICS = os.path.join(_REPO, 'results', 'data_metrics', 'data_metrics.csv')
PUB_FILE_1 = 'Expira_biorep_t0A_filtered.csv'    # Exponential
PUB_FILE_2 = 'VapC_biorep_t2A_filtered.csv'      # VapC: 2h
VAPC_UMAP = os.path.join(_REPO, 'scanpy', 'umap_coordinates_vapc.csv')

# ── Palette (matches figure5.py) ─────────────────────────────────────────────

C_EXP = '#4393c3'      # exponential / sub-population A
C_T2 = '#f4a582'       # VapC 2h
C_T5 = '#d6604d'       # VapC 5h
C_TON = '#b2182b'      # VapC 24h / sub-population B
C_SIM = '#4393c3'
C_DATA = '#b2182b'
BOX = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.72, linewidth=0.8)


def _load(path):
    with open(path, encoding='utf8') as fh:
        return json.load(fh)


def _subsampling_mixing_values(log):
    """
    shared gene panel, n=900 randomly sampled cells, p=2000 highest fano factor of the shared gene panel

    Returns:
    dictionary with values from run "dataset_mixing_ratio_20260830_122400.csv"
    """
    per = pd.DataFrame(log['per_ratio']).set_index('ratio_1')
    return dict(mixture=per.loc[0.5, 'gmp_cor_mean'],
                sd_mixture=per.loc[0.5, 'gmp_cor_std'],
                d_gmp=per.loc[0.5, 'd_gmp_mean'],
                n_mixture=int(per.loc[0.5, 'n_cells_mean']),
                p_mixture=int(per.loc[0.5, 'p_kept_mean']),
                n_repeats=log['params']['repeats'],
                pure_a=per.loc[1.0, 'gmp_cor_mean'],
                pure_b=per.loc[0.0, 'gmp_cor_mean'])


def _experimental_mixing_values(log):
    """
    Panel C annotation values.

    The two pure populations are the PUBLISHED per-sample GMP-Cor from
    results/data_metrics/data_metrics.csv (column `sum_denoised_ev`) -- the same
    numbers reported per sample elsewhere in the manuscript -- so that the panel is
    comparable to those.

    The mixture is computed by simulations/dataset_mixing_ratio_run.py on the same two
    data_for_paper matrices, using the FULL union gene panel and ALL cells (no
    subsampling, no Fano cut) under the paper's own exact-case reporter drop list.

    In that frame the pure populations reproduce their published values, so the three
    numbers are comparable in a way the earlier top-2000-Fano version was not. The
    mixture still has more genes and more cells than either pure point, both of which
    raise GMP-Cor on their own: see documents/figure_s12_details.md, panel C.
    """
    metrics = pd.read_csv(DATA_METRICS).set_index('file_name')['sum_denoised_ev']
    per = pd.DataFrame(log['per_ratio']).set_index('ratio_1')
    return dict(pure_a=float(metrics.loc[PUB_FILE_1]),
                pure_b=float(metrics.loc[PUB_FILE_2]),
                mixture=per.loc[0.5, 'gmp_cor_mean'],
                sd_mixture=per.loc[0.5, 'gmp_cor_std'],
                d_gmp=per.loc[0.5, 'd_gmp_mean'],
                n_mixture=int(per.loc[0.5, 'n_cells_mean']),
                p_mixture=int(per.loc[0.5, 'p_kept_mean']),
                n_repeats=log['params']['repeats'],
                mix_pure_a=per.loc[1.0, 'gmp_cor_mean'],
                mix_pure_b=per.loc[0.0, 'gmp_cor_mean'])


def _mixing_values(log, sigma_mode='shared'):
    """Pure sub-population and 50/50 mixture GMP-Cor, averaged over repeats."""
    df = pd.DataFrame(log['records'])
    df = df[df['sigma_mode'] == sigma_mode]
    df = df[df['ratio_a'] != 'dysregulated_reference']
    df['ratio_a'] = df['ratio_a'].astype(float)
    g = df.groupby('ratio_a')['gmp_cor'].mean()
    # ratio_a is the fraction of cells drawn from sub-population A
    return dict(pure_a=g.loc[1.0], pure_b=g.loc[0.0], mixture=g.loc[0.5])


def _make_room(ax, xy, xpad=0.18, ypad=0.30):
    """
    Widen the axes limits so annotation boxes can sit in the margin rather than on
    top of the points. Returns (x0, x1, y0, y1, xr, yr) of the DATA extent.
    """
    x, y = np.asarray(xy[:, 0]), np.asarray(xy[:, 1])
    x0, x1, y0, y1 = x.min(), x.max(), y.min(), y.max()
    xr, yr = x1 - x0, y1 - y0
    ax.set_xlim(x0 - xpad * xr, x1 + xpad * xr)
    ax.set_ylim(y0 - 0.05 * yr, y1 + ypad * yr)
    return x0, x1, y0, y1, xr, yr


def _annotate_mixture(ax, pt_a, pt_b, mix_pt, vals, label_a, label_b,
                      color_a, color_b):
    """
    Label each cluster with its own GMP-Cor and draw a box above them, joined to
    both by a line, carrying the GMP-Cor of the mixture. All three positions are
    given by the caller in data coordinates, so the boxes can be placed in the
    margin created by `_make_room` instead of over the point cloud.
    """
    for pt in (pt_a, pt_b):
        ax.annotate('', xy=pt, xytext=mix_pt, zorder=2,
                    arrowprops=dict(arrowstyle='-', color='0.45', lw=0.8,
                                    shrinkA=2, shrinkB=2))

    ax.text(*pt_a, f'{label_a}\nGMP-Cor = {vals["pure_a"]:.1f}',
            fontsize=fsize - 3, ha='center', va='center', zorder=3,
            bbox=dict(**BOX, edgecolor=color_a))
    ax.text(*pt_b, f'{label_b}\nGMP-Cor = {vals["pure_b"]:.1f}',
            fontsize=fsize - 3, ha='center', va='center', zorder=3,
            bbox=dict(**BOX, edgecolor=color_b))
    ax.text(*mix_pt, f'mixture\nGMP-Cor = {vals["mixture"]:.1f}',
            fontsize=fsize - 3, ha='center', va='center', zorder=3,
            fontweight='bold', bbox=dict(**BOX, edgecolor='0.35'))


def _bare(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])


# ── Panel A — simulated sub-population mixture ───────────────────────────────

def _simulated_umap(ax, log, title):
    ex = log['example']
    emb = np.asarray(ex['umap'])
    lab = np.asarray(ex['labels'])
    vals = _mixing_values(log)

    for k, color in ((0, C_EXP), (1, C_TON)):
        ax.scatter(emb[lab == k, 0], emb[lab == k, 1], color=color, alpha=.8,
                   s=.5, zorder=1)

    _bare(ax)
    x0, x1, y0, y1, xr, yr = _make_room(ax, emb, xpad=0.30, ypad=0.34)

    # each cluster's box goes in the margin on its own side, clear of the points
    cen = {k: emb[lab == k].mean(axis=0) for k in (0, 1)}
    a_is_left = cen[0][0] < cen[1][0]
    pt_a = ((x0 - 0.16 * xr) if a_is_left else (x1 + 0.16 * xr), cen[0][1])
    pt_b = ((x1 + 0.16 * xr) if a_is_left else (x0 - 0.16 * xr), cen[1][1])
    mix_pt = ((x0 + x1) / 2, y1 + 0.20 * yr)

    _annotate_mixture(ax, pt_a, pt_b, mix_pt, vals, 'population A', 'population B',
                      C_EXP, C_TON)
    ax.set_title(title, fontsize=fsize - 1, pad=2)
    return vals


def panel_A(ax):
    return _simulated_umap(ax, _load(REG_LOG),
                           'Regulated sub-populations ($\\chi$ = 0.7)')


# ── Panel B — experimental VapC UMAP (reproduced from figure5.py panel C) ────

def panel_B(ax):
    data = pd.read_csv(VAPC_UMAP, index_col=0, header=0)
    exp_data = data[data['batch'] == 'exp']
    t2_data = data[data['batch'] == 'T2']
    t5_data = data[np.logical_or(data['batch'] == 'T5A', data['batch'] == 'T5B')]
    ton_data = data[data['batch'] == 'TON']

    # the two mixed populations carry their colors; the other time points stay as a
    # transparent grey backdrop so the annotated clusters read on their own
    for sub in (t5_data, ton_data):
        ax.scatter(sub.UMAP_1, sub.UMAP_2, color='0.6', alpha=.1, s=.5, zorder=1)
    for sub, color in ((exp_data, C_EXP), (t2_data, C_T2)):
        ax.scatter(sub.UMAP_1, sub.UMAP_2, color=color, alpha=.8, s=.5, zorder=2)

    _bare(ax)
    xy = data[['UMAP_1', 'UMAP_2']].to_numpy(dtype=float)
    x0, x1, y0, y1, xr, yr = _make_room(ax, xy, xpad=0.26, ypad=0.34)

    # Exponential sits bottom-left and VapC-2h upper-right in this embedding, so
    # each box goes into the margin on its own side.
    pt_exp = (x0 - 0.14 * xr, exp_data.UMAP_2.mean())
    pt_t2 = (x1 - 0.12 * xr, t2_data.UMAP_2.mean() + 0.10 * yr)
    mix_pt = (x0 + 0.30 * xr, y1 + 0.20 * yr)
    vals = _subsampling_mixing_values(_load(MIX_LOG_SUB))
    _annotate_mixture(ax, pt_exp, pt_t2, mix_pt, vals, 'Exponential', 'VapC: 2h',
                      C_EXP, C_T2)
    ax.set_title('Experimental (VapC)', fontsize=fsize - 1, pad=2)


# ── Panel C — cell subsampling, experiment vs simulation ─────────────────────

def panel_C(ax):
    log = _load(CELL_LOG)
    per_size = pd.DataFrame(log['per_size'])
    ref = log['reference_size']

    arms = [('experimental (sample_2b)', C_DATA, 'o', 'Experimental (regulated)'),
            ('simulated (calibrated, rho=0.7)', C_SIM, 's', 'Simulation ($\\chi$ = 0.7)')]
    for arm, color, marker, label in arms:
        g = per_size[per_size['arm'] == arm].sort_values('n_cells')
        yerr = g['gmp_cor_std'] / g['gmp_cor_mean'] * g['frac_of_reference']
        ax.errorbar(g['n_cells'], g['frac_of_reference'], yerr=yerr, marker=marker,
                    color=color, capsize=2.5, markersize=4, lw=1.2, label=label)

    ax.axhline(1.0, color='0.6', ls=':', lw=0.8)
    ax.set_xlabel('# of Cells, Genes fixed', fontsize=fsize - 2, labelpad=1)
    ax.set_ylabel(f'GMP-Cor / GMP-Cor at n = {ref}', fontsize=fsize - 2, labelpad=2)
    ax.set_title('Cell subsampling', fontsize=fsize - 1, pad=3)
    ax.tick_params(axis='both', labelsize=fsize - 3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=fsize - 3, frameon=False, loc='lower right',
              handlelength=1.4, borderpad=0.1, labelspacing=0.25)


# ── Panel D — gene subsampling at fixed cell:gene ratio (extensivity) ────────

def panel_D(ax):
    log = _load(GENE_LOG)
    per = pd.DataFrame(log['per_size']).sort_values('n_genes')

    ax.errorbar(per['n_genes'], per['mean'], yerr=per['std'], marker='o',
                color=C_SIM, capsize=2.5, markersize=4, lw=1.2, ls='none',
                label='GMP-Cor')
    slope = float(np.sum(per['n_genes'] * per['mean']) / np.sum(per['n_genes'] ** 2))
    xs = np.array([0, per['n_genes'].max() * 1.05])
    ax.plot(xs, slope * xs, color=C_SIM, lw=1.0, ls='--', alpha=0.8,
            label=f'{slope:.4f} $\\times$ genes')
    ax.set_xlim(0, per['n_genes'].max() * 1.08)
    ax.set_ylim(0, per['mean'].max() * 1.25)
    ax.set_xlabel('# of Genes (cell:gene ratio fixed)', fontsize=fsize - 2, labelpad=1)
    ax.set_ylabel('GMP-Cor', fontsize=fsize - 2, labelpad=2, color=C_SIM)
    ax.tick_params(axis='y', labelsize=fsize - 3, colors=C_SIM)
    ax.tick_params(axis='x', labelsize=fsize - 3)
    ax.set_title('GMP-Cor scales with gene number', fontsize=fsize - 1, pad=3)
    ax.spines['top'].set_visible(False)

    ax2 = ax.twinx()
    ax2.errorbar(per['n_genes'], per['mean_per_gene'], yerr=per['std_per_gene'],
                 marker='^', color='0.35', capsize=2.5, markersize=4, lw=1.2,
                 label='GMP-Cor / gene')
    ax2.set_ylabel('GMP-Cor per gene', fontsize=fsize - 2, labelpad=2, color='0.35')
    ax2.set_ylim(0, max(per['mean_per_gene'] + per['std_per_gene']) * 1.9)
    ax2.tick_params(axis='y', labelsize=fsize - 3, colors='0.35')
    ax2.spines['top'].set_visible(False)

    handles = (ax.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0])
    labels = (ax.get_legend_handles_labels()[1] + ax2.get_legend_handles_labels()[1])
    ax.legend(handles, labels, fontsize=fsize - 3, frameon=False, loc='upper left',
              handlelength=1.4, borderpad=0.1, labelspacing=0.25)


# ------------------------------------------------------------------
# BUILD FIGURE
# ------------------------------------------------------------------

pf = PanelFigure(figsize=(9, 7.2), label_offset=(-0.045, 0.035))

pf.add_panel([0.075, 0.575, 0.37, 0.345], label='A', draw_func=panel_A)
pf.add_panel([0.565, 0.575, 0.37, 0.345], label='B', draw_func=panel_B)
pf.add_panel([0.095, 0.085, 0.345, 0.355], label='C', draw_func=panel_C)
pf.add_panel([0.595, 0.085, 0.345, 0.355], label='D', draw_func=panel_D)

pf.fig.text(0.5, 0.535, 'Mixing two distinct sub-populations raises GMP-Cor',
            ha='center', va='center', fontsize=fsize - 1, style='italic', color='0.35')

pf.save(os.path.join(_HERE, 'figure_s12.pdf'), dpi=300)
pf.save(os.path.join(_HERE, 'figure_s12_preview.png'), dpi=200)
print('Saved figure_s12.pdf and figure_s12_preview.png')

for name, log in (('regulated', REG_LOG),):
    v = _mixing_values(_load(log))
    print(f'{name:<13} pure A={v["pure_a"]:7.2f}  pure B={v["pure_b"]:7.2f}  '
          f'mixture={v["mixture"]:7.2f}  ({v["mixture"] / max(v["pure_a"], v["pure_b"]):.2f}x '
          f'the larger pure population)')
v = _subsampling_mixing_values(_load(MIX_LOG_SUB))
print(f'experimental  Exp={v["pure_a"]:7.2f}  VapC-2h={v["pure_b"]:7.2f}   [published per-sample]')
print(f'              mixture={v["mixture"]:7.2f} +/- {v["sd_mixture"]:.2f}  '
      f'(full union panel, all cells: n={v["n_mixture"]}, p={v["p_mixture"]}, '
      f'{v["n_repeats"]} repeats, dGMP={v["d_gmp"]:.2f})')
#print(f'              same-frame endpoints from that run: Exp={v["mix_pure_a"]:.2f}, '
#      f'VapC-2h={v["mix_pure_b"]:.2f}  (both reproduce their published values)')
