"""
Supplementary Figure S3 — scRNA-seq / bulk RNA-seq concordance (Reviewer #1, Comment 3).

Panel layout (9 × 8 in):
  A  Three stacked mean-expression scatter plots (sc vs bulk, log-log scale).
     Top: exponential (control); middle: Dis-Arrest; bottom: Reg-Arrest.
  B  Fold-change scatter: sc pseudobulk DESeq2 vs bulk DESeq2 (Dis-Arrest vs ctrl).
  C  Fold-change scatter: sc pseudobulk DESeq2 vs bulk DESeq2 (Reg-Arrest vs ctrl).
  D  GO terms significant in both conditions: paired -log10(FDR) bar chart.

Run from any directory; uses absolute paths derived from __file__.
"""
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('default')
from scipy.stats import pearsonr, spearmanr

_HERE    = os.path.dirname(os.path.abspath(__file__))
_REPO    = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_REPO, 'scripts', 'figures'))

from figure_functions import PanelFigure

META       = os.path.join(_REPO, 'metadata')
DATA       = os.path.join(_REPO, 'data_for_paper')
BULK_XLSX  = os.path.join(_REPO, 'bulk_data', 'All bulk data.xlsx')
BULK2_CSV = os.path.join(_REPO, 'bulk_data', 'exp0224_molecules_per_cell.csv')
BULK_DESEQ = os.path.join(_REPO, 'results', 'deseq_results', 'from counts')
SC_DESEQ   = os.path.join(_REPO, 'results', 'deseq_results', 'sc_pseudobulk')
GO_TABLE   = os.path.join(_REPO, 'results', 'GO_results', 'sc_pseudobulk',
                          'sc_GO_FDR_comparison_down_N500_table.csv')

SPIKE_INS = {'gfp', 'mcherry', 'tetr', 'laci', 'ampr', 'lelobekk'}


def _fmt_p(p):
    """Format a p-value rounded to the nearest decade, e.g. p = 10⁻¹⁶ (mathtext)."""
    if not np.isfinite(p) or p <= 0:
        return r'$p < 10^{-300}$'
    exp = int(round(np.log10(p)))
    return rf'$p = 10^{{{exp}}}$'


# rounded annotation box with a transparent fill (grid shows through)
_ANNO_BBOX = dict(boxstyle='round,pad=0.35', facecolor='none', edgecolor='0.6',
                  linewidth=0.6)


def _style_ax(ax, grid=True):
    """White background, full black box, optional light grid."""
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(0.8)
    if grid:
        ax.grid(True, color='lightgrey', linewidth=0.5, zorder=0)
    else:
        ax.grid(False)


# ── gene-name helpers (inlined from sc_pseudobulk_go / sc_bulk_correlation) ──

def get_gene_synonyms():
    gtf = pd.read_csv(os.path.join(META, 'genomic.gtf'), sep='\t',
                      comment='#', header=None)
    syn = {}
    for i in range(len(gtf)):
        if gtf.iloc[i, 2] == 'gene':
            attr = gtf.iloc[i, 8]
            m = re.search(r'gene "([^"]+)"', attr)
            if m:
                primary = m.group(1).lower()
                for s in re.findall(r'gene_synonym "([^"]+)"', attr):
                    syn[s.lower()] = primary
    return syn


def harmonize(genes, syn):
    out = []
    for g in genes:
        g = g.replace('LELOBEKK_', '').replace('LELOBEKK', '').lower()
        out.append(syn.get(g, g))
    return out


def map_locus_to_name():
    gff = pd.read_csv(os.path.join(META, 'gffFile_combined.gff'),
                      sep='\t', comment='#', header=None)
    out = {}
    for i in range(len(gff)):
        if gff.iloc[i, 2] == 'gene':
            attr = gff.iloc[i, 8]
            nm = re.search(r'Name=([^;]+)', attr)
            lt = re.search(r'locus_tag=([^;]+)', attr)
            if lt and nm:
                out[lt.group(1)] = nm.group(1)
    return out


def harmonize_bulk(genes, syn, locus2name):
    mapped = [locus2name.get(g, g) for g in genes]
    return harmonize(mapped, syn)


# ── panel functions ────────────────────────────────────────────────────────────

def mean_expression_corr(syn, locus2name, bulk_df, sc_file, bulk_cols_prefix, title, ax):
    sc = pd.read_csv(sc_file, index_col=0)
    sc_mean = sc.mean(axis=0)
    sc_mean.index = harmonize(list(sc_mean.index), syn)
    sc_mean = sc_mean[~sc_mean.index.isin(SPIKE_INS)].groupby(level=0).sum()

    bulk_cols = [c for c in bulk_df.columns if bulk_cols_prefix in str(c).lower()]
    bulk_mean = bulk_df[bulk_cols].mean(axis=1)
    bulk_mean.index = harmonize_bulk(list(bulk_mean.index), syn, locus2name)
    bulk_mean = bulk_mean.groupby(level=0).sum()

    common    = sorted(set(sc_mean.index) & set(bulk_mean.index))
    x = bulk_mean.reindex(common).to_numpy()
    y = sc_mean.reindex(common).to_numpy()
    keep = (x > 1e-3) & (y > 5e-3)
    x, y = x[keep], y[keep]
    r, p = pearsonr(np.log10(x), np.log10(y))
    s, _ = spearmanr(x, y)

    ax.scatter(x, y, s=5, alpha=0.4, color='#3182bd', edgecolors='none')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('bulk', fontsize=8)
    ax.set_ylabel('scRNA-seq', fontsize=8)
    ax.text(0.05, 0.97, title, transform=ax.transAxes,
            fontsize=9, va='top', fontweight='bold')
    ax.text(0.05, 0.82, f'$r = {r:.2f}$\n{_fmt_p(p)}\n$n = {len(x)}$',
            transform=ax.transAxes, fontsize=7.5, va='top', bbox=_ANNO_BBOX)
    ax.tick_params(labelsize=7)
    _style_ax(ax, grid=True)
    return r, s, len(x)


def lfc_corr(syn, locus2name, study, ax, y_label = True):
    sc = pd.read_csv(os.path.join(SC_DESEQ, f'deseq2_results_{study}_vs_control.csv'),
                     index_col=0)
    sc_lfc = sc['log2FoldChange']
    sc_lfc.index = harmonize(list(sc_lfc.index), syn)
    sc_lfc = sc_lfc[~sc_lfc.index.isin(SPIKE_INS)].groupby(level=0).mean()

    bulk = pd.read_csv(os.path.join(BULK_DESEQ, f'deseq2_results_{study}.csv'), index_col=0)
    bulk_lfc = bulk['log2FoldChange']
    bulk_lfc.index = harmonize_bulk(list(bulk_lfc.index), syn, locus2name)
    bulk_lfc = bulk_lfc.groupby(level=0).mean()

    common = sorted(set(sc_lfc.index) & set(bulk_lfc.index))
    x = bulk_lfc.reindex(common).to_numpy()
    y = sc_lfc.reindex(common).to_numpy()
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    r, p = pearsonr(x, y)
    s, _ = spearmanr(x, y)

    label = 'Dis-Arrest' if study == 'disrupted' else 'Reg-Arrest'
    ax.scatter(x, y, s=6, alpha=0.4, color='#de2d26', edgecolors='none')
    ax.axhline(0, color='grey', lw=0.6)
    ax.axvline(0, color='grey', lw=0.6)
    ax.set_xlabel('bulk log₂FC', fontsize=10)
    if y_label:
        ax.set_ylabel('scRNA-seq log₂FC', fontsize=10)
    ax.set_title(f'{label} vs control', fontsize=9, pad=2)
    ax.set_yticks([-10,-8,-6,-4,-2,0,2,4,6,8])
    ax.text(0.05, 0.94, f'$r = {r:.2f}$\n{_fmt_p(p)}\n$n = {len(x)}$',
            transform=ax.transAxes, fontsize=9, va='top', bbox=_ANNO_BBOX)
    _style_ax(ax, grid=True)
    return r, s, len(x)


def go_bar(ax):
    comp = pd.read_csv(GO_TABLE, index_col=0)
    dd = -np.log10(comp['disrupted_FDR'].to_numpy())
    dr = -np.log10(comp['regulated_FDR'].to_numpy())
    xpos = np.arange(len(comp))
    w = 0.38
    ax.bar(xpos - w / 2, dd, width=w, label='Dis-Arrest', color='#de2d26')
    ax.bar(xpos + w / 2, dr, width=w, label='Reg-Arrest', color='#9ecae1')
    ax.set_ylabel(r'$-\log_{10}(\mathrm{FDR})$', fontsize=11)
    ax.set_xticks(xpos)
    ax.set_xticklabels(comp['Term'], rotation=45, ha='right', fontsize=8)
    ax.set_title('Down-regulated GO terms — scRNA-seq pseudobulk', fontsize=10, pad=2)
    ax.legend(fontsize=10)
    _style_ax(ax, grid=False)


# ── main ───────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    syn        = get_gene_synonyms()
    locus2name = map_locus_to_name()
    bulk_df    = pd.read_excel(BULK_XLSX, index_col=0)
    bulk_exp_df = pd.read_csv(BULK2_CSV, index_col=0)

    fig = PanelFigure(figsize=(7, 7), label_offset=(-0.02, 0.04))

    # Panel A — three stacked mean-expression scatter plots (full left column)
    axs_A = fig.add_grid_panel([0.02, 0.04, 0.25, 0.90],
                               nrows=3, ncols=1, label='A', hspace=0.55)
    sc_specs = [
        ('Expira_biorep_t0A_filtered.csv', 'exp',       'Exponential'),
        ('sample_13a_filtered.csv',     'disrupted',  'Dis-Arrest'),
        ('sample_15b_filtered.csv',     'casp',       'Reg-Arrest'),
    ]
    for ax, (fname, prefix, title) in zip(axs_A[:, 0], sc_specs):
        if prefix == 'exp':
            r, s, n = mean_expression_corr(
                syn, locus2name, bulk_exp_df,
                os.path.join(DATA, fname), prefix, title, ax)
        else:
            r, s, n = mean_expression_corr(
                syn, locus2name, bulk_df,
                os.path.join(DATA, fname), prefix, title, ax)
        print(f'Mean expr {title:12s}: Pearson r={r:.2f}, Spearman ρ={s:.2f}, n={n}')

    # Panel B — LFC scatter: Dis-Arrest
    ax_B = fig.add_panel([0.4, 0.6, 0.25, 0.34], label='B')
    r, s, n = lfc_corr(syn, locus2name, 'disrupted', ax_B)
    print(f'LFC Dis-Arrest:      Pearson r={r:.2f}, Spearman ρ={s:.2f}, n={n}')

    # Panel C — LFC scatter: Reg-Arrest
    ax_C = fig.add_panel([0.72, 0.6, 0.25, 0.34], label='C')
    r, s, n = lfc_corr(syn, locus2name, 'regulated', ax_C, y_label=False)
    print(f'LFC Reg-Arrest:      Pearson r={r:.2f}, Spearman ρ={s:.2f}, n={n}')

    # Panel D — GO bar chart
    ax_D = fig.add_panel([0.38, 0.08, 0.59, 0.4], label='D')
    go_bar(ax_D)

    out = os.path.join(_HERE, 'figure_s3')
    fig.save(out + '.pdf', format='pdf', bbox_inches='tight')
    fig.save(out + '.png', dpi=200, bbox_inches='tight')
    print('wrote', out + '.svg')
