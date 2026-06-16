"""
scRNA-seq vs bulk RNA-seq correlation (Reviewer #1, Comment 3, part 3).

Two panels:
  (A) mean-expression concordance: per-gene mean expression, Dis-Arrest sc vs
      Disrupted bulk (reproduces scripts/bulk_correlations.ipynb).
  (B) fold-change concordance: sc pseudobulk DESeq2 log2FC vs bulk DESeq2 log2FC
      (disrupted vs control), the data that actually feed GO enrichment.

Run from scripts/ so os.path.dirname(os.getcwd()) is the repo root.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

from sc_pseudobulk_go import REPO, META, DATA, get_gene_synonyms, harmonize, SPIKE_INS

BULK_XLSX = os.path.join(REPO, 'bulk_data', 'All bulk data.xlsx')
BULK_DESEQ = os.path.join(REPO, 'results', 'deseq_results', 'ercc_norm')
SC_DESEQ = os.path.join(REPO, 'results', 'deseq_results', 'sc_pseudobulk')
OUT = os.path.join(REPO, 'results', 'GO_results', 'sc_pseudobulk')


def map_locus_to_name():
    """locus_tag -> gene Name, from gffFile_combined.gff (for bulk locus-tag indices)."""
    gff = pd.read_csv(os.path.join(META, 'gffFile_combined.gff'), sep='\t', comment='#', header=None)
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
    """bulk index -> canonical lowercase names (resolve locus tags first)."""
    mapped = [locus2name.get(g, g) for g in genes]
    return harmonize(mapped, syn)


def mean_expression_corr(syn, locus2name, ax):
    # sc Dis-Arrest mean expression (sample_15a)
    sc = pd.read_csv(os.path.join(DATA, 'sample_15a_filtered.csv'), index_col=0)
    sc_mean = sc.mean(axis=0)
    sc_mean.index = harmonize(list(sc_mean.index), syn)
    sc_mean = sc_mean[~sc_mean.index.isin(SPIKE_INS)].groupby(level=0).sum()
    # bulk Disrupted mean expression
    bulk = pd.read_excel(BULK_XLSX, index_col=0)
    dis_cols = [c for c in bulk.columns if str(c).lower().startswith('disrupted')]
    bulk_mean = bulk[dis_cols].mean(axis=1)
    bulk_mean.index = harmonize_bulk(list(bulk_mean.index), syn, locus2name)
    bulk_mean = bulk_mean.groupby(level=0).sum()
    # normalize each to fractions (counts differ in scale) and intersect
    sc_frac = sc_mean / sc_mean.sum()
    bulk_frac = bulk_mean / bulk_mean.sum()
    common = sorted(set(sc_frac.index) & set(bulk_frac.index))
    x = bulk_frac.reindex(common).to_numpy()
    y = sc_frac.reindex(common).to_numpy()
    keep = (x > 1e-6) & (y > 1e-6)
    x, y = x[keep], y[keep]
    r, _ = pearsonr(np.log10(x), np.log10(y))
    s, _ = spearmanr(x, y)
    ax.scatter(x, y, s=8, alpha=0.4, color='#3182bd', edgecolors='none')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('bulk mean expression (fraction)', fontsize=12)
    ax.set_ylabel('scRNA-seq mean expression (fraction)', fontsize=12)
    ax.set_title('Mean expression (Dis-Arrest)', fontsize=12)
    ax.text(0.05, 0.92, f'Pearson r = {r:.2f}\nSpearman ρ = {s:.2f}\nn = {len(x)} genes',
            transform=ax.transAxes, fontsize=11, va='top')
    return r, s, len(x)


def lfc_corr(syn, locus2name, study, ax):
    sc = pd.read_csv(os.path.join(SC_DESEQ, f'deseq2_results_{study}_vs_control.csv'), index_col=0)
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
    r, _ = pearsonr(x, y)
    s, _ = spearmanr(x, y)
    label = 'Dis-Arrest' if study == 'disrupted' else 'Reg-Arrest'
    ax.scatter(x, y, s=8, alpha=0.4, color='#de2d26', edgecolors='none')
    ax.axhline(0, color='grey', lw=0.6); ax.axvline(0, color='grey', lw=0.6)
    ax.set_xlabel('bulk log2 fold-change', fontsize=12)
    ax.set_ylabel('scRNA-seq log2 fold-change', fontsize=12)
    ax.set_title(f'Fold-change ({label} vs control)', fontsize=12)
    ax.text(0.05, 0.92, f'Pearson r = {r:.2f}\nSpearman ρ = {s:.2f}\nn = {len(x)} genes',
            transform=ax.transAxes, fontsize=11, va='top')
    return r, s, len(x)


if __name__ == '__main__':
    syn = get_gene_synonyms()
    locus2name = map_locus_to_name()
    plt.style.use('ggplot')
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.13, top=0.9, wspace=0.3)
    print('Mean expression:', mean_expression_corr(syn, locus2name, axes[0]))
    print('LFC disrupted: ', lfc_corr(syn, locus2name, 'disrupted', axes[1]))
    print('LFC regulated: ', lfc_corr(syn, locus2name, 'regulated', axes[2]))
    svg = os.path.join(OUT, 'sc_bulk_correlation.svg')
    fig.savefig(svg, format='svg'); fig.savefig(svg.replace('.svg', '.png'), dpi=200)
    print('wrote', svg)
