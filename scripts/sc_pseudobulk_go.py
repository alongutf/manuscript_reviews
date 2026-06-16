"""
scRNA-seq pseudobulk GO analysis (Reviewer #1, Comment 3).

Aggregates single-cell count matrices into per-replicate pseudobulk profiles,
runs DESeq2 (Dis-Arrest / Reg-Arrest vs control) and GO enrichment, reproducing
the bulk Section 2.1 comparison (weaker enrichment in Dis-Arrest) in single-cell
data. Gene-name harmonization follows scripts/bulk_correlations.ipynb.

Run from the scripts/ directory so that os.path.dirname(os.getcwd()) is the repo root.
"""
import os
import re
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.getcwd()) if os.path.basename(os.getcwd()) == 'scripts' else os.getcwd()
META = os.path.join(REPO, 'metadata')
# data_for_umap holds the FULL gene panel (~3900 shared genes); prefer it.
# adam_matrix_filtered lives only in data_for_paper, so resolve() falls back there.
DATA_UMAP = os.path.join(REPO, 'data_for_umap')
DATA_PAPER = os.path.join(REPO, 'data_for_paper')
DATA = DATA_UMAP  # default lookup dir for downstream scripts

# ---- condition -> replicate files (per user correction; see memory) ----
CONDITIONS = {
    'control':   ['EXP_biorep_t0A_filtered.csv', 'adam_matrix_filtered.csv'],
    'disrupted': ['sample_13a_filtered.csv', 'sample_15a_filtered.csv'],
    'regulated': ['sample_13b_filtered.csv', 'sample_15b_filtered.csv'],
}
SPIKE_INS = {'gfp', 'mcherry', 'tetr', 'laci', 'ampr', 'lelobekk'}


def resolve(fname):
    """Prefer the full-panel data_for_umap copy; fall back to data_for_paper (adam)."""
    p = os.path.join(DATA_UMAP, fname)
    return p if os.path.exists(p) else os.path.join(DATA_PAPER, fname)


# ---- gene-name harmonization (adapted from bulk_correlations.ipynb) ----
def get_gene_synonyms():
    gtf = pd.read_csv(os.path.join(META, 'genomic.gtf'), sep='\t', comment='#', header=None)
    syn = {}
    for i in range(len(gtf)):
        if gtf.iloc[i, 2] == 'gene':
            attr = gtf.iloc[i, 8]
            m = re.search(r'gene "([^"]+)"', attr)
            if m:
                primary = m.group(1)
                for s in re.findall(r'gene_synonym "([^"]+)"', attr):
                    syn[s.lower()] = primary.lower()
    return syn


def harmonize(genes, syn):
    """Return canonical lowercase gene names for a list of sc column names."""
    out = []
    for g in genes:
        g = g.replace('LELOBEKK_', '').replace('LELOBEKK', '')
        g = g.lower()
        g = syn.get(g, g)
        out.append(g)
    return out


def pseudobulk_vector(path, syn):
    """Sum counts across cells -> Series indexed by harmonized gene name."""
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)   # drop stray non-numeric cols
    totals = df.sum(axis=0)                       # gene -> summed count
    names = harmonize(list(totals.index), syn)
    s = pd.Series(totals.values, index=names)
    s = s[~s.index.isin(SPIKE_INS)]
    s = s[s.index != '']
    s = s.groupby(level=0).sum()                  # collapse duplicate names
    return s


def build_pseudobulk():
    syn = get_gene_synonyms()
    cols, meta_rows = {}, []
    for cond, files in CONDITIONS.items():
        for f in files:
            name = f.replace('_filtered.csv', '').replace('.csv', '')
            cols[name] = pseudobulk_vector(resolve(f), syn)
            meta_rows.append((name, cond))
    # intersection of genes across all replicates
    common = set.intersection(*(set(s.index) for s in cols.values()))
    common = sorted(common)
    count = pd.DataFrame({name: s.reindex(common) for name, s in cols.items()})
    count = count.round().astype(int)             # DESeq2 needs integer counts
    meta = pd.DataFrame(meta_rows, columns=['sample', 'condition']).set_index('sample')
    meta = meta.loc[count.columns]                # align order
    return count, meta


def run_deseq_contrasts(count, meta, out_dir):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    counts_T = count.T                            # samples x genes
    dds = DeseqDataSet(counts=counts_T, metadata=meta, design_factors='condition', quiet=True)
    dds.deseq2()
    results = {}
    for study in ('disrupted', 'regulated'):
        ds = DeseqStats(dds, contrast=['condition', study, 'control'], quiet=True)
        ds.summary()
        res = ds.results_df
        res.to_csv(os.path.join(out_dir, f'deseq2_results_{study}_vs_control.csv'))
        results[study] = res
    return results


def run_go(deg_df, fold, out_csv, n_study=None, p_cutoff=0.05, lfc_cutoff=1.0):
    """GO enrichment on DE genes (goatools), mirroring src.bulk_functions.run_go_enrichment.

    n_study: if given, the study set is the n_study genes with the most extreme LFC
    (lowest for fold='down'), matched between conditions; otherwise padj/LFC cutoffs.
    """
    from goatools import obo_parser
    from goatools.associations import read_gaf
    from goatools.go_enrichment import GOEnrichmentStudy
    import sys
    sys.path.insert(0, REPO)
    from src.bulk_functions import get_ID_conversion, remove_unidentified_genes, GTF_FILE

    go_dag = obo_parser.GODag(os.path.join(META, 'go-basic.obo'))
    geneid2gos = read_gaf(os.path.join(META, 'ecocyc.gaf'))
    gene_id_name = get_ID_conversion(GTF_FILE)

    deg_df = remove_unidentified_genes(deg_df, gene_id_name)
    if n_study is not None:
        # matched study set: the n_study genes with the most extreme LFC
        ranked = deg_df.sort_values('log2FoldChange', ascending=(fold == 'down'))
        sel = ranked.index[:n_study]
    elif fold == 'up':
        sel = deg_df.index[(deg_df['padj'] < p_cutoff) & (deg_df['log2FoldChange'] > lfc_cutoff)]
    else:
        sel = deg_df.index[(deg_df['padj'] < p_cutoff) & (deg_df['log2FoldChange'] < -lfc_cutoff)]

    to_id = lambda gs: [gene_id_name[g.lower()] for g in gs if g.lower() in gene_id_name]
    study_ids = to_id(set(sel))
    bg_ids = to_id(set(deg_df.index))
    print(f'  {fold}: {len(sel)} DE genes ({len(study_ids)} mapped), background {len(bg_ids)}')

    goea = GOEnrichmentStudy(bg_ids, geneid2gos, go_dag, propagate_counts=False,
                             alpha=0.05, methods=['fdr_bh'])
    res = goea.run_study(study_ids)
    res = [r for r in res if r.enrichment == 'e']     # keep enriched direction only
    full = pd.DataFrame({
        'GO_ID': [r.GO for r in res], 'Term': [r.name for r in res],
        'Category': [r.NS for r in res], 'p-value': [r.p_uncorrected for r in res],
        'FDR': [r.p_fdr_bh for r in res], 'Fold Enrichment': [r.enrichment for r in res],
        'Ratio in Study': [r.ratio_in_study for r in res],
        'Ratio in Population': [r.ratio_in_pop for r in res],
    })
    full.to_csv(out_csv, index=False)
    n_sig = (full['FDR'] < 0.05).sum()
    print(f'  -> {len(full)} terms tested, {n_sig} significant (FDR<0.05); written to {os.path.basename(out_csv)}')
    return full


def _both_sig_table(go_dis, go_reg):
    """Down-regulated GO terms significant (FDR<0.05) in BOTH conditions."""
    d = go_dis.set_index('GO_ID'); r = go_reg.set_index('GO_ID')
    both = sorted(set(d.index[d['FDR'] < 0.05]) & set(r.index[r['FDR'] < 0.05]))
    return pd.DataFrame({
        'Term': [d.loc[g, 'Term'] for g in both],
        'disrupted_FDR': [d.loc[g, 'FDR'] for g in both],
        'regulated_FDR': [r.loc[g, 'FDR'] for g in both],
    }, index=both).sort_values('regulated_FDR')


def make_comparison_figure(go, sizes, out_svg):
    """Down-regulated GO comparison for matched study sets of size `sizes`.

    `go` is {(study, fold, N): full_go_table}. For each N, one panel showing the
    -log10(FDR) of GO terms significant (FDR<0.05) in BOTH conditions, Dis vs Reg.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import wilcoxon

    plt.style.use('ggplot')
    fig, axes = plt.subplots(1, len(sizes), figsize=(7.5 * len(sizes), 6.5), squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(top=0.90, bottom=0.40, left=0.06, right=0.98, wspace=0.22)
    summary = []
    for ax, N in zip(axes, sizes):
        comp = _both_sig_table(go[('disrupted', 'down', N)], go[('regulated', 'down', N)])
        dd = -np.log10(comp['disrupted_FDR'].to_numpy())
        dr = -np.log10(comp['regulated_FDR'].to_numpy())
        x = np.arange(len(comp)); w = 0.4
        ax.bar(x - w/2, dd, width=w, label='Dis-Arrest', color='#de2d26')
        ax.bar(x + w/2, dr, width=w, label='Reg-Arrest', color='#9ecae1')
        ax.set_ylabel(r'$-\log_{10}(FDR)$', fontsize=13)
        ax.set_xticks(x); ax.set_xticklabels(comp['Term'], rotation=45, ha='right', fontsize=8)
        try:
            pw = wilcoxon(dr, dd, alternative='greater')[1]
        except ValueError:
            pw = np.nan
        ax.set_title(f'Down-regulated, study set = {N} genes\n'
                     f'{len(comp)} terms sig. in both; Wilcoxon p={pw:.3f}', fontsize=12)
        ax.legend(fontsize=11)
        comp.to_csv(out_svg.replace('.svg', f'_down_N{N}_table.csv'))
        summary.append((N, len(comp), np.median(dd), np.median(dr), pw))

    fig.savefig(out_svg, format='svg'); fig.savefig(out_svg.replace('.svg', '.png'), dpi=200)
    print('  N    | both-sig terms | median -log10FDR Dis | Reg | Wilcoxon(Reg>Dis) p')
    for N, n, md, mr, pw in summary:
        print(f'  {N:<4} | {n:^14} | {md:>19.2f} | {mr:>3.2f} | {pw:.3f}')
    return summary


if __name__ == '__main__':
    count, meta = build_pseudobulk()
    print('pseudobulk matrix:', count.shape, '(genes x samples)')
    print(meta)
    print('library sizes:\n', count.sum(axis=0))
    deseq_out = os.path.join(REPO, 'results', 'deseq_results', 'sc_pseudobulk')
    go_out = os.path.join(REPO, 'results', 'GO_results', 'sc_pseudobulk')
    os.makedirs(deseq_out, exist_ok=True)
    os.makedirs(go_out, exist_ok=True)
    count.to_csv(os.path.join(deseq_out, 'sc_pseudobulk_counts.csv'))
    meta.to_csv(os.path.join(deseq_out, 'sc_pseudobulk_metadata.csv'))

    print('\n=== DESeq2 ===')
    res = run_deseq_contrasts(count, meta, deseq_out)
    for study, r in res.items():
        sig = r[(r['padj'] < 0.05) & (r['log2FoldChange'].abs() > 1)]
        print(f'{study}: {len(sig)} DE genes (padj<0.05, |LFC|>1); '
              f'down={sum(sig["log2FoldChange"]<0)}, up={sum(sig["log2FoldChange"]>0)}')

    print('\n=== GO enrichment (down-regulated, matched study set) ===')
    sizes = (500,)
    go = {}
    for study in ('disrupted', 'regulated'):
        for N in sizes:
            print(f'{study} / down / N={N}')
            go[(study, 'down', N)] = run_go(
                res[study], 'down', os.path.join(go_out, f'GO_enrichment_{study}_down_N{N}.csv'),
                n_study=N)

    print('\n=== Section 2.1 comparison (terms significant in BOTH conditions) ===')
    make_comparison_figure(go, sizes, os.path.join(go_out, 'sc_GO_FDR_comparison.svg'))
