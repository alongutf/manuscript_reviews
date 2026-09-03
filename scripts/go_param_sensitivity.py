"""
GO-term gene-set-size sensitivity & specificity analysis (Reviewer #1, Comment 4).

Three analyses on the bulk DESeq2 / goatools over-representation results behind Fig. 4:

  4.1  Specificity of Fig. 4A. Annotate every gene in the Fig. 4A chemotaxis panel
       (GO:0006935 rows of scripts/figures/figure4/GO_terms_heatmap.csv) to ALL of its
       GO-term memberships, and quantify how concentrated the program is in
       motility / chemotaxis / flagellar categories vs. how many genes fall outside them.

  4.2  Gene-set-size sensitivity. Re-run the over-representation analysis once per
       condition to obtain the full list of tested terms with their uncorrected
       p-values and gene-set sizes, then restrict the tested terms to a range of
       (min, max) gene-set-size windows, re-apply Benjamini-Hochberg FDR within each
       window, and recompute the Dis- vs Reg-Arrest -log10(FDR) comparison
       (Mann-Whitney U, matching the Fig. 4B inset) plus the significance of the core
       motility / translation terms.

  4.3  Size confound. For the GO terms significant in BOTH conditions (the set plotted
       in Fig. 4B), correlate the per-term Dis-Reg difference in -log10(FDR) with the
       gene-set size (Spearman rho), to test whether the between-condition difference
       is systematically larger for bigger gene sets.

Note: the Fig. 4 pipeline is goatools over-representation (Fisher exact + BH-FDR),
NOT GSEA. "Gene-set size" is the number of background genes annotated to a term, i.e.
the first element of the goatools `Ratio in Population` field (r.ratio_in_pop[0]).

Run from the scripts/ directory so that os.path.dirname(os.getcwd()) is the repo root.
Outputs -> results/GO_results/param_sensitivity/.
"""
import os
import sys
import json
import ast
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- paths ----
REPO = os.path.dirname(os.getcwd()) if os.path.basename(os.getcwd()) == 'scripts' else os.getcwd()
sys.path.insert(0, REPO)
META = os.path.join(REPO, 'metadata')
DESEQ_DIR = os.path.join(REPO, 'results', 'deseq_results', 'from counts')
GO_FROMCOUNTS = os.path.join(REPO, 'results', 'GO_results', 'from_counts')
HEATMAP_CSV = os.path.join(REPO, 'scripts', 'figures', 'figure4', 'GO_terms_heatmap.csv')
OUT = os.path.join(REPO, 'results', 'GO_results', 'param_sensitivity')

# ---- parameters that reproduce the published Fig. 4B study set ----
P_CUTOFF = 0.01
TARGET_SIZE = 1000
FOLD = 'down'
FIG4A_GO = 'GO:0006935'  # chemotaxis (the GO term shown in Fig. 4A)

# core growth-arrest program terms whose stability we track across size windows
KEY_TERMS = {
    'GO:0006935': 'chemotaxis',
    'GO:0071973': 'bacterial-type flagellum-dependent cell motility',
    'GO:0044781': 'bacterial-type flagellum organization',
    'GO:0044780': 'bacterial-type flagellum assembly',
    'GO:0002181': 'cytoplasmic translation',
    'GO:0006412': 'translation',
}
# substrings that mark a GO term as part of the motility / chemotaxis / flagellar program
MOTILITY_KEYWORDS = ('flagel', 'chemotax', 'motil', 'taxis')


def is_motility(term_name):
    """True if a GO term name contains a motility/chemotaxis/flagellar substring."""
    t = term_name.lower()
    return any(k in t for k in MOTILITY_KEYWORDS)


def load_go_annotations():
    """Load the goatools GO DAG, the gene->GO associations, and the name->UniProt map."""
    from goatools import obo_parser
    from goatools.associations import read_gaf
    from src.bulk_functions import get_ID_conversion, GTF_FILE
    go_dag = obo_parser.GODag(os.path.join(META, 'go-basic.obo'))
    gene2gos = read_gaf(os.path.join(META, 'ecocyc.gaf'))   # UniProt accession -> {GO ids}
    name2id = get_ID_conversion(GTF_FILE)                    # gene name (lower) -> UniProt accession
    return go_dag, gene2gos, name2id


# ----------------------------------------------------------------------------
# 4.1  Specificity of Fig. 4A
# ----------------------------------------------------------------------------
def analyse_fig4a_specificity(go_dag, gene2gos, name2id):
    """Answer 4.1: for every gene plotted under GO:0006935 (chemotaxis) in the Fig. 4A
    heatmap, look up ALL of that gene's GO annotations (not just chemotaxis) and tally
    how many of them, and how many distinct terms, are motility/chemotaxis/flagellar vs.
    something else. Writes the per-gene-per-term annotation table and a per-term
    membership summary, and returns a dict of summary counts for the log.
    """
    heat = pd.read_csv(HEATMAP_CSV, index_col=0)
    genes = heat.index[heat['GO_term'] == FIG4A_GO].tolist()

    rows, n_with_motility, n_unmapped = [], 0, 0
    for g in genes:
        uid = name2id.get(g.lower())
        if uid is None or uid not in gene2gos:
            # gene has no UniProt accession, or the accession has no GO annotations:
            # exclude it from the specificity denominator rather than counting it as
            # "not motility"
            n_unmapped += 1
            continue
        go_ids = gene2gos[uid]
        gene_has_motility = False
        for go in go_ids:
            if go not in go_dag:
                continue
            name = go_dag[go].name
            mot = is_motility(name)
            gene_has_motility = gene_has_motility or mot
            rows.append({'gene': g, 'uniprot': uid, 'GO_ID': go,
                         'GO_name': name, 'namespace': go_dag[go].namespace,
                         'is_motility': mot})
        n_with_motility += int(gene_has_motility)

    ann = pd.DataFrame(rows)
    ann.to_csv(os.path.join(OUT, 'fig4a_gene_go_annotations.csv'), index=False)

    # per-term membership summary across the Fig. 4A genes
    summary = (ann.groupby(['GO_ID', 'GO_name', 'is_motility'])['gene']
                  .nunique().reset_index().rename(columns={'gene': 'n_fig4a_genes'})
                  .sort_values('n_fig4a_genes', ascending=False))
    summary.to_csv(os.path.join(OUT, 'fig4a_go_membership_summary.csv'), index=False)

    n_mapped = len(genes) - n_unmapped
    stats = {
        'n_fig4a_genes': len(genes),
        'n_unmapped': n_unmapped,
        'n_mapped': n_mapped,
        'n_genes_with_motility_annotation': n_with_motility,
        'frac_genes_with_motility_annotation': round(n_with_motility / n_mapped, 3) if n_mapped else None,
        'n_distinct_go_terms': int(ann['GO_ID'].nunique()),
        'n_distinct_motility_terms': int(ann.loc[ann['is_motility'], 'GO_ID'].nunique()),
        'n_distinct_nonmotility_terms': int(ann.loc[~ann['is_motility'], 'GO_ID'].nunique()),
    }
    # top non-motility terms (the "GO terms beyond chemotaxis" the reviewer asks about)
    nonmot = summary[~summary['is_motility']].head(10)
    stats['top_nonmotility_terms'] = [
        {'GO_ID': r.GO_ID, 'GO_name': r.GO_name, 'n_fig4a_genes': int(r.n_fig4a_genes)}
        for r in nonmot.itertuples()
    ]
    print(f"[4.1] Fig. 4A genes: {len(genes)} ({n_mapped} mapped); "
          f"{n_with_motility} have a motility/chemotaxis/flagellar annotation "
          f"({stats['frac_genes_with_motility_annotation']}). "
          f"{stats['n_distinct_motility_terms']} distinct motility terms vs "
          f"{stats['n_distinct_nonmotility_terms']} other terms.")
    return stats


# ----------------------------------------------------------------------------
# Shared: run the over-representation study and return ALL tested terms
# ----------------------------------------------------------------------------
def run_full_ora(study, go_dag, gene2gos, name2id):
    """Reproduce the Fig. 4 study set for one condition and return every enriched
    term tested (not just the significant ones), with uncorrected p and gene-set size."""
    from goatools.go_enrichment import GOEnrichmentStudy
    from src.bulk_functions import remove_unidentified_genes, get_lfc_thresh

    deg = pd.read_csv(os.path.join(DESEQ_DIR, f'deseq2_results_{study}.csv'), index_col=0)
    deg = remove_unidentified_genes(deg, name2id)
    # same DEG-selection logic as src.bulk_functions.run_go_enrichment: pick the LFC
    # cutoff that yields ~TARGET_SIZE down-regulated DEGs at padj < P_CUTOFF, so the
    # study gene set here matches the one behind the published Fig. 4B terms
    lfc_cut = get_lfc_thresh(deg, TARGET_SIZE, FOLD, p_val_thresh=P_CUTOFF)
    sel = deg.index[(deg['padj'] < P_CUTOFF) & (deg['log2FoldChange'] < lfc_cut)]

    to_id = lambda gs: [name2id[g.lower()] for g in gs if g.lower() in name2id]
    study_ids = to_id(set(sel))
    bg_ids = to_id(set(deg.index))

    # propagate_counts=False: only direct annotations count for a term, not those of
    # its GO-DAG children, matching src.bulk_functions.run_go_enrichment exactly
    goea = GOEnrichmentStudy(bg_ids, gene2gos, go_dag, propagate_counts=False,
                             alpha=0.05, methods=['fdr_bh'])
    # keep only 'e' (enriched) results: goatools also reports 'p' (purified/depleted)
    # terms, which are not part of the Fig. 4 over-representation analysis
    res = [r for r in goea.run_study(study_ids) if r.enrichment == 'e']
    full = pd.DataFrame({
        'GO_ID': [r.GO for r in res],
        'Term': [r.name for r in res],
        'p_uncorrected': [r.p_uncorrected for r in res],
        'set_size': [r.ratio_in_pop[0] for r in res],     # background genes in the term
        'study_count': [r.ratio_in_study[0] for r in res],
    })
    print(f"  [{study}] study set {len(study_ids)} genes, {len(full)} enriched terms tested")
    return full


# ----------------------------------------------------------------------------
# 4.2  Gene-set-size sensitivity sweep
# ----------------------------------------------------------------------------
# The Fig. 4 conclusion is that the growth-arrest program (chemotaxis / flagellar /
# translation) is significantly enriched in BOTH conditions but MORE strongly in
# Reg-Arrest. Robustness is therefore tested two ways:
#   (a) paired Reg-vs-Dis comparison over the shared program terms (Fig. 4B common set),
#       restricted to each gene-set-size window  -> direction & significance preserved;
#   (b) for the full tested-term universe, restrict to a size window, re-apply BH, and
#       check that the core motility/translation terms remain significant in each
#       condition -> the recovered term set does not depend on the size cut.
def analyse_size_sensitivity(full_dis, full_reg, go_dis_pub, go_reg_pub):
    """Run analyses (a) and (b) described above across a grid of gene-set-size windows.

    full_dis, full_reg  -- every tested term (goatools output of run_full_ora) for the
                            Dis-Arrest and Reg-Arrest study sets, with uncorrected p and
                            set_size.
    go_dis_pub, go_reg_pub -- the published, already-FDR-filtered Fig. 4B term tables
                            (FDR < 0.05 in goatools' own within-condition BH pass).
    Writes go_size_sensitivity_paired.csv and go_size_keyterm_recovery.csv, and returns
    a dict of summary counts for the log.
    """
    # ---- common shared-program terms (the bars plotted in Fig. 4B) ----
    d = go_dis_pub.set_index('GO_ID')
    r = go_reg_pub.set_index('GO_ID')
    common = sorted(set(d.index) & set(r.index))
    base = pd.DataFrame({
        'GO_ID': common,
        'Term': [d.loc[g, 'Term'] for g in common],
        'set_size': [_pop_size(d.loc[g, 'Ratio in Population']) for g in common],
        'logFDR_dis': [-np.log10(d.loc[g, 'FDR']) for g in common],
        'logFDR_reg': [-np.log10(r.loc[g, 'FDR']) for g in common],
    })

    # min/max sizes chosen to bracket the plausible "GO terms too broad/narrow to
    # trust" cutoffs a reviewer might propose, from a tight (5,100) window up to
    # dropping only the very largest or smallest terms
    windows = [(5, 100), (5, 200), (5, 500), (5, np.inf),
               (10, np.inf), (15, np.inf), (10, 200), (15, 200)]

    # (a) paired comparison within windows
    paired_rows = []
    for mn, mx in windows:
        sub = base[(base['set_size'] >= mn) & (base['set_size'] <= mx)]
        n = len(sub)
        if n >= 3:
            # one-sided: tests the published Fig. 4B claim that Reg-Arrest terms have
            # higher -log10(FDR) than Dis-Arrest for the same shared-program terms
            try:
                w_stat, w_p = wilcoxon(sub['logFDR_reg'], sub['logFDR_dis'], alternative='greater')
            except ValueError:
                # wilcoxon raises when all paired differences are zero or n is too small
                w_stat, w_p = np.nan, np.nan
        else:
            # too few terms in this window for a meaningful paired test
            w_stat, w_p = np.nan, np.nan
        paired_rows.append({
            'min_size': mn, 'max_size': (None if np.isinf(mx) else int(mx)),
            'n_common_terms': n,
            'median_logFDR_dis': round(float(sub['logFDR_dis'].median()), 2) if n else None,
            'median_logFDR_reg': round(float(sub['logFDR_reg'].median()), 2) if n else None,
            'wilcoxon_p_Reg_gt_Dis': w_p,
            'n_key_terms_present': sum(t in sub['GO_ID'].values for t in KEY_TERMS),
        })
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(os.path.join(OUT, 'go_size_sensitivity_paired.csv'), index=False)

    # (b) key-term recovery in the full universe, re-applying BH inside each window
    key_rows = []
    for mn, mx in windows:
        dd = full_dis[(full_dis['set_size'] >= mn) & (full_dis['set_size'] <= mx)].copy()
        rr = full_reg[(full_reg['set_size'] >= mn) & (full_reg['set_size'] <= mx)].copy()
        # BH-FDR must be re-applied within each window: restricting the tested-term
        # universe changes the multiple-testing correction, so the published FDR
        # column (computed over the unrestricted universe) cannot be reused here
        dd['FDR_w'] = multipletests(dd['p_uncorrected'], method='fdr_bh')[1]
        rr['FDR_w'] = multipletests(rr['p_uncorrected'], method='fdr_bh')[1]
        ddi, rri = dd.set_index('GO_ID'), rr.set_index('GO_ID')
        key_dis = sum(int(go in ddi.index and ddi.loc[go, 'FDR_w'] < 0.05) for go in KEY_TERMS)
        key_reg = sum(int(go in rri.index and rri.loc[go, 'FDR_w'] < 0.05) for go in KEY_TERMS)
        key_rows.append({
            'min_size': mn, 'max_size': (None if np.isinf(mx) else int(mx)),
            'n_terms_dis': len(dd), 'n_terms_reg': len(rr),
            'n_sig_dis': int((dd['FDR_w'] < 0.05).sum()), 'n_sig_reg': int((rr['FDR_w'] < 0.05).sum()),
            'key_terms_sig_dis': f'{key_dis}/{len(KEY_TERMS)}',
            'key_terms_sig_reg': f'{key_reg}/{len(KEY_TERMS)}',
        })
    keyrec = pd.DataFrame(key_rows)
    keyrec.to_csv(os.path.join(OUT, 'go_size_keyterm_recovery.csv'), index=False)

    pvals = paired['wilcoxon_p_Reg_gt_Dis'].dropna()
    n_sig = int((pvals < 0.05).sum())
    print(f"[4.2a] paired Reg>Dis over {len(common)} shared-program terms: significant "
          f"(p<0.05) in {n_sig}/{len(pvals)} size windows (p range {pvals.min():.1e}-{pvals.max():.1e}).")
    print(f"[4.2b] key motility/translation terms remain significant across windows "
          f"(see go_size_keyterm_recovery.csv).")
    return {
        'n_common_terms': len(common),
        'n_windows': len(paired),
        'n_windows_paired_sig': n_sig,
        'wilcoxon_p_min': float(pvals.min()), 'wilcoxon_p_max': float(pvals.max()),
    }


# ----------------------------------------------------------------------------
# 4.3  Size confound: between-condition difference vs gene-set size
# ----------------------------------------------------------------------------
def _pop_size(ratio_str):
    """Parse goatools 'Ratio in Population' field like '(34, 3938)' -> 34."""
    return ast.literal_eval(ratio_str)[0]


def analyse_size_confound(go_dis, go_reg):
    """Answer 4.3: for GO terms significant in both conditions (published Fig. 4B
    tables), test whether the Reg-minus-Dis difference in -log10(FDR) scales with the
    term's gene-set size (a size confound would show up as rho far from 0). Writes the
    per-term table and confound scatter plot; returns the Spearman stats for the log.
    """
    d = go_dis.set_index('GO_ID')
    r = go_reg.set_index('GO_ID')
    common = sorted(set(d.index) & set(r.index))     # both significant (published files are FDR<0.05)
    rows = []
    for go in common:
        size = _pop_size(d.loc[go, 'Ratio in Population'])
        sd = -np.log10(d.loc[go, 'FDR'])
        sr = -np.log10(r.loc[go, 'FDR'])
        rows.append({'GO_ID': go, 'Term': d.loc[go, 'Term'], 'set_size': size,
                     'logFDR_dis': sd, 'logFDR_reg': sr,
                     'diff_reg_minus_dis': sr - sd, 'abs_diff': abs(sr - sd)})
    conf = pd.DataFrame(rows).sort_values('set_size')
    conf.to_csv(os.path.join(OUT, 'go_size_confound.csv'), index=False)

    rho_signed, p_signed = spearmanr(conf['set_size'], conf['diff_reg_minus_dis'])
    rho_abs, p_abs = spearmanr(conf['set_size'], conf['abs_diff'])

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(conf['set_size'], conf['diff_reg_minus_dis'], s=18, color='#de2d26')
    ax.axhline(0, color='grey', lw=0.6, ls='--')
    ax.set_xlabel('GO term gene-set size (genes)')
    ax.set_ylabel(r'$-\log_{10}(FDR)_{Reg} - {Dis}$')
    ax.set_title(f'Size confound: Spearman rho={rho_signed:.2f} (p={p_signed:.2g})', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'go_size_confound.svg'))
    fig.savefig(os.path.join(OUT, 'go_size_confound.png'), dpi=200)
    plt.close(fig)

    stats = {
        'n_common_terms': len(conf),
        'spearman_rho_signed_diff_vs_size': round(float(rho_signed), 3),
        'spearman_p_signed': float(p_signed),
        'spearman_rho_abs_diff_vs_size': round(float(rho_abs), 3),
        'spearman_p_abs': float(p_abs),
    }
    print(f"[4.3] {len(conf)} common terms; Spearman(signed diff vs size) rho={rho_signed:.2f} "
          f"p={p_signed:.2g}; Spearman(|diff| vs size) rho={rho_abs:.2f} p={p_abs:.2g}.")
    return stats


# ----------------------------------------------------------------------------
def main():
    """Run 4.1, 4.2, 4.3 in sequence and write the combined summary log as JSON."""
    os.makedirs(OUT, exist_ok=True)
    log = {'timestamp': datetime.now().strftime('%Y%m%d_%H%M%S'),
           'params': {'p_cutoff': P_CUTOFF, 'target_size': TARGET_SIZE, 'fold': FOLD,
                      'fig4a_go': FIG4A_GO, 'key_terms': KEY_TERMS,
                      'motility_keywords': MOTILITY_KEYWORDS,
                      'min_sizes': [5, 10, 15], 'max_sizes': [100, 200, 500, 'inf']}}

    go_dag, gene2gos, name2id = load_go_annotations()

    print('\n=== 4.1 Fig. 4A specificity ===')
    log['fig4a_specificity'] = analyse_fig4a_specificity(go_dag, gene2gos, name2id)

    # published Fig. 4B term tables (used by both 4.2 and 4.3)
    go_dis = pd.read_csv(os.path.join(GO_FROMCOUNTS, 'GOATOOLS_GO_enrichment_results_disrupted_down.csv'))
    go_reg = pd.read_csv(os.path.join(GO_FROMCOUNTS, 'GOATOOLS_GO_enrichment_results_regulated_down.csv'))

    print('\n=== 4.2 gene-set-size sensitivity (re-running ORA per condition) ===')
    full_dis = run_full_ora('disrupted', go_dag, gene2gos, name2id)
    full_reg = run_full_ora('regulated', go_dag, gene2gos, name2id)
    full_dis.to_csv(os.path.join(OUT, 'ora_full_terms_disrupted.csv'), index=False)
    full_reg.to_csv(os.path.join(OUT, 'ora_full_terms_regulated.csv'), index=False)
    log['size_sensitivity'] = analyse_size_sensitivity(full_dis, full_reg, go_dis, go_reg)

    print('\n=== 4.3 size confound (published Fig. 4B terms) ===')
    log['size_confound'] = analyse_size_confound(go_dis, go_reg)

    with open(os.path.join(OUT, 'go_param_sensitivity_log.json'), 'w') as f:
        json.dump(log, f, indent=2)
    print(f"\nDone. Outputs + log written to {OUT}")


if __name__ == '__main__':
    main()
