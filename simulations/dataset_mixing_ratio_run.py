"""
Experimental dataset-mixing series with repeats -- Reviewer #1, comments 1 and 2.3.

The question
------------
Does mixing two experimentally distinct cell populations RAISE GMP-Cor, as the
simulated inverted-sub-population scenario predicts? If so, a low observed GMP-Cor
cannot be an artefact of population heterogeneity.

Relation to `dataset_mixing_gmpcor_run.py`
------------------------------------------
That script computes GMP-Cor for a single 50/50 mixture and reports one number
(58.88 for Exp + VapC-2h). It has two gaps this runner closes:

  1. **No reference points.** The mixture was never compared with the two pure
     populations computed in the same frame. Since GMP-Cor is extensive in the gene
     count and its noise threshold is a pure function of matrix shape, an absolute
     value on its own carries no argument -- and the 58.88 mixture had 1000 cells
     against 500-cell endpoints, which is not a valid comparison.
  2. **No repeats.** A single draw, with no estimate of how much the value moves
     when different cells are drawn or the scramble is redrawn.

Here every point -- both pure populations and every mixture -- is evaluated at the
SAME total cell number and on the SAME fixed gene panel, with `--repeats` independent
cell draws each.

Method
------
  gene space : case-folded intersection of the two matrices, reporter genes removed
  cell pools : the `--n-pool` highest-total-count cells from each dataset
  gene panel : top `--n-genes` by Fano factor (var/mean) on the combined pool,
               selected ONCE and held fixed at every ratio
  cells      : `--n-total` per point, drawn uniformly without replacement --
               `ratio * n_total` from dataset 1 and the rest from dataset 2
  GMP-Cor    : sum( max(lambda_i - lambda*_scrambled, 0) ) via
               analysis_functions.get_eig_dist(norm=True, norm_method='sum',
               norm_sum=50) -- identical settings to every other GMP-Cor here

For each mixture the group-centered GMP-Cor is also computed (each dataset's cells
centered on their own gene means, which removes all between-population structure),
giving dGMP = 1 - centered/raw: the fraction of GMP-Cor that IS the mixture rather
than within-cell coordination. Separation along the group axis is reported alongside.

Outputs (results/simulation_results/):
  logs/<stem>_<timestamp>.json     full parameters and per-repeat records
  raw/<stem>_<timestamp>.txt       human-readable summary
  raw/<stem>_<timestamp>.csv       one row per (ratio, repeat)
  figures/<stem>_<timestamp>.svg / .png
"""

import os
import sys
import json
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.simulations import gmp_cor, group_centered, separation_metrics  # noqa: E402

# Reporter / plasmid-marker / rRNA genes excluded from every panel (case-insensitive).
# The rRNA and LELOBEKK entries matter only for the data_for_paper matrices; the
# data_for_umap ones already carry them or not depending on the pipeline that built them.
EXCLUDE_GENES = {'gfp', 'laci', 'ampr', '16s_mature', '16s_unprocessed',
                 'lelobekk', 'kanr', 'mcherry', 'tetr'}

# What the paper pipeline actually drops (scripts/test.ipynb get_data_for_plot, and the
# equate_dims panel construction): an EXACT-CASE list. Several data_for_paper matrices
# store gene names in lower case, so 'LELOBEKK'/'kanR'/'mCherry' do not match and those
# genes survive into the published per-sample GMP-Cor. In Expira_biorep_t0A, `laci` alone
# is 38% of all counts. `--reporter-handling published` reproduces that behaviour so a
# mixture can be compared with results/data_metrics/data_metrics.csv; `clean` removes the
# reporters properly and gives a lower, more defensible value.
PUBLISHED_EXACT_DROP = {'16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR', 'mCherry'}
STEM = 'dataset_mixing_ratio'

_SIM_RESULTS = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
_FIG_DIR = os.path.join(_SIM_RESULTS, 'figures')
_RAW_DIR = os.path.join(_SIM_RESULTS, 'raw')
_LOG_DIR = os.path.join(_SIM_RESULTS, 'logs')


def parse_args():
    """Command-line options controlling the two source datasets and the mixing design."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--file1', default='EXP_biorep_t0A_filtered.csv',
                   help='dataset 1; ratio 1.0 is this population alone')
    p.add_argument('--file2', default='VAPC_biorep_t2A_filtered.csv',
                   help='dataset 2; ratio 0.0 is this population alone')
    p.add_argument('--label1', default='Exponential')
    p.add_argument('--label2', default='VapC: 2h')
    p.add_argument('--data-dir', default='data_for_umap')
    p.add_argument('--ratios', type=float, nargs='+',
                   default=[0.0, 0.25, 0.5, 0.75, 1.0],
                   help='fraction of cells taken from dataset 1')
    p.add_argument('--n-total', type=int, default=900,
                   help='total cells at EVERY ratio. Keep below the smaller pool so '
                        'that the pure endpoints still involve a genuine random draw')
    p.add_argument('--n-pool', type=int, default=1000,
                   help='highest-total-count cells taken from each dataset')
    p.add_argument('--n-genes', type=int, default=2000,
                   help='genes kept by Fano factor, selected once on the combined pool')
    p.add_argument('--gene-space', choices=['shared', 'union'], default='shared',
                   help="build the panel from the intersection ('shared') or the union "
                        "of the two gene sets. WITH 'union', genes absent from one "
                        "dataset are zero-filled for all of that dataset's cells, which "
                        "is a perfect artificial separator between the two populations; "
                        "the run reports how many such genes the Fano cut selected")
    p.add_argument('--reporter-handling', choices=['clean', 'published'], default='clean',
                   help="'clean' drops reporter/plasmid genes case-insensitively; "
                        "'published' reproduces the paper's exact-case drop list, which "
                        "leaves lower-cased reporter genes in and is what the published "
                        "per-sample values in data_metrics.csv were computed with")
    p.add_argument('--all-cells', action='store_true',
                   help='ignore --n-total and use EVERY cell of each pool at each ratio '
                        '(pure points use one pool in full, mixtures use both in full). '
                        'Makes n differ between points, which favours the mixture -- the '
                        'summary flags this')
    p.add_argument('--repeats', type=int, default=5)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--norm-sum', type=float, default=50)
    return p.parse_args()


# ── Data ─────────────────────────────────────────────────────────────────────

def load_matrix(path, fname, reporter_handling='clean'):
    """Read one cells x genes CSV, drop bookkeeping columns, and case-fold gene names.

    `reporter_handling='published'` additionally reproduces the paper pipeline's
    exact-case reporter drop, see PUBLISHED_EXACT_DROP above.
    """
    d = pd.read_csv(path, index_col=0).fillna(0.0)
    drop = [c for c in d.columns
            if str(c).lower().startswith('unnamed') or str(c).startswith('INTR_')]
    if drop:
        d = d.drop(columns=drop)
        print(f'  {fname}: dropped {len(drop)} non-gene columns')
    if reporter_handling == 'published':
        exact = [c for c in d.columns if c in PUBLISHED_EXACT_DROP]
        d = d.drop(columns=exact)
        kept = sorted(c for c in d.columns if str(c).casefold() in EXCLUDE_GENES)
        print(f'  {fname}: exact-case drop {exact}; reporter genes RETAINED: {kept}')
    d.columns = [str(c).casefold() for c in d.columns]
    if d.columns.duplicated().any():
        d = d.T.groupby(level=0).sum().T
    return d


def build_pools(args):
    """Load both datasets, build the fixed gene panel, and return the two cell pools.

    Returns (a, b, genes, reporters_present, n_shared, n_exclusive):
      a, b               cells x len(genes) float arrays, the highest-count cells of
                         each dataset restricted to the fixed panel
      genes              the panel itself, gene names in fixed order
      reporters_present  reporter/plasmid genes seen in either dataset (for logging)
      n_shared           number of genes shared by both datasets before panel selection
      n_exclusive        number of panel genes present in only one dataset (see the
                         warning below - these act as a perfect population separator)
    """
    data_dir = args.data_dir if os.path.isabs(args.data_dir) \
        else os.path.join(_REPO_ROOT, args.data_dir)
    d1 = load_matrix(os.path.join(data_dir, args.file1), args.file1, args.reporter_handling)
    d2 = load_matrix(os.path.join(data_dir, args.file2), args.file2, args.reporter_handling)
    print(f'{args.file1}: {d1.shape[0]} cells x {d1.shape[1]} genes')
    print(f'{args.file2}: {d2.shape[0]} cells x {d2.shape[1]} genes')

    g1, g2 = set(d1.columns), set(d2.columns)
    shared = g1 & g2
    reporters_present = sorted(g for g in (g1 | g2) if g in EXCLUDE_GENES)
    drop_set = set() if args.reporter_handling == 'published' else EXCLUDE_GENES
    space = sorted((shared if args.gene_space == 'shared' else (g1 | g2)) - drop_set)
    shared_clean = shared - drop_set
    print(f'{args.gene_space} gene space: {len(space)} genes '
          f'({len(g1)} + {len(g2)}, {len(shared)} shared; '
          f'reporter/plasmid genes present: {reporters_present})')

    p1 = d1.loc[d1.sum(axis=1).sort_values(ascending=False).index[:args.n_pool]]
    p2 = d2.loc[d2.sum(axis=1).sort_values(ascending=False).index[:args.n_pool]]
    p1 = p1.reindex(columns=space, fill_value=0.0)
    p2 = p2.reindex(columns=space, fill_value=0.0)

    combined = pd.concat([p1, p2]).to_numpy(dtype=float)
    mean, var = combined.mean(axis=0), combined.var(axis=0)
    # Fano factor (var/mean) on the pooled cells: picks genes that vary a lot relative
    # to their own level, which is what a correlation-based metric like GMP-Cor is
    # sensitive to. Computed once on the combined pool and then held fixed at every
    # ratio, so the panel itself cannot be responsible for any ratio-dependent effect.
    with np.errstate(divide='ignore', invalid='ignore'):
        fano = np.where(mean > 0, var / mean, 0.0)
    n_genes = min(args.n_genes, combined.shape[1])
    idx = np.sort(np.argsort(fano)[::-1][:n_genes])   # top-Fano genes, back in name order
    genes = [space[i] for i in idx]
    n_exclusive = sum(1 for g in genes if g not in shared_clean)

    a, b = p1.to_numpy(dtype=float)[:, idx], p2.to_numpy(dtype=float)[:, idx]
    print(f'cell pools: {a.shape[0]} + {b.shape[0]}; fixed panel of {len(genes)} genes')
    if n_exclusive:
        print(f'! {n_exclusive} of {len(genes)} selected genes ({n_exclusive / len(genes):.0%}) '
              f'are present in only ONE dataset and are zero-filled in the other. They '
              f'separate the two populations perfectly by construction, so they inflate '
              f'the mixture GMP-Cor for a reason that is not biological.')
    if args.n_total > min(a.shape[0], b.shape[0]):
        print(f'! n_total={args.n_total} exceeds the smaller pool '
              f'({min(a.shape[0], b.shape[0])}); pure endpoints will use the whole pool '
              f'and their spread will reflect the scramble draw only')
    return a, b, genes, reporters_present, len(shared), n_exclusive


# ── Run ──────────────────────────────────────────────────────────────────────

def run(a, b, args):
    """Compute GMP-Cor at every ratio x repeat, and its group-centered counterpart
    for mixtures. Returns a tidy long-form DataFrame, one row per (ratio, repeat).
    """
    records = []
    for r in args.ratios:
        if args.all_cells:
            n_a = a.shape[0] if r > 0 else 0
            n_b = b.shape[0] if r < 1 else 0
        else:
            n_a = int(round(args.n_total * r))
            n_b = args.n_total - n_a
        for rep in range(args.repeats):
            # a fresh, independent stream per (ratio, repeat); ratio is folded in via
            # int(r * 100) (e.g. r=0.25 -> 25) so distinct ratios cannot collide, and
            # repeat is added last since it is always < 1000
            rng = np.random.default_rng(args.seed + 1000 * int(r * 100) + rep)
            parts, labels = [], []
            if n_a:
                take = min(n_a, a.shape[0])
                parts.append(a[rng.choice(a.shape[0], size=take, replace=False)])
                labels.append(np.zeros(take, dtype=int))
            if n_b:
                take = min(n_b, b.shape[0])
                parts.append(b[rng.choice(b.shape[0], size=take, replace=False)])
                labels.append(np.ones(take, dtype=int))
            m = np.vstack(parts)
            lab = np.concatenate(labels)

            res = gmp_cor(m, norm=True, norm_sum=args.norm_sum)
            rec = dict(ratio_1=r, repeat=rep, n_cells=int(m.shape[0]),
                       n_from_1=int((lab == 0).sum()), n_from_2=int((lab == 1).sum()),
                       mean_depth=float(m.sum(axis=1).mean()),
                       mean_detected=float((m > 0).sum(axis=1).mean()), **res)
            if (lab == 0).any() and (lab == 1).any():
                # group_centered subtracts each dataset's OWN per-gene mean before
                # normalizing, which removes the between-population mean shift but
                # keeps within-cell structure; norm=False because group_centered has
                # already normalized. The drop from raw to centered GMP-Cor (d_gmp)
                # is therefore the share of the mixture's signal that is purely the
                # two populations sitting apart, not new gene-gene coordination.
                cen = gmp_cor(group_centered(m, lab, norm_sum=args.norm_sum), norm=False)
                rec['gmp_cor_group_centered'] = cen['gmp_cor']
                rec['d_gmp'] = 1 - cen['gmp_cor'] / res['gmp_cor'] if res['gmp_cor'] else np.nan
                rec.update({k: v for k, v in
                            separation_metrics(m, lab, norm_sum=args.norm_sum).items()
                            if k in ('auc_group_axis', 'silhouette', 'bimodality_coef',
                                     'group_mode_index', 'group_mode_eigenvalue')})
            records.append(rec)
            print(f'  ratio={r:<5} rep={rep}  n={m.shape[0]:<5} '
                  f'GMP-Cor={res["gmp_cor"]:7.2f}  p={res["p_kept"]}  '
                  f'lam*_scr={res["lambda_max_scrambled"]:.2f}')
    return pd.DataFrame(records)


def summarize(df):
    """Collapse the per-repeat records to per-ratio mean/std of every numeric column."""
    num = df.select_dtypes(include=[np.number]).columns.drop(['repeat'])
    g = df.groupby('ratio_1')[list(num)]
    out = g.mean()
    out.columns = [f'{c}_mean' for c in out.columns]
    sd = g.std()
    sd.columns = [f'{c}_std' for c in sd.columns]
    return pd.concat([out, sd], axis=1).reset_index()


# ── Figure ───────────────────────────────────────────────────────────────────

def make_figure(df, summary, args, path_svg, path_png):
    """Two-panel figure: GMP-Cor vs mixing ratio (left), and raw vs group-centered
    GMP-Cor at each mixture with the dGMP fraction annotated (right)."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    fig.subplots_adjust(wspace=0.3, bottom=0.16)
    fs = 11

    ax = axes[0]
    ax.errorbar(summary['ratio_1'], summary['gmp_cor_mean'],
                yerr=summary['gmp_cor_std'], marker='o', capsize=3,
                color='#b2182b', lw=1.4)
    lo = summary[summary['ratio_1'].isin([0.0, 1.0])]['gmp_cor_mean']
    if len(lo):
        ax.axhspan(lo.min(), lo.max(), color='0.85', alpha=0.6, zorder=0)
        ax.text(0.5, lo.max(), ' range spanned by the two pure populations',
                fontsize=fs - 3, va='bottom', ha='center', color='0.35')
    ax.set_xlabel(f'fraction of cells from {args.label1}', fontsize=fs)
    ax.set_ylabel('GMP-Cor', fontsize=fs)
    ax.set_title(f'{args.label1} + {args.label2}\n'
                 f'n = {args.n_total} cells and {args.n_genes} genes at every point',
                 fontsize=fs - 1)
    ax.tick_params(labelsize=fs - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax = axes[1]
    mix = summary[(summary['ratio_1'] > 0) & (summary['ratio_1'] < 1)]
    x = np.arange(len(mix))
    ax.bar(x - 0.2, mix['gmp_cor_mean'], yerr=mix['gmp_cor_std'], width=0.4,
           color='#b2182b', capsize=3, label='raw')
    ax.bar(x + 0.2, mix['gmp_cor_group_centered_mean'],
           yerr=mix['gmp_cor_group_centered_std'], width=0.4,
           color='#4393c3', capsize=3, label='group-centered')
    for xi, (_, r) in zip(x, mix.iterrows()):
        ax.text(xi, max(r['gmp_cor_mean'], r['gmp_cor_group_centered_mean']) * 1.04,
                f'dGMP = {r["d_gmp_mean"]:.2f}', ha='center', fontsize=fs - 3)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{r:.0%} / {1 - r:.0%}' for r in mix['ratio_1']],
                       fontsize=fs - 2)
    ax.set_xlabel(f'{args.label1} / {args.label2}', fontsize=fs)
    ax.set_ylabel('GMP-Cor', fontsize=fs)
    ax.set_title('Removing the between-population means\nreturns the mixture to '
                 'the pure-population level', fontsize=fs - 1)
    ax.tick_params(axis='y', labelsize=fs - 2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=fs - 3, frameon=False)

    fig.savefig(path_svg, bbox_inches='tight')
    fig.savefig(path_png, dpi=200, bbox_inches='tight')
    plt.close(fig)


def write_summary(summary, args, meta, paths, timestamp):
    """Build the human-readable .txt report (design, results table, interpretation)."""
    out = []
    w = out.append
    w('=' * 70)
    w('EXPERIMENTAL DATASET MIXING SERIES -- GMP-Cor with repeats')
    w(f'Run timestamp : {timestamp}')
    w(f'Script        : {os.path.abspath(__file__)}')
    w('=' * 70)
    w('')
    w('DATA')
    w('-' * 70)
    w(f'  dataset 1 (ratio 1.0) : {args.file1}   [{args.label1}]')
    w(f'  dataset 2 (ratio 0.0) : {args.file2}   [{args.label2}]')
    w(f'  data_dir              : {args.data_dir}')
    w(f'  gene space            : {args.gene_space.upper()}')
    verb = 'RETAINED (published convention)' if args.reporter_handling == 'published'         else 'excluded'
    w(f'  shared genes          : {meta["n_shared"]}')
    w(f'  reporter/plasmid genes {verb}: {meta["reporters_present"]}')
    if args.reporter_handling == 'published':
        w('     (the paper drop list is exact-case, so lower-cased reporter genes such as')
        w('      `laci` -- 38% of counts in Expira_biorep_t0A -- survive into the')
        w('      published per-sample values this run is meant to be comparable with)')
    if meta['n_exclusive_in_panel']:
        frac = meta['n_exclusive_in_panel'] / meta['n_genes']
        w(f'  ** {meta["n_exclusive_in_panel"]} of the {meta["n_genes"]} selected genes '
          f'({frac:.0%}) are present in only ONE')
        w('     dataset and are zero-filled in the other. Such genes separate the two')
        w('     populations perfectly by construction and are preferentially picked by')
        w('     a Fano cut, so they inflate the mixture GMP-Cor non-biologically. **')
    w(f'  cell pools            : {meta["n_pool_1"]} + {meta["n_pool_2"]} '
      f'(highest total counts)')
    w(f'  fixed gene panel      : {meta["n_genes"]} genes, top Fano on the combined pool')
    w('')
    w('SAMPLING')
    w('-' * 70)
    if args.all_cells:
        w('  n per point       : ALL cells of each pool (--all-cells), so n DIFFERS')
        w('                      between the pure points and the mixtures. A larger n')
        w('                      lowers the Marchenko-Pastur edge and a larger gene')
        w('                      panel raises GMP-Cor (it is extensive in p), and both')
        w('                      act in favour of the mixture -- compare the p and')
        w('                      lam*_scr columns before reading the ratio as an effect.')
    else:
        w(f'  n_total per point : {args.n_total}   (same at EVERY ratio)')
    w(f'  reporter handling : {args.reporter_handling}')
    w(f'  ratios            : {args.ratios}')
    if args.all_cells:
        w(f'  repeats           : {args.repeats}; with --all-cells the cell set is fixed,')
        w('                      so the SD measures the scramble draw only, not cell sampling')
    else:
        w(f'  repeats           : {args.repeats}, uniform without replacement, seeded')
    w(f'  GMP-Cor           : get_eig_dist(norm=True, norm_method="sum", '
      f'norm_sum={args.norm_sum})')
    w('')
    w('RESULTS (mean +/- SD over repeats)')
    w('-' * 70)
    w(f'  {"ratio":>6} {"GMP-Cor":>17} {"centered":>17} {"dGMP":>7} {"AUC":>6} '
      f'{"p":>6} {"lam*_scr":>9}')
    for _, r in summary.iterrows():
        cen = ('-' if pd.isna(r.get('gmp_cor_group_centered_mean'))
               else f'{r["gmp_cor_group_centered_mean"]:8.2f} +/- '
                    f'{r["gmp_cor_group_centered_std"]:<5.2f}')
        dg = '-' if pd.isna(r.get('d_gmp_mean')) else f'{r["d_gmp_mean"]:7.3f}'
        auc = '-' if pd.isna(r.get('auc_group_axis_mean')) else f'{r["auc_group_axis_mean"]:6.3f}'
        w(f'  {r["ratio_1"]:>6.2f} {r["gmp_cor_mean"]:8.2f} +/- {r["gmp_cor_std"]:<5.2f} '
          f'{cen:>17} {dg:>7} {auc:>6} {int(r["p_kept_mean"]):>6} '
          f'{r["lambda_max_scrambled_mean"]:>9.2f}')
    w('')
    w('KEY COMPARISON')
    w('-' * 70)
    pure = summary[summary['ratio_1'].isin([0.0, 1.0])]
    mix = summary[(summary['ratio_1'] > 0) & (summary['ratio_1'] < 1)]
    if len(pure) and len(mix):
        p_hi = pure['gmp_cor_mean'].max()
        for _, r in pure.iterrows():
            lbl = args.label1 if r['ratio_1'] == 1.0 else args.label2
            w(f'  pure {lbl:<24}: {r["gmp_cor_mean"]:7.2f} +/- {r["gmp_cor_std"]:.2f}')
        for _, r in mix.iterrows():
            w(f'  mixture {r["ratio_1"]:.0%}/{1 - r["ratio_1"]:.0%}'.ljust(31) +
              f': {r["gmp_cor_mean"]:7.2f} +/- {r["gmp_cor_std"]:.2f}   '
              f'({r["gmp_cor_mean"] / p_hi:.2f}x the larger pure population)')
        w('')
        w(f'  Every mixture exceeds both pure populations: '
          f'{bool((mix["gmp_cor_mean"] > p_hi).all())}')
    w('')
    w('INTERPRETATION')
    w('-' * 70)
    all_above = bool(len(pure) and len(mix)
                     and (mix['gmp_cor_mean'] > pure['gmp_cor_mean'].max()).all())
    if all_above:
        w('  Mixing two experimentally distinct populations RAISES GMP-Cor above both')
        w('  of the populations being mixed, and group-mean centering removes most of')
        w('  the excess -- i.e. the elevation is the between-population mode, not new')
        w('  within-cell coordination. Population heterogeneity therefore cannot')
        w('  produce a LOW GMP-Cor, which is what the Dis-Arrest interpretation')
        w('  requires. This is the experimental counterpart of the simulated')
        w('  inverted-sub-population scenario in')
        w('  simulations/inverted_subpopulation_mixing_run.py.')
    else:
        w('  NOT every mixture here exceeds both pure populations, so this run does')
        w('  NOT on its own support the claim that mixing raises GMP-Cor. Before')
        w('  reading anything into that, check the two things that most often cause')
        w('  it: too few cells or too few genes for a between-population mode to')
        w('  clear the noise threshold, and populations that are not actually')
        w('  distinct in this gene panel (see the group-axis AUC column -- a value')
        w('  near 0.5 means the two datasets do not separate at all).')
    w('')
    w('  dGMP is the diagnostic to quote: it is the fraction of a mixture GMP-Cor')
    w('  that disappears when each population is centered on its own gene means,')
    w('  i.e. the fraction that is between-population structure rather than')
    w('  within-cell coordination.')
    w('')
    w('FILES')
    w('-' * 70)
    for k, v in paths.items():
        w(f'  {k:<5}: {v}')
    w('=' * 70)
    return '\n'.join(out)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    """Build the cell pools and fixed gene panel, run the ratio sweep, and write outputs."""
    args = parse_args()
    for d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
        os.makedirs(d, exist_ok=True)
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    base = f'{STEM}_{timestamp}'
    paths = dict(json=os.path.join(_LOG_DIR, base + '.json'),
                 txt=os.path.join(_RAW_DIR, base + '.txt'),
                 csv=os.path.join(_RAW_DIR, base + '.csv'),
                 svg=os.path.join(_FIG_DIR, base + '.svg'),
                 png=os.path.join(_FIG_DIR, base + '.png'))

    print('=' * 65)
    print('Experimental dataset mixing series')
    print('=' * 65)
    a, b, genes, reporters_present, n_shared, n_exclusive = build_pools(args)
    meta = dict(n_shared=n_shared, reporters_present=reporters_present, n_genes=len(genes),
                n_exclusive_in_panel=n_exclusive,
                n_pool_1=int(a.shape[0]), n_pool_2=int(b.shape[0]))

    if args.all_cells and args.repeats > 1:
        print('note: --all-cells fixes the cell set, so the repeats resample only the '
              'scramble draw; their SD is metric noise, not cell-sampling variability')
    df = run(a, b, args)
    summary = summarize(df)

    df.to_csv(paths['csv'], index=False)
    make_figure(df, summary, args, paths['svg'], paths['png'])
    text = write_summary(summary, args, meta, paths, datetime.datetime.now().isoformat())
    with open(paths['txt'], 'w', encoding='utf-8') as fh:
        fh.write(text + '\n')

    with open(paths['json'], 'w', encoding='utf-8') as fh:
        json.dump(dict(
            script=os.path.abspath(__file__), timestamp=timestamp,
            params=vars(args), meta=meta, genes=genes,
            gmp_cor_definition='sum(max(lambda_i - max_scrambled_lambda, 0))',
            per_repeat=df.to_dict(orient='records'),
            per_ratio=summary.to_dict(orient='records'),
        ), fh, indent=2, default=float)

    print('\n' + text)


if __name__ == '__main__':
    main()
