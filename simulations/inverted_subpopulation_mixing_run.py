"""
Inverted-subpopulation mixing runner -- Reviewer #1, comments 1 and 2.3.

The question
------------
Could the low GMP-Cor we attribute to dysregulation instead be produced by a
mixture of distinct, internally-regulated subpopulations (e.g. growing and
non-growing cells)?

Why a new scenario was needed
-----------------------------
The earlier runners (`subpopulation_mixing_run.py`, `subpopulation_mixing_rho09_run.py`)
gave each sub-population its own hub network but left their marginal expression
profiles as two independent draws from the same prior. Those populations are not
actually distinct in the sense that matters: they barely separate in PCA, so no
between-population mode forms, and mixing them merely dilutes each network with
the other. That understates the reviewer's scenario rather than testing it.

Here sub-population B is built by **inverting sub-population A's expression
ranking** (`invert_gene_means`): A's most lowly-expressed gene is B's most highly
expressed, and so on down the ranking. The multiset of gene means -- and hence the
marginal expression distribution, dynamic range and sparsity of each population --
is exactly preserved, so the two populations are transcriptomically opposite while
remaining statistically identical in every global property. Any separation between
them is therefore a genuine difference in *which* genes are expressed and cannot be
a depth artefact.

Two configurations are run:
  - `shared`   -- both populations use the same hub network, differing only in which
                 genes they express. This is the experimental case: two states of the
                 same organism share regulatory architecture.
  - `distinct` -- each population additionally gets its own network topology. The two
                 networks dilute one another, so this is the conservative case.

Calibration
-----------
`inv_gamma_scale = 0.04` (rather than the 0.01 used elsewhere in this module) is set
so simulated cells detect ~85 genes each, matching the experimental matrices. This
matters: at 0.01 the data are so sparse that only ~6 genes carry any between-population
difference, no separating mode can form, and the scenario cannot be tested at all.

`rho_high = 0.7` is chosen so that a pure sub-population lands at GMP-Cor ~42, inside
the range observed experimentally for regulated samples (~30-50). At rho = 0.9 the
simulated populations sit near 130, far off the experimental scale.

Every ratio is evaluated at the same total cell number, because the scrambled
threshold is a pure function of matrix shape -- comparing a 500-cell pure population
against a 1000-cell mixture is invalid.

Only the 50/50 mixture is run, together with both pure populations (ratios 0 and 1)
as its reference points.

What is measured, at each mixing ratio
--------------------------------------
  GMP-Cor                  = sum( max(lambda_i - lambda*_scrambled, 0) )
  GMP-Cor (group-centered) = the same after subtracting each population's own gene
                             means, which removes all between-population structure
  dGMP                     = 1 - centered/raw: the fraction of GMP-Cor that IS the
                             mixture, rather than within-cell coordination
  separation               = AUC and bimodality along the group axis, and which mode
                             of the spectrum carries it
A single dysregulated population (rho_low) of the same size is the reference for
genuine loss of coordination.

Scenarios (--scenario)
----------------------
  regulated    -- both sub-populations coordinated (rho 0.7 -> pure GMP-Cor ~43);
                 shared and distinct networks; stem `inverted_subpopulation_mixing`
  dysregulated -- both sub-populations with essentially no coupling (rho 0.1 -> pure
                 GMP-Cor ~4); shared network only; stem
                 `inverted_subpopulation_mixing_dysregulated`
In each, `rho_low` sets the single-population reference line to the opposite state.

Outputs (results/simulation_results/, `<stem>_<timestamp>`):
  logs/<stem>_<timestamp>.json        full parameters, per-repeat records, gene-mean
                                      profiles and the example UMAP/PC coordinates
  raw/<stem>_<timestamp>.txt          human-readable summary
  raw/<stem>_<timestamp>.csv          one row per (sigma_mode, repeat, ratio)
  raw/<stem>_<timestamp>_umap.csv     per-cell UMAP/PC/group-axis coordinates
  figures/<stem>_<timestamp>.svg/.png

Re-plotting without re-simulating
---------------------------------
The `_umap.csv` holds every per-cell coordinate of the representative mixture, and the
`.json` holds the same plus the gene-mean profiles and eigenvalues. Panels A-D can all
be rebuilt from those two files alone:

    df   = pd.read_csv('..._<ts>.csv')
    log  = json.load(open('..._<ts>.json'))
    make_figure(df, log['example'], log['gene_means'], 'out.svg', 'out.png')
"""

import os
import sys
import json
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Allow running directly from simulations/ or from repo root
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.simulations import inverted_subpopulation_mixing  # noqa: E402

# ── Parameters ───────────────────────────────────────────────────────────────

PARAMS = dict(
    n_cells=1000,          # total cells at EVERY ratio (pools hold this many each)
    n_genes=2000,          # matches the experimental gene panels
    ratios=(0.0, 0.5, 1.0),   # the 50/50 mixture, plus both pure populations as reference
    rho_high=0.7,          # both sub-populations are internally regulated; 0.7 puts a pure
                           # population at GMP-Cor ~42, inside the experimental range
                           # observed for regulated samples (~30-50)
    rho_low=0.1,           # the dysregulated reference population
    seed_a=20,             # hub-network seed for sub-population A
    seed_b=21,             # hub-network seed for B (used only when sigma_mode='distinct')
    mu_seed=7,             # seed of the expression profile that B inverts
    count_seed_a=0,
    count_seed_b=1,
    dropout_rate=1.0,
    shape=1.5,             # Pareto shape for cluster sizes
    hub_probability=0.2,
    inv_gamma_scale=0.04,  # calibrated to ~85 detected genes/cell (see module docstring)
    norm_sum=50,
)
# Two scenarios, selected with --scenario. They differ only in how strongly each
# sub-population is internally regulated; everything else is held fixed.
#   regulated    — both sub-populations are coordinated (rho 0.7 → pure GMP-Cor ~43)
#   dysregulated — both sub-populations have essentially no gene-gene coupling
#                  (rho 0.1 → pure GMP-Cor ~4, the metric's floor at this data scale:
#                   rho = 0.0 gives the same value, so ~4 is not residual coupling)
# In each case `rho_low` sets the single-population reference line: the opposite state.
SCENARIOS = {
    'regulated': dict(
        rho_high=0.7, rho_low=0.1, sigma_modes=('shared', 'distinct'),
        stem='inverted_subpopulation_mixing',
        ref_label='dysregulated single population'),
    'dysregulated': dict(
        rho_high=0.1, rho_low=0.7, sigma_modes=('shared',),
        stem='inverted_subpopulation_mixing_dysregulated',
        ref_label='regulated single population'),
}
SCENARIO = 'regulated'     # overridden by --scenario
REF_LABEL = SCENARIOS[SCENARIO]['ref_label']

SIGMA_MODES = ('shared', 'distinct')
N_REPEATS = 5
EXAMPLE_RATIO = 0.5        # the mixture shown in the UMAP / bimodality panel
EXAMPLE_MODE = 'shared'

# Experimental reference points, computed on data_for_umap at matched n=1000 and the
# same 2000-gene shared panel (Exp = EXP_biorep_t0A, VapC-2h = VAPC_biorep_t2A).
EXPERIMENTAL = dict(pure_exp=16.52, pure_vapc2h=28.46, mixture_50_50=41.22, d_gmp=0.55)

_SIM_RESULTS = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
_FIG_DIR = os.path.join(_SIM_RESULTS, 'figures')
_RAW_DIR = os.path.join(_SIM_RESULTS, 'raw')
_LOG_DIR = os.path.join(_SIM_RESULTS, 'logs')
STEM = 'inverted_subpopulation_mixing'


# ── Run ──────────────────────────────────────────────────────────────────────

def run():
    """
    Execute the full sweep: for every sigma_mode ('shared'/'distinct' hub network)
    and every repeat, call `inverted_subpopulation_mixing` (src/simulations.py) once,
    which itself evaluates every ratio in PARAMS['ratios'] plus the single-population
    dysregulated reference.

    Returns (df, example, gene_means):
      df         — one row per (sigma_mode, repeat, ratio), plus one reference row per
                   (sigma_mode, repeat) tagged ratio_a='dysregulated_reference'
      example    — per-cell UMAP/PC/group-axis coordinates for the one repeat picked
                   as EXAMPLE_MODE/EXAMPLE_RATIO (used for panel B), else None
      gene_means — the mu_a/mu_b profiles and their Spearman correlation for that same
                   example repeat (used for panel A), else None
    """
    records, example, gene_means = [], None, None
    for mode in SIGMA_MODES:
        for rep in range(N_REPEATS):
            print(f'--- sigma_mode={mode}  repeat {rep + 1}/{N_REPEATS}')
            want_example = (mode == EXAMPLE_MODE and rep == 0)
            res = inverted_subpopulation_mixing(
                sigma_mode=mode, repeat=rep,
                example_ratio=EXAMPLE_RATIO if want_example else None,
                **PARAMS)
            if want_example:
                example = res['example']
                gene_means = res['gene_means']
            for rec in res['ratios']:
                records.append(dict(sigma_mode=mode, repeat=rep, **rec))
            # single-population reference gets its own row, tagged by a sentinel
            # string in the otherwise-numeric ratio_a column so it can be filtered
            # out of / into df alongside the real mixing ratios; 'rho' is dropped
            # since its value differs from ratio_a's per-ratio rho_high/rho_low use
            records.append(dict(sigma_mode=mode, repeat=rep, ratio_a='dysregulated_reference',
                                **{k: v for k, v in res['reference'].items() if k != 'rho'}))
    return pd.DataFrame(records), example, gene_means


def summarize(df):
    """Mean ± SD of every numeric column over repeats, per sigma_mode × ratio."""
    num = df.select_dtypes(include=[np.number]).columns.drop('repeat')
    g = df.groupby(['sigma_mode', 'ratio_a'], dropna=False)[list(num)]
    out = g.agg(['mean', 'std', 'count'])
    return out


# ── Figure ───────────────────────────────────────────────────────────────────

_MODE_COLOR = {'shared': 'steelblue', 'distinct': 'darkorange'}


def make_figure(df, example, gene_means, path_svg, path_png):
    """
    Build the 4-panel summary figure and write it to `path_svg`/`path_png`.

    A: log-log scatter of sub-population A's gene means vs B's, showing the
       inversion (and, via Spearman, that it is an exact rank reversal).
    B: UMAP of the representative 50/50 mixture (`example`), coloured by
       sub-population, with an inset histogram of each population's projection
       onto the group axis (the axis separating the two labels) showing bimodality.
    C: GMP-Cor vs mixing ratio, one line per sigma_mode, against the single-
       population dysregulated-reference line.
    D: dGMP (fraction of GMP-Cor attributable to the between-population mode) at
       the interior ratio(s), compared with the experimental Exp+VapC-2h mixture.
    """
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.26)
    fs = 11

    mix = df[df['ratio_a'] != 'dysregulated_reference'].copy()
    mix['ratio_a'] = mix['ratio_a'].astype(float)
    ref = df[df['ratio_a'] == 'dysregulated_reference']['gmp_cor'].astype(float)

    # A — the inverted expression profile
    ax = fig.add_subplot(gs[0, 0])
    mu_a = np.asarray(gene_means['mu_a'])
    mu_b = np.asarray(gene_means['mu_b'])
    ax.loglog(mu_a, mu_b, '.', color='steelblue', markersize=3, alpha=0.5)
    ax.set_xlabel('gene mean, sub-population A', fontsize=fs)
    ax.set_ylabel('gene mean, sub-population B', fontsize=fs)
    ax.set_title(f"A  Inverted expression profile\n(Spearman = "
                 f"{gene_means['spearman_a_vs_b']:.2f}; identical marginal distribution)",
                 fontsize=fs)
    ax.tick_params(labelsize=fs - 2)

    # B — the mixture separates in UMAP, bimodally along the group axis
    ax = fig.add_subplot(gs[0, 1])
    lab = np.asarray(example['labels'])
    emb = np.asarray(example['umap'])
    for k, (name, color) in enumerate([('sub-pop A', 'steelblue'), ('sub-pop B', 'indianred')]):
        ax.plot(emb[lab == k, 0], emb[lab == k, 1], '.', color=color, markersize=3,
                alpha=0.6, label=name)
    up = example['umap_params']
    ax.set_xlabel('UMAP 1', fontsize=fs)
    ax.set_ylabel('UMAP 2', fontsize=fs)
    ax.set_title(f'B  50/50 mixture separates in UMAP ({EXAMPLE_MODE} network)\n'
                 f"n_neighbors={up['n_neighbors']}, min_dist={up['min_dist']}, "
                 f"{up['n_pcs']} PCs", fontsize=fs)
    ax.legend(fontsize=fs - 2, markerscale=3)
    ax.tick_params(labelsize=fs - 2)

    ins = ax.inset_axes([0.04, 0.10, 0.33, 0.24])
    proj = np.asarray(example['group_axis_projection'])
    bins = np.linspace(proj.min(), proj.max(), 40)
    ins.hist(proj[lab == 0], bins=bins, color='steelblue', alpha=0.7)
    ins.hist(proj[lab == 1], bins=bins, color='indianred', alpha=0.7)
    # label inside the inset — an x-label here would collide with the host axis ticks
    ins.text(0.5, 0.92, 'group axis', transform=ins.transAxes, fontsize=fs - 4,
             ha='center', va='top')
    ins.set_xticks([])
    ins.set_yticks([])

    # C — GMP-Cor vs mixing ratio
    ax = fig.add_subplot(gs[1, 0])
    for mode in SIGMA_MODES:
        color = _MODE_COLOR[mode]
        sub = mix[mix['sigma_mode'] == mode].groupby('ratio_a')['gmp_cor']
        m, sd = sub.mean(), sub.std()
        ax.errorbar(m.index, m.values, yerr=sd.values, marker='o', capsize=3,
                    color=color, label=f'{mode} network')
    ax.axhline(ref.mean(), color='k', ls='--', alpha=0.7)
    ax.text(0.5, ref.mean(), f' {REF_LABEL}', fontsize=fs - 3,
            va='bottom', ha='center', color='k', alpha=0.8)
    ax.set_xlabel('fraction of cells from sub-population A', fontsize=fs)
    ax.set_ylabel('GMP-Cor', fontsize=fs)
    ax.set_title(f'C  GMP-Cor of the 50/50 mixture vs the two\npure {SCENARIO} '
                 f'sub-populations', fontsize=fs)
    ax.legend(fontsize=fs - 2)
    ax.tick_params(labelsize=fs - 2)

    # D — how much of GMP-Cor is the mixture itself
    ax = fig.add_subplot(gs[1, 1])
    interior = mix[(mix['ratio_a'] > 0) & (mix['ratio_a'] < 1)]
    ratios_interior = sorted(interior['ratio_a'].unique())

    if len(ratios_interior) == 1:
        # a single mixing ratio: compare both configurations against the experiment
        labels, values, errors, colors = [], [], [], []
        for mode in SIGMA_MODES:
            color = _MODE_COLOR[mode]
            g = interior[interior['sigma_mode'] == mode]['d_gmp']
            labels.append(f'simulated\n{mode} network')
            values.append(g.mean())
            errors.append(g.std())
            colors.append(color)
        labels.append('experimental\nExp + VapC-2h')
        values.append(EXPERIMENTAL['d_gmp'])
        errors.append(np.nan)
        colors.append('crimson')
        x = np.arange(len(labels))
        ax.bar(x, values, yerr=errors, width=0.6, color=colors, capsize=4)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=fs - 2)
        pct = int(round(ratios_interior[0] * 100))
        ax.set_xlabel(f'{pct}/{100 - pct} mixture', fontsize=fs)
    else:
        width = 0.06
        for k, mode in enumerate(SIGMA_MODES):
            color = _MODE_COLOR[mode]
            sub = interior[interior['sigma_mode'] == mode].groupby('ratio_a')['d_gmp']
            ax.bar(np.asarray(sub.mean().index) + (k - 0.5) * width, sub.mean().values,
                   yerr=sub.std().values, width=width, color=color, capsize=3,
                   label=f'{mode} network')
        ax.axhline(EXPERIMENTAL['d_gmp'], color='crimson', ls='--', alpha=0.8)
        ax.text(0.5, EXPERIMENTAL['d_gmp'], ' experimental Exp+VapC-2h mixture',
                fontsize=fs - 3, va='bottom', ha='center', color='crimson')
        ax.set_xlabel('fraction of cells from sub-population A', fontsize=fs)
        ax.legend(fontsize=fs - 2)
    ax.set_ylabel(r'dGMP = 1 $-$ centered / raw', fontsize=fs)
    ax.set_title('D  Fraction of GMP-Cor carried by the\nbetween-population mode', fontsize=fs)
    ax.tick_params(labelsize=fs - 2)

    fig.savefig(path_svg, bbox_inches='tight')
    fig.savefig(path_png, dpi=200, bbox_inches='tight')
    plt.close(fig)


# ── Text summary ─────────────────────────────────────────────────────────────

def write_summary(df, summary, gene_means, paths, timestamp):
    """Render the human-readable .txt report: parameters, GMP-Cor definition, the
    per-ratio results table, key summary comparisons, the experimental reference
    numbers and a plain-language interpretation section (which differs depending on
    whether SCENARIO is 'regulated' or 'dysregulated'). Returns the report as a
    single string; `summary` (from `summarize`) is accepted but not used directly —
    the per-ratio numbers here are recomputed straight from `df`."""
    mix = df[df['ratio_a'] != 'dysregulated_reference'].copy()
    mix['ratio_a'] = mix['ratio_a'].astype(float)
    ref = df[df['ratio_a'] == 'dysregulated_reference']['gmp_cor'].astype(float)

    lines = []
    w = lines.append
    w('=' * 70)
    w('INVERTED-SUBPOPULATION MIXING — Reviewer #1, comments 1 and 2.3')
    w(f'Run timestamp : {timestamp}')
    w(f'Script        : {os.path.abspath(__file__)}')
    w('=' * 70)
    w('')
    w('PARAMETERS')
    w('-' * 70)
    for k, v in PARAMS.items():
        w(f'  {k:<20}: {v}')
    w(f'  {"sigma_modes":<20}: {SIGMA_MODES}')
    w(f'  {"n_repeats":<20}: {N_REPEATS}')
    w('')
    w('  Sub-population B is sub-population A with the gene expression ranking')
    w('  inverted (lowest-expressed <-> highest-expressed). The multiset of gene')
    w('  means is preserved exactly, so both populations have identical marginal')
    w(f'  expression distributions (Spearman of the two profiles = '
      f'{gene_means["spearman_a_vs_b"]:.3f}).')
    w('')
    w('GMP-COR DEFINITION')
    w('-' * 70)
    w('  GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) )')
    w('  dGMP    = 1 - GMP-Cor(group-centered) / GMP-Cor(raw)')
    w('            i.e. the fraction of GMP-Cor that is between-population structure')
    w('')
    w('RESULTS — GMP-Cor by mixing ratio (mean +/- SD over repeats)')
    w('-' * 70)
    for mode in SIGMA_MODES:
        sub = mix[mix['sigma_mode'] == mode]
        w('')
        w(f'  sigma_mode = {mode}')
        w(f'  {"ratio A":>8} {"GMP-Cor":>16} {"centered":>16} {"dGMP":>8} '
          f'{"AUC":>6} {"bimod":>7} {"mode":>5}')
        for r, g in sub.groupby('ratio_a'):
            auc = g['auc_group_axis'].mean() if 'auc_group_axis' in g else np.nan
            bim = g['bimodality_coef'].mean() if 'bimodality_coef' in g else np.nan
            gm = g['group_mode_index'].mean() if 'group_mode_index' in g else np.nan
            w(f'  {r:>8.2f} {g["gmp_cor"].mean():>9.2f} +/- {g["gmp_cor"].std():<4.1f} '
              f'{g["gmp_cor_group_centered"].mean():>9.2f} +/- {g["gmp_cor_group_centered"].std():<4.1f} '
              f'{g["d_gmp"].mean():>8.3f} {auc:>6.3f} {bim:>7.3f} '
              f'{("-" if np.isnan(gm) else f"{gm + 1:.0f}"):>5}')
    w('')
    w(f'  {REF_LABEL.upper()} (rho={PARAMS["rho_low"]}) : '
      f'{ref.mean():.2f} +/- {ref.std():.2f}')
    w('')
    w('KEY COMPARISONS')
    w('-' * 70)
    for mode in SIGMA_MODES:
        sub = mix[mix['sigma_mode'] == mode]
        pure = sub[sub['ratio_a'].isin([0.0, 1.0])]['gmp_cor']
        inter = sub[(sub['ratio_a'] > 0) & (sub['ratio_a'] < 1)]['gmp_cor']
        w(f'  {mode:<9} pure populations   : {pure.mean():7.2f} +/- {pure.std():.2f}')
        w(f'  {mode:<9} mixtures           : {inter.mean():7.2f} +/- {inter.std():.2f}')
        w(f'  {mode:<9} mixture / pure     : {inter.mean() / pure.mean():7.2f}x')
        w(f'  {mode:<9} lowest mixture     : {inter.min():7.2f}  (the "mixture floor")')
        w(f'  {mode:<9} floor / reference  : {inter.min() / ref.mean():7.2f}x')
        w('')
    w('EXPERIMENTAL REFERENCE (data_for_umap, matched n=1000, shared 2000-gene panel)')
    w('-' * 70)
    w(f'  pure Exp            : {EXPERIMENTAL["pure_exp"]:.2f}')
    w(f'  pure VapC-2h        : {EXPERIMENTAL["pure_vapc2h"]:.2f}')
    w(f'  50/50 mixture       : {EXPERIMENTAL["mixture_50_50"]:.2f}  '
      f'(above BOTH pure populations)')
    w(f'  dGMP of that mixture: {EXPERIMENTAL["d_gmp"]:.2f}')
    w('')
    w('INTERPRETATION')
    w('-' * 70)
    w(f'  Two {SCENARIO} sub-populations that are transcriptomically opposite')
    w('  separate cleanly in UMAP and are bimodal along the group axis - the')
    w('  scenario the reviewer raises.')
    w('')
    if SCENARIO == 'regulated':
        w('  Under GMP-Cor that mixture does not look dysregulated: the separating')
        w('  mode ADDS eigenvalue mass, so the mixture sits ABOVE both pure')
        w('  populations, and far above a dysregulated one. The lowest value')
        w('  reachable by mixing (the mixture floor, above) is the quantity to')
        w('  compare the experimental Dis-Arrest value against.')
    else:
        w('  Both sub-populations are individually at the noise floor, yet their')
        w('  mixture is strongly elevated - the separating mode alone accounts for')
        w('  almost all of it. Two conclusions follow:')
        w('')
        w('  1. Mixing NEVER lowers GMP-Cor. Even mixing two populations that have')
        w('     no internal coordination at all raises the index well above what')
        w('     either shows alone. A low observed GMP-Cor therefore cannot be')
        w('     produced by population heterogeneity - which is exactly the')
        w("     reviewer's concern, answered from the opposite direction.")
        w('  2. The converse does NOT hold, and this bounds the claim: the raw')
        w('     index alone cannot separate a coordinated single population from a')
        w('     mixture of uncoordinated ones, since the two land in the same')
        w('     range. dGMP does separate them (see below), which is why it, and')
        w('     not the scalar, is the diagnostic to report for a clustered sample.')
    w('')
    w('  Group-mean centering separates the two contributions: it removes the')
    w('  between-population mode while leaving within-cell coordination intact.')
    w('  dGMP is therefore a direct test on real data of whether an observed GMP-Cor')
    w('  is mixture structure or genuine coordination — a test that requires')
    w('  single-cell data and cannot be performed on bulk RNA-seq.')
    w('')
    w('FILES')
    w('-' * 70)
    for k, v in paths.items():
        w(f'  {k:<5}: {v}')
    w('=' * 70)
    return '\n'.join(lines)


# ── Main ─────────────────────────────────────────────────────────────────────

def write_umap_csv(example, path):
    """
    Per-cell coordinates of the representative mixture, so the UMAP panel can be
    re-plotted (or re-coloured, or re-clustered) without re-simulating anything.

    One row per cell: sub-population label, UMAP 1/2, the first two principal
    components, and the projection on the group axis (the bimodality inset).
    """
    emb = np.asarray(example['umap'])
    pcs = np.asarray(example['pc_scores'])
    lab = np.asarray(example['labels'])
    pd.DataFrame({
        'cell': np.arange(lab.size),
        'subpopulation': np.where(lab == 0, 'A', 'B'),
        'label': lab,
        'umap_1': emb[:, 0],
        'umap_2': emb[:, 1],
        'pc_1': pcs[:, 0],
        'pc_2': pcs[:, 1],
        'group_axis_projection': np.asarray(example['group_axis_projection']),
    }).to_csv(path, index=False)


def main():
    # entry point: parse --scenario/--repeats, override the module-level scenario
    # settings and PARAMS accordingly, run the full sweep, and write every output
    # (csv, per-cell umap csv, figure, text summary, full json log) under a shared
    # <stem>_<timestamp> basename in results/simulation_results/
    global SCENARIO, REF_LABEL, SIGMA_MODES, STEM, N_REPEATS
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scenario', choices=sorted(SCENARIOS), default='regulated',
                    help='how strongly each sub-population is internally regulated')
    ap.add_argument('--repeats', type=int, default=N_REPEATS)
    args = ap.parse_args()

    SCENARIO = args.scenario
    cfg = SCENARIOS[SCENARIO]
    REF_LABEL = cfg['ref_label']
    SIGMA_MODES = cfg['sigma_modes']
    STEM = cfg['stem']
    N_REPEATS = args.repeats
    PARAMS['rho_high'] = cfg['rho_high']
    PARAMS['rho_low'] = cfg['rho_low']

    for d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
        os.makedirs(d, exist_ok=True)
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    base = f'{STEM}_{timestamp}'
    paths = dict(json=os.path.join(_LOG_DIR, base + '.json'),
                 txt=os.path.join(_RAW_DIR, base + '.txt'),
                 csv=os.path.join(_RAW_DIR, base + '.csv'),
                 umap=os.path.join(_RAW_DIR, base + '_umap.csv'),
                 svg=os.path.join(_FIG_DIR, base + '.svg'),
                 png=os.path.join(_FIG_DIR, base + '.png'))

    print(f'scenario = {SCENARIO}  (rho_high={PARAMS["rho_high"]}, '
          f'rho_low={PARAMS["rho_low"]}, sigma_modes={SIGMA_MODES})')

    df, example, gene_means = run()
    summary = summarize(df)

    df.to_csv(paths['csv'], index=False)
    write_umap_csv(example, paths['umap'])
    make_figure(df, example, gene_means, paths['svg'], paths['png'])
    text = write_summary(df, summary, gene_means, paths,
                         datetime.datetime.now().isoformat())
    with open(paths['txt'], 'w', encoding='utf-8') as fh:
        fh.write(text + '\n')

    with open(paths['json'], 'w', encoding='utf-8') as fh:
        json.dump(dict(
            script=os.path.abspath(__file__),
            timestamp=timestamp,
            scenario=SCENARIO,
            params=PARAMS,
            sigma_modes=list(SIGMA_MODES),
            n_repeats=N_REPEATS,
            example_ratio=EXAMPLE_RATIO,
            example_mode=EXAMPLE_MODE,
            gmp_cor_definition='sum(max(lambda_i - max_scrambled_lambda, 0))',
            experimental_reference=EXPERIMENTAL,
            gene_mean_spearman_a_vs_b=gene_means['spearman_a_vs_b'],
            gene_means=gene_means,
            example=example,   # UMAP / PC coordinates, so the figure can be re-plotted
            records=df.to_dict(orient='records'),
        ), fh, indent=2, default=float)

    print(text)


if __name__ == '__main__':
    main()
