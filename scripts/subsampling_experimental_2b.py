"""
Cell-subsampling scaling of GMP-Cor on experimental data (Reviewer #1, comment 2.3).

Question
--------
The subsampling robustness analysis in `documents/subsampling_analysis.md` was done
entirely on synthetic data, where the ground-truth correlation structure is known.
Does a real dataset show the SAME scaling with cell number, or a different one?

Design — mirrors Experiment 4 of the synthetic analysis
-------------------------------------------------------
Experiment 4 ("all genes retained, cells only") is the clean test: the complete gene
panel is used at every size, so there is no gene-selection axis and the ONLY variable
is the number of cells. This script reproduces that design on
`data_for_paper/sample_2b_filtered.csv` -- the published exponential-growth sample,
1041 cells x 2071 genes, rRNA already removed by the paper pipeline.

  - gene panel : all 2071 genes, held fixed at every subsample size
  - cells      : 500, 600, 700, 800, 900, 1000, drawn uniformly without replacement
  - repeats    : 5 per size, seeded
  - GMP-Cor    : sum( max(lambda_i - lambda*_scrambled, 0) ), via
                 analysis_functions.get_eig_dist(norm=True, norm_method='sum', norm_sum=50)
                 -- identical settings to every other GMP-Cor in this project

Two simulated arms are computed for comparison:

  simulated (calibrated) -- rho 0.7, inv_gamma_scale 0.04, matched to this dataset's
                 dimensions (n_pool = 1041, p = 2071) and to its sparsity. This is the
                 like-for-like control.
  simulated (published)  -- the rho = 0.9 curve from
                 `results/simulation_results/logs/subsampling_allgenes_rho09_*.json`,
                 read from disk, not re-run. Note it used a much sparser expression
                 prior (inv_gamma_scale 0.01, ~24 detected genes per cell versus ~62
                 here), so it is NOT a like-for-like comparison -- it is included
                 because it is the curve the existing write-up reports.

Because GMP-Cor is extensive in the number of genes and its noise threshold is a pure
function of matrix shape, the three arms are compared as ratios to their own value at
n = 1000, not in absolute units.

Outputs (results/subsampling_experimental/):
  logs/subsampling_experimental_2b_<timestamp>.json
  raw/subsampling_experimental_2b_<timestamp>.txt
  raw/subsampling_experimental_2b_<timestamp>.csv
  figures/subsampling_experimental_2b_<timestamp>.svg / .png
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

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.simulations import (gmp_cor, generate_gram_hub_matrix,  # noqa: E402
                             simulate_scRNA_data, draw_gene_means)

# ── Parameters ───────────────────────────────────────────────────────────────

DATA_FILE = os.path.join(_REPO_ROOT, 'data_for_paper', 'sample_2b_filtered.csv')
CELL_SIZES = [500, 600, 700, 800, 900, 1000]
N_REPEATS = 5
SUBSAMPLE_SEED = 1000       # base seed for the cell draws
NORM_SUM = 50
REFERENCE_SIZE = 1000       # size every arm is normalised to

# Non-gene columns to exclude if present (rRNA is already gone from this matrix)
EXCLUDE = {'16s_mature', '16s_unprocessed', 'kanr', 'mcherry', 'gfp', 'laci', 'ampr'}

# Calibrated simulation arm — matched to the experimental matrix
SIM_RHO = 0.7
SIM_INV_GAMMA_SCALE = 0.04
SIM_DROPOUT = 1.0
SIM_SHAPE = 1.5
SIM_HUB_PROB = 0.2
SIM_SIGMA_SEED = 20
SIM_COUNT_SEED = 0
SIM_MU_SEED = 7

# Published synthetic curve to overlay (read, not re-run)
PUBLISHED_SIM_LOG = os.path.join(
    _REPO_ROOT, 'results', 'simulation_results', 'logs',
    'subsampling_allgenes_rho09_20260615_171556.json')

_OUT = os.path.join(_REPO_ROOT, 'results', 'subsampling_experimental')
_FIG_DIR = os.path.join(_OUT, 'figures')
_RAW_DIR = os.path.join(_OUT, 'raw')
_LOG_DIR = os.path.join(_OUT, 'logs')
STEM = 'subsampling_experimental_2b'


# ── Data ─────────────────────────────────────────────────────────────────────

def load_matrix():
    """Read sample_2b_filtered.csv and drop non-gene / reporter columns."""
    d = pd.read_csv(DATA_FILE, index_col=0).fillna(0.0)
    drop = [c for c in d.columns
            if str(c).lower().startswith('unnamed') or str(c).startswith('INTR_')
            or str(c).casefold() in EXCLUDE]
    if drop:
        d = d.drop(columns=drop)
        print(f'dropped {len(drop)} columns: {drop}')
    return d


def subsample_curve(pool, label, sizes, n_repeats, seed_base):
    """GMP-Cor over cell subsamples of a fixed-gene matrix.

    `pool` is cells x genes; the gene panel is never subsampled here, only cells,
    so the resulting curve isolates the effect of cell number alone (Experiment 4
    of the synthetic subsampling analysis, see module docstring).
    """
    records = []
    for n in sizes:
        for rep in range(n_repeats):
            # same collision-free scheme as the synthetic runner: multiplying the
            # size by more than the repeat count keeps every (size, repeat) distinct
            rng = np.random.default_rng(seed_base + 1000 * n + rep)
            idx = rng.choice(pool.shape[0], size=n, replace=False)
            sub = pool[idx]
            res = gmp_cor(sub, norm=True, norm_sum=NORM_SUM)
            records.append(dict(
                arm=label, n_cells=n, repeat=rep,
                mean_depth=float(sub.sum(axis=1).mean()),
                mean_detected=float((sub > 0).sum(axis=1).mean()),
                **res))
            print(f'  {label:<22} n={n:<5} rep={rep}  GMP-Cor={res["gmp_cor"]:8.2f} '
                  f'(p={res["p_kept"]}, lam*_scr={res["lambda_max_scrambled"]:.2f})')
    return records


def simulated_pool(n_cells, n_genes):
    """Generate the calibrated (rho=0.7) synthetic arm, matched to the experimental
    matrix's cell/gene dimensions so its scaling is directly comparable."""
    sigma = generate_gram_hub_matrix(n_genes, SIM_RHO, SIM_SHAPE, SIM_HUB_PROB,
                                     seed=SIM_SIGMA_SEED)
    mu = draw_gene_means(n_genes, seed=SIM_MU_SEED, inv_gamma_scale=SIM_INV_GAMMA_SCALE)
    _, obs = simulate_scRNA_data(n_cells=n_cells, n_genes=n_genes, sigma=sigma,
                                 dropout_rate=SIM_DROPOUT, seed=SIM_COUNT_SEED,
                                 gene_mu=mu)
    return obs


def published_sim_curve():
    """The rho=0.9 synthetic curve from the existing run log.

    Read from disk rather than regenerated, so this arm reflects exactly the run the
    existing write-up cites; it is NOT dimension- or sparsity-matched to sample_2b
    (see the module docstring caveat), so treat it as an illustrative third curve,
    not a controlled comparison.
    """
    if not os.path.exists(PUBLISHED_SIM_LOG):
        print(f'! published sim log not found: {PUBLISHED_SIM_LOG}')
        return []
    log = json.load(open(PUBLISHED_SIM_LOG, encoding='utf8'))
    return [dict(arm='simulated (published, rho=0.9)', n_cells=r['n_cells'],
                 repeat=r['repeat'], gmp_cor=r['gmp_cor'], p_kept=r['n_genes'])
            for r in log['per_repeat']]


# ── Summary ──────────────────────────────────────────────────────────────────

def summarize(df):
    """Per-arm, per-size summary, each arm normalised to its own value at n=1000.

    GMP-Cor is extensive in gene count and its noise threshold depends on matrix
    shape, so the three arms (different gene panels, different generative models)
    are only comparable as a ratio to their own reference point, not in raw units.
    """
    rows = []
    for arm, g_arm in df.groupby('arm', sort=False):
        ref = g_arm[g_arm['n_cells'] == REFERENCE_SIZE]['gmp_cor'].mean()
        for n, g in g_arm.groupby('n_cells'):
            m, sd = g['gmp_cor'].mean(), g['gmp_cor'].std()
            rows.append(dict(
                arm=arm, n_cells=n, n_repeats=len(g),
                gmp_cor_mean=m, gmp_cor_std=sd,
                cv=sd / m if m else np.nan,
                frac_of_reference=m / ref if ref else np.nan,
                p_kept=g['p_kept'].mean() if 'p_kept' in g else np.nan,
                gmp_cor_per_gene=(g['gmp_cor'] / g['p_kept']).mean() if 'p_kept' in g else np.nan,
                lambda_max_scrambled=g['lambda_max_scrambled'].mean() if 'lambda_max_scrambled' in g else np.nan,
                mean_detected=g['mean_detected'].mean() if 'mean_detected' in g else np.nan,
            ))
    return pd.DataFrame(rows)


_ARM_STYLE = {
    'experimental (sample_2b)': dict(color='crimson', marker='o'),
    'simulated (calibrated, rho=0.7)': dict(color='steelblue', marker='s'),
    'simulated (published, rho=0.9)': dict(color='darkorange', marker='^'),
}


def make_figure(summary, path_svg, path_png):
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))
    fig.subplots_adjust(hspace=0.33, wspace=0.27)
    fs = 11

    # A — absolute GMP-Cor, experimental only (arms are not comparable in absolute units)
    ax = axes[0, 0]
    e = summary[summary['arm'] == 'experimental (sample_2b)']
    ax.errorbar(e['n_cells'], e['gmp_cor_mean'], yerr=e['gmp_cor_std'],
                marker='o', capsize=3, color='crimson')
    ax.set_xlabel('cells subsampled', fontsize=fs)
    ax.set_ylabel('GMP-Cor', fontsize=fs)
    ax.set_title('A  Experimental GMP-Cor vs cell number\n(sample_2b, all 2071 genes retained)',
                 fontsize=fs)
    ax.tick_params(labelsize=fs - 2)

    # B — the comparison: each arm scaled to its own n=1000 value
    ax = axes[0, 1]
    for arm, g in summary.groupby('arm', sort=False):
        st = _ARM_STYLE.get(arm, {})
        yerr = (g['gmp_cor_std'] / g['gmp_cor_mean'] * g['frac_of_reference']).values
        ax.errorbar(g['n_cells'], g['frac_of_reference'], yerr=yerr, capsize=3,
                    label=arm, **st)
    ax.axhline(1.0, color='k', ls=':', alpha=0.5)
    ax.set_xlabel('cells subsampled', fontsize=fs)
    ax.set_ylabel(f'GMP-Cor / GMP-Cor at n={REFERENCE_SIZE}', fontsize=fs)
    ax.set_title('B  Scaling with cell number, each arm\nnormalised to its own n=1000',
                 fontsize=fs)
    ax.legend(fontsize=fs - 3)
    ax.tick_params(labelsize=fs - 2)

    # C — run-to-run variability
    ax = axes[1, 0]
    for arm, g in summary.groupby('arm', sort=False):
        st = _ARM_STYLE.get(arm, {})
        ax.plot(g['n_cells'], g['cv'], label=arm, **st)
    ax.set_xlabel('cells subsampled', fontsize=fs)
    ax.set_ylabel('CV over repeats', fontsize=fs)
    ax.set_title('C  Variability across subsample draws', fontsize=fs)
    ax.legend(fontsize=fs - 3)
    ax.tick_params(labelsize=fs - 2)

    # D — what moves with n: the noise threshold and the surviving gene count
    ax = axes[1, 1]
    for arm, g in summary.groupby('arm', sort=False):
        if g['lambda_max_scrambled'].notna().any():
            st = dict(_ARM_STYLE.get(arm, {}))
            ax.plot(g['n_cells'], g['lambda_max_scrambled'], label=arm, **st)
    ax.set_xlabel('cells subsampled', fontsize=fs)
    ax.set_ylabel(r'$\lambda^{\mathrm{scr}}_{\mathrm{max}}$ (noise threshold)', fontsize=fs)
    ax.set_title('D  The scrambled threshold rises as cells are\nremoved (pure matrix shape)',
                 fontsize=fs)
    ax.legend(fontsize=fs - 3)
    ax.tick_params(labelsize=fs - 2)

    fig.savefig(path_svg, bbox_inches='tight')
    fig.savefig(path_png, dpi=200, bbox_inches='tight')
    plt.close(fig)


def write_summary(summary, pool_shape, paths, timestamp):
    """Build the human-readable .txt report (design, per-arm tables, scaling comparison)."""
    lines, w = [], None
    out = []
    w = out.append
    w('=' * 70)
    w('CELL-SUBSAMPLING SCALING OF GMP-COR -- EXPERIMENTAL vs SIMULATED')
    w(f'Run timestamp : {timestamp}')
    w(f'Script        : {os.path.abspath(__file__)}')
    w('=' * 70)
    w('')
    w('DATA')
    w('-' * 70)
    w(f'  file            : {DATA_FILE}')
    w(f'  matrix          : {pool_shape[0]} cells x {pool_shape[1]} genes '
      f'(gene panel held fixed at every size)')
    w(f'  cell sizes      : {CELL_SIZES}')
    w(f'  repeats         : {N_REPEATS}, uniform without replacement, seeded')
    w(f'  GMP-Cor         : get_eig_dist(norm=True, norm_method="sum", norm_sum={NORM_SUM})')
    w('')
    w('CALIBRATED SIMULATION ARM')
    w('-' * 70)
    w(f'  rho={SIM_RHO}, inv_gamma_scale={SIM_INV_GAMMA_SCALE}, dropout={SIM_DROPOUT}, '
      f'shape={SIM_SHAPE}, hub_p={SIM_HUB_PROB}')
    w(f'  sigma_seed={SIM_SIGMA_SEED}, count_seed={SIM_COUNT_SEED}, mu_seed={SIM_MU_SEED}')
    w('  dimensions matched to the experimental matrix')
    w('')
    w('RESULTS')
    w('-' * 70)
    for arm, g in summary.groupby('arm', sort=False):
        w('')
        w(f'  {arm}')
        w(f'  {"cells":>6} {"GMP-Cor":>18} {"CV":>6} {"frac of n=1000":>15} '
          f'{"p":>6} {"lam*_scr":>9} {"detected":>9}')
        for _, r in g.iterrows():
            lam = '-' if pd.isna(r['lambda_max_scrambled']) else f'{r["lambda_max_scrambled"]:.2f}'
            det = '-' if pd.isna(r['mean_detected']) else f'{r["mean_detected"]:.1f}'
            w(f'  {int(r["n_cells"]):>6} {r["gmp_cor_mean"]:>10.2f} +/- {r["gmp_cor_std"]:<5.2f} '
              f'{r["cv"]:>6.3f} {r["frac_of_reference"]:>15.3f} {int(r["p_kept"]):>6} '
              f'{lam:>9} {det:>9}')
    w('')
    w('SCALING COMPARISON -- GMP-Cor as a fraction of the same arm at n=1000')
    w('-' * 70)
    arms = list(summary['arm'].unique())
    common = sorted(set.intersection(*[set(summary[summary['arm'] == a]['n_cells'])
                                       for a in arms]))
    w(f'  {"arm":<34}' + ''.join(f'{f"n={n}":>10}' for n in common))
    for arm in arms:
        g = summary[summary['arm'] == arm].set_index('n_cells')
        w(f'  {arm:<34}' + ''.join(f'{g.loc[n, "frac_of_reference"]:>10.3f}'
                                   for n in common))
    w('')
    w('  Each arm at its own smallest size:')
    for arm in arms:
        g = summary[summary['arm'] == arm]
        lo = g[g['n_cells'] == g['n_cells'].min()]
        w(f'    {arm:<34} n={int(lo["n_cells"].iloc[0]):<5} '
          f'{lo["frac_of_reference"].iloc[0]:.3f}')
    w('')
    w('FILES')
    w('-' * 70)
    for k, v in paths.items():
        w(f'  {k:<5}: {v}')
    w('=' * 70)
    return '\n'.join(out)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    """Run all three arms (experimental, calibrated simulation, published simulation),
    summarise, plot, and write the json/txt/csv/svg outputs."""
    for d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
        os.makedirs(d, exist_ok=True)
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    base = f'{STEM}_{timestamp}'
    paths = dict(json=os.path.join(_LOG_DIR, base + '.json'),
                 txt=os.path.join(_RAW_DIR, base + '.txt'),
                 csv=os.path.join(_RAW_DIR, base + '.csv'),
                 svg=os.path.join(_FIG_DIR, base + '.svg'),
                 png=os.path.join(_FIG_DIR, base + '.png'))

    d = load_matrix()
    pool = d.to_numpy(dtype=float)
    print(f'experimental pool: {pool.shape[0]} cells x {pool.shape[1]} genes')

    records = subsample_curve(pool, 'experimental (sample_2b)', CELL_SIZES,
                              N_REPEATS, SUBSAMPLE_SEED)

    print('\nsimulating calibrated arm matched to the experimental dimensions...')
    sim_pool = simulated_pool(pool.shape[0], pool.shape[1])
    records += subsample_curve(sim_pool, 'simulated (calibrated, rho=0.7)', CELL_SIZES,
                               N_REPEATS, SUBSAMPLE_SEED)

    records += published_sim_curve()

    df = pd.DataFrame(records)
    summary = summarize(df)

    df.to_csv(paths['csv'], index=False)
    make_figure(summary, paths['svg'], paths['png'])
    text = write_summary(summary, pool.shape, paths, datetime.datetime.now().isoformat())
    with open(paths['txt'], 'w', encoding='utf-8') as fh:
        fh.write(text + '\n')

    with open(paths['json'], 'w', encoding='utf-8') as fh:
        json.dump(dict(
            script=os.path.abspath(__file__), timestamp=timestamp,
            data_file=DATA_FILE, pool_shape=list(pool.shape),
            cell_sizes=CELL_SIZES, n_repeats=N_REPEATS,
            subsample_seed=SUBSAMPLE_SEED, norm_sum=NORM_SUM,
            reference_size=REFERENCE_SIZE,
            excluded_columns=sorted(EXCLUDE),
            simulated_arm=dict(rho=SIM_RHO, inv_gamma_scale=SIM_INV_GAMMA_SCALE,
                               dropout_rate=SIM_DROPOUT, shape=SIM_SHAPE,
                               hub_probability=SIM_HUB_PROB, sigma_seed=SIM_SIGMA_SEED,
                               count_seed=SIM_COUNT_SEED, mu_seed=SIM_MU_SEED),
            published_sim_log=PUBLISHED_SIM_LOG,
            gmp_cor_definition='sum(max(lambda_i - max_scrambled_lambda, 0))',
            per_repeat=df.to_dict(orient='records'),
            per_size=summary.to_dict(orient='records'),
        ), fh, indent=2, default=float)

    print('\n' + text)


if __name__ == '__main__':
    main()
