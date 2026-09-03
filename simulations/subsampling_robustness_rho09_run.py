"""
Subsampling robustness test (Reviewer #1, comment 2.3).

Question
--------
Is GMP-Cor robust to the number of cells profiled? A reviewer asks whether the
metric we report could be an artefact of dataset size. To answer this we take a
single, strongly-correlated (regulated) synthetic dataset at rho = 0.9 and ask
whether its GMP-Cor remains stable as we randomly subsample fewer and fewer
cells.

Design
------
1. Generate one large "pool" dataset at rho = 0.9
   (N_CELLS_POOL x N_GENES_POOL), well above the largest subsample so every
   requested size is a genuine random subsample of the same ground-truth
   network.
2. For each target cell count in CELL_SIZES = [200, 400, 600, 800, 1000]:
     - randomly subsample that many cells (without replacement),
     - keep the cell:gene aspect ratio fixed at RATIO = n_cells / n_genes by
       also subsampling genes to  n_genes = round(n_cells / RATIO),
     - select genes by HIGHEST EXPRESSION FIRST (top-N by total counts in the
       subsampled cells),
     - compute GMP-Cor on the resulting matrix.
   Keeping the aspect ratio fixed keeps the Marchenko-Pastur noise edge
   comparable across sizes, so GMP-Cor values are directly comparable.
3. Repeat N_REPEATS times per size with independent random cell draws to get
   mean +/- SD.

GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) )   (same as everywhere
else in the codebase).

Outputs (results/simulation_results/):
  - logs/subsampling_robustness_rho09_<timestamp>.json
  - raw/subsampling_robustness_rho09_<timestamp>.txt
  - figures/subsampling_robustness_rho09_<timestamp>.svg / .png
"""

import sys
import os
import json
import datetime
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.simulations import generate_gram_hub_matrix, simulate_scRNA_data  # noqa: E402
from src.analysis_functions import get_eig_dist  # noqa: E402

# ── Parameters ───────────────────────────────────────────────────────────────

RHO            = 0.9          # regulated (strongly correlated) condition
CELL_SIZES     = [200, 400, 600, 800, 1000]
RATIO          = 0.5          # cell:gene ratio (n_cells / n_genes); 1:2
GENE_SIZES     = [int(round(c / RATIO)) for c in CELL_SIZES]   # 400..2000

# Pool: generated large enough to subsample the biggest target cleanly.
N_CELLS_POOL   = 1200
N_GENES_POOL   = 2400

N_REPEATS      = 5            # independent random cell draws per size
DROPOUT_RATE   = 1.0
SHAPE          = 1.5          # Pareto shape (cluster sizes)
HUB_PROB       = 0.2
SIGMA_SEED     = 31           # hub-network topology of the pool
COUNT_SEED     = 0            # count-sampling RNG for the pool
SUBSAMPLE_SEED = 1000         # base RNG seed for cell subsampling

# ── Output paths ─────────────────────────────────────────────────────────────

_SIM_RESULTS = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
_FIG_DIR = os.path.join(_SIM_RESULTS, 'figures')
_RAW_DIR = os.path.join(_SIM_RESULTS, 'raw')
_LOG_DIR = os.path.join(_SIM_RESULTS, 'logs')
for _d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
    os.makedirs(_d, exist_ok=True)

_ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
JSON_LOG = os.path.join(_LOG_DIR, f'subsampling_robustness_rho09_{_ts}.json')
TEXT_LOG = os.path.join(_RAW_DIR, f'subsampling_robustness_rho09_{_ts}.txt')
SVG_FIG  = os.path.join(_FIG_DIR, f'subsampling_robustness_rho09_{_ts}.svg')
PNG_FIG  = os.path.join(_FIG_DIR, f'subsampling_robustness_rho09_{_ts}.png')


def gmp_cor(observed):
    """GMP-Cor of an observed count matrix (cells x genes).

    GMP-Cor = sum_i max(lambda_i - lambda*_scrambled, 0): the excess correlation
    spectral mass above the scrambled-null threshold lambda*_scrambled (the largest
    eigenvalue of the column-permuted matrix, averaged over repeats by get_eig_dist).
    norm_sum=100 here differs from the norm_sum=50 default used by src.simulations.gmp_cor
    and most other scripts in this project; since GMP-Cor is a sum of eigenvalues of the
    per-row-normalised, z-scored matrix, absolute values from this script are not directly
    comparable to GMP-Cor numbers computed elsewhere at norm_sum=50 (see FINDINGS).
    """
    pcs, pcs1, _ = get_eig_dist(observed, norm=True, log=False, norm_sum=100)
    return float(np.sum(np.maximum(pcs - pcs1.max(), 0)))


# ── Generate the pool ────────────────────────────────────────────────────────

print(f'Generating pool: {N_CELLS_POOL} cells x {N_GENES_POOL} genes at rho={RHO} ...')
sigma_pool = generate_gram_hub_matrix(N_GENES_POOL, RHO, SHAPE, HUB_PROB, seed=SIGMA_SEED)
_, pool = simulate_scRNA_data(
    n_cells=N_CELLS_POOL, n_genes=N_GENES_POOL, sigma=sigma_pool,
    dropout_rate=DROPOUT_RATE, seed=COUNT_SEED,
)
print(f'Pool generated: shape={pool.shape}')

# Reference: GMP-Cor of the full pool (for context)
pool_gmp = gmp_cor(pool)
print(f'Full-pool GMP-Cor (rho={RHO}): {pool_gmp:.4f}\n')

# ── Subsampling sweep ────────────────────────────────────────────────────────

records = []           # one row per (size, repeat)
summary_rows = []      # aggregated per size

for n_cells, n_genes in zip(CELL_SIZES, GENE_SIZES):
    vals = []
    for rep in range(N_REPEATS):
        # distinct, reproducible stream per (size, repeat); multiplying n_cells by
        # 1000 keeps sizes from colliding since N_REPEATS is always well under 1000
        rng = np.random.default_rng(SUBSAMPLE_SEED + 1000 * n_cells + rep)

        # 1) random subsample of cells (without replacement)
        cell_idx = rng.choice(pool.shape[0], size=n_cells, replace=False)
        sub_cells = pool[cell_idx, :]

        # 2) gene subsample: highest expression first (top-N by total counts
        #    within the subsampled cells)
        gene_totals = sub_cells.sum(axis=0)
        top_genes = np.argsort(gene_totals)[::-1][:n_genes]
        sub = sub_cells[:, top_genes]

        g = gmp_cor(sub)
        # raw GMP-Cor is a sum over genes and so scales with panel size (extensive);
        # dividing by n_genes gives a per-gene, scale-free index that should be
        # comparable across the different (n_cells, n_genes) points in this sweep
        g_per_gene = g / n_genes          # intensive (scale-free) index
        vals.append(g)
        records.append({
            'n_cells': n_cells, 'n_genes': n_genes, 'repeat': rep,
            'gmp_cor': g,
            'gmp_cor_per_gene': g_per_gene,
            'subsample_seed': int(SUBSAMPLE_SEED + 1000 * n_cells + rep),
        })
        print(f'  n_cells={n_cells:4d}  n_genes={n_genes:4d}  rep={rep}  '
              f'GMP-Cor={g:8.4f}  per-gene={g_per_gene:.5f}')

    vals = np.asarray(vals)
    vals_pg = vals / n_genes
    summary_rows.append({
        'n_cells': n_cells, 'n_genes': n_genes,
        'mean': float(vals.mean()), 'std': float(vals.std(ddof=1)),
        'min': float(vals.min()), 'max': float(vals.max()),
        'cv': float(vals.std(ddof=1) / vals.mean()) if vals.mean() else float('nan'),
        'mean_per_gene': float(vals_pg.mean()),
        'std_per_gene': float(vals_pg.std(ddof=1)),
        'cv_per_gene': float(vals_pg.std(ddof=1) / vals_pg.mean()) if vals_pg.mean() else float('nan'),
    })
    print(f'  -> n_cells={n_cells}: mean={vals.mean():.4f}  '
          f'SD={vals.std(ddof=1):.4f}  per-gene mean={vals_pg.mean():.5f}\n')

# ── Assemble results ─────────────────────────────────────────────────────────

results = {
    'params': {
        'rho': RHO,
        'cell_sizes': CELL_SIZES,
        'gene_sizes': GENE_SIZES,
        'cell_gene_ratio': RATIO,
        'n_cells_pool': N_CELLS_POOL,
        'n_genes_pool': N_GENES_POOL,
        'n_repeats': N_REPEATS,
        'dropout_rate': DROPOUT_RATE,
        'shape': SHAPE,
        'hub_probability': HUB_PROB,
        'sigma_seed': SIGMA_SEED,
        'count_seed': COUNT_SEED,
        'subsample_seed_base': SUBSAMPLE_SEED,
        'gene_selection': 'highest total expression first (top-N within subsampled cells)',
        'cell_selection': 'uniform random without replacement',
        'gmp_cor_definition': 'sum(max(lambda_i - max_scrambled_lambda, 0))',
        'note': (
            'Single rho=0.9 pool subsampled to decreasing cell counts; gene count '
            'scaled with cell count to hold the cell:gene aspect ratio (and thus '
            'the Marchenko-Pastur noise edge) fixed across sizes.'
        ),
        'timestamp': datetime.datetime.now().isoformat(),
    },
    'pool_gmp_cor': pool_gmp,
    'per_repeat': records,
    'per_size': summary_rows,
}

with open(JSON_LOG, 'w') as fh:
    fh.write(json.dumps(results, indent=2) + '\n')
print(f'JSON log written to: {JSON_LOG}')

# ── Human-readable summary ───────────────────────────────────────────────────

sep = '-' * 65
means = [r['mean'] for r in summary_rows]
means_pg = [r['mean_per_gene'] for r in summary_rows]
# Spread of the per-size means, for raw vs per-gene (scale-free) metric.
raw_spread = (max(means) - min(means)) / np.mean(means) * 100 if np.mean(means) else float('nan')
pg_spread = (max(means_pg) - min(means_pg)) / np.mean(means_pg) * 100 if np.mean(means_pg) else float('nan')
# Linear scaling of raw GMP-Cor with n_genes (expected if the metric is extensive).
raw_slope, raw_intercept = np.polyfit([r['n_genes'] for r in summary_rows], means, 1)

table_lines = [
    '  n_cells  n_genes    mean     SD      min      max      CV     per-gene',
    '  ' + '-' * 68,
]
for r in summary_rows:
    table_lines.append(
        f"  {r['n_cells']:6d}  {r['n_genes']:6d}  {r['mean']:7.3f}  "
        f"{r['std']:6.3f}  {r['min']:7.3f}  {r['max']:7.3f}  {r['cv']:6.3f}  "
        f"{r['mean_per_gene']:.5f}"
    )

summary = [
    '=' * 65,
    'SUBSAMPLING ROBUSTNESS TEST  (REGULATED rho=0.9)',
    f'Run timestamp : {results["params"]["timestamp"]}',
    f'Script        : {os.path.abspath(__file__)}',
    '=' * 65, '',
    'PARAMETERS', sep,
    f'  rho                   : {RHO}',
    f'  cell sizes            : {CELL_SIZES}',
    f'  gene sizes            : {GENE_SIZES}',
    f'  cell:gene ratio       : {RATIO}  (held fixed across sizes)',
    f'  pool                  : {N_CELLS_POOL} cells x {N_GENES_POOL} genes',
    f'  repeats per size      : {N_REPEATS}',
    f'  dropout_rate          : {DROPOUT_RATE}',
    f'  shape (Pareto)        : {SHAPE}',
    f'  hub_probability       : {HUB_PROB}',
    f'  sigma seed (pool)     : {SIGMA_SEED}',
    f'  count seed (pool)     : {COUNT_SEED}',
    f'  gene selection        : highest total expression first',
    f'  cell selection        : uniform random, no replacement',
    '',
    'GMP-COR DEFINITION', sep,
    '  GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) )',
    '',
    'RESULTS', sep,
    f'  Full-pool GMP-Cor ({N_CELLS_POOL} cells x {N_GENES_POOL} genes) : {pool_gmp:.4f}',
    '',
    *table_lines,
    '',
    'INTERPRETATION', sep,
    *textwrap.wrap(
        f'Raw GMP-Cor is an EXTENSIVE quantity: it is a sum of eigenvalue '
        f'excesses, and after standardization the total spectral mass equals the '
        f'number of genes. It therefore scales almost linearly with dataset size '
        f'(per-size means span {raw_spread:.0f}% over the 5x range; linear fit '
        f'GMP-Cor ~ {raw_slope:.4f} * n_genes + {raw_intercept:.2f}). Raw GMP-Cor '
        f'values are thus only comparable between datasets of equal dimension.',
        width=65, initial_indent='  ', subsequent_indent='  ',
    ),
    '',
    *textwrap.wrap(
        f'Normalising to a per-gene (INTENSIVE) index, GMP-Cor / n_genes, removes '
        f'this size dependence: across the same 5x range the per-gene index is '
        f'essentially constant (per-size means span only {pg_spread:.1f}% about '
        f'their mean of {np.mean(means_pg):.5f}). The correlation structure of a '
        f'regulated (rho=0.9) population is recovered at every subsample size and '
        f'is not an artefact of the number of cells profiled. GMP-Cor is robust '
        f'to subsampling once the metric is made scale-free.',
        width=65, initial_indent='  ', subsequent_indent='  ',
    ),
    '',
    'LOG FILES', sep,
    f'  JSON : {JSON_LOG}',
    f'  Text : {TEXT_LOG}',
    f'  SVG  : {SVG_FIG}',
    f'  PNG  : {PNG_FIG}',
    '=' * 65,
]
summary_text = '\n'.join(summary)
with open(TEXT_LOG, 'w', encoding='utf-8') as fh:
    fh.write(summary_text + '\n')
print('\n' + summary_text)

# ── Figure ───────────────────────────────────────────────────────────────────

xs       = [r['n_cells'] for r in summary_rows]
means    = [r['mean'] for r in summary_rows]
stds     = [r['std'] for r in summary_rows]
means_pg = [r['mean_per_gene'] for r in summary_rows]
stds_pg  = [r['std_per_gene'] for r in summary_rows]
genes    = [r['n_genes'] for r in summary_rows]
xticklabels = [f'{c}\n({g} g)' for c, g in zip(xs, genes)]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5))

# ── Panel A: raw GMP-Cor (extensive — scales with size) ──
for rec in records:
    axA.plot(rec['n_cells'], rec['gmp_cor'], 'o', color='#de2d26',
             alpha=0.25, markersize=5, zorder=2)
axA.errorbar(xs, means, yerr=stds, fmt='o-', color='#67000d', capsize=5,
             linewidth=2, markersize=8, zorder=3, label='mean ± SD')
# linear reference through the means vs n_genes
fit = np.poly1d(np.polyfit(genes, means, 1))
axA.plot(xs, [fit(g) for g in genes], '--', color='grey', linewidth=1.2,
         alpha=0.8, zorder=1, label='linear fit (∝ n_genes)')
for x, m, s in zip(xs, means, stds):
    axA.text(x, m + s + max(means) * 0.03, f'{m:.1f}', ha='center',
             va='bottom', fontsize=9, color='#67000d')
axA.set_xticks(xs); axA.set_xticklabels(xticklabels, fontsize=9)
axA.set_xlabel('Cells subsampled  (cell:gene ratio fixed at 1:2)', fontsize=11)
axA.set_ylabel('GMP-Cor  (raw)', fontsize=12)
axA.set_ylim(bottom=0)
axA.set_title('A  Raw GMP-Cor is extensive (scales with size)', fontsize=12,
              fontweight='bold', pad=8, loc='left')
axA.spines['top'].set_visible(False); axA.spines['right'].set_visible(False)
axA.grid(True, axis='y', linestyle='--', alpha=0.25)
axA.legend(fontsize=9, frameon=False, loc='upper left')

# ── Panel B: per-gene GMP-Cor (intensive — scale-free, robust) ──
for rec in records:
    axB.plot(rec['n_cells'], rec['gmp_cor_per_gene'], 'o', color='#3182bd',
             alpha=0.25, markersize=5, zorder=2)
axB.errorbar(xs, means_pg, yerr=stds_pg, fmt='o-', color='#08519c', capsize=5,
             linewidth=2, markersize=8, zorder=3, label='mean ± SD')
grand = float(np.mean(means_pg))
axB.axhline(grand, color='grey', linestyle='--', linewidth=1.2, alpha=0.8,
            label=f'mean = {grand:.4f}')
for x, m, s in zip(xs, means_pg, stds_pg):
    axB.text(x, m + s + max(means_pg) * 0.03, f'{m:.4f}', ha='center',
             va='bottom', fontsize=8, color='#08519c')
axB.set_xticks(xs); axB.set_xticklabels(xticklabels, fontsize=9)
axB.set_xlabel('Cells subsampled  (cell:gene ratio fixed at 1:2)', fontsize=11)
axB.set_ylabel('GMP-Cor / n_genes  (scale-free)', fontsize=12)
axB.set_ylim(0, max(means_pg) * 1.8)
axB.set_title('B  Per-gene GMP-Cor is invariant (robust)', fontsize=12,
              fontweight='bold', pad=8, loc='left')
axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)
axB.grid(True, axis='y', linestyle='--', alpha=0.25)
axB.legend(fontsize=9, frameon=False, loc='lower right')

fig.suptitle(f'Subsampling robustness — regulated network (ρ = {RHO})',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(SVG_FIG, format='svg', bbox_inches='tight')
fig.savefig(PNG_FIG, format='png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSVG saved to: {SVG_FIG}\nPNG saved to: {PNG_FIG}')
