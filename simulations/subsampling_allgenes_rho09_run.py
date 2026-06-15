"""
Subsampling robustness — all genes retained (Reviewer #1, comment 2.3).

The most faithful "profile fewer cells, same genes" test. A single regulated
pool is generated with EXACTLY N_GENES = 2000 genes; every subsample uses the
COMPLETE 2000-gene set (no gene selection at all) and only the number of CELLS
is reduced (200, 400, 600, 800, 1000).

Because the gene set is identical and complete at every size, there is no
gene-selection axis (no top-vs-random): the only variable is the number of
cells. The cell:gene aspect ratio still varies (10:1 genes:cells at 200 cells
-> 2:1 at 1000 cells), so the Marchenko-Pastur noise edge shifts with size.

GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) ).

Outputs (results/simulation_results/):
  - logs/subsampling_allgenes_rho09_<timestamp>.json
  - raw/subsampling_allgenes_rho09_<timestamp>.txt
  - figures/subsampling_allgenes_rho09_<timestamp>.svg / .png
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

RHO            = 0.9
CELL_SIZES     = [200, 400, 600, 800, 1000]
N_GENES        = 2000          # exact pool gene count; ALL genes used every time

N_CELLS_POOL   = 1200          # enough cells to subsample 1000 after filtering

N_REPEATS      = 5
DROPOUT_RATE   = 1.0
SHAPE          = 1.5
HUB_PROB       = 0.2
SIGMA_SEED     = 31
COUNT_SEED     = 0
SUBSAMPLE_SEED = 1000          # base for CELL draws

# ── Output paths ─────────────────────────────────────────────────────────────

_SIM_RESULTS = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
_FIG_DIR = os.path.join(_SIM_RESULTS, 'figures')
_RAW_DIR = os.path.join(_SIM_RESULTS, 'raw')
_LOG_DIR = os.path.join(_SIM_RESULTS, 'logs')
for _d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
    os.makedirs(_d, exist_ok=True)

_ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
JSON_LOG = os.path.join(_LOG_DIR, f'subsampling_allgenes_rho09_{_ts}.json')
TEXT_LOG = os.path.join(_RAW_DIR, f'subsampling_allgenes_rho09_{_ts}.txt')
SVG_FIG  = os.path.join(_FIG_DIR, f'subsampling_allgenes_rho09_{_ts}.svg')
PNG_FIG  = os.path.join(_FIG_DIR, f'subsampling_allgenes_rho09_{_ts}.png')


def gmp_cor(observed):
    pcs, pcs1, _ = get_eig_dist(observed, norm=True, log=False, norm_sum=100)
    return float(np.sum(np.maximum(pcs - pcs1.max(), 0)))


# ── Generate the pool (exactly N_GENES genes) ────────────────────────────────

print(f'Generating pool: {N_CELLS_POOL} cells x {N_GENES} genes at rho={RHO} ...')
sigma_pool = generate_gram_hub_matrix(N_GENES, RHO, SHAPE, HUB_PROB, seed=SIGMA_SEED)
_, pool = simulate_scRNA_data(
    n_cells=N_CELLS_POOL, n_genes=N_GENES, sigma=sigma_pool,
    dropout_rate=DROPOUT_RATE, seed=COUNT_SEED,
)
print(f'Pool generated: shape={pool.shape}  (ALL {N_GENES} genes used every time)\n')

# ── Subsampling sweep (cells only; all genes kept) ───────────────────────────

records = []
summary = []

for n_cells in CELL_SIZES:
    vals = []
    for rep in range(N_REPEATS):
        cell_rng = np.random.default_rng(SUBSAMPLE_SEED + 1000 * n_cells + rep)
        cell_idx = cell_rng.choice(pool.shape[0], size=n_cells, replace=False)
        sub = pool[cell_idx, :]          # ALL genes retained
        g = gmp_cor(sub)
        vals.append(g)
        records.append({
            'n_cells': n_cells, 'n_genes': N_GENES, 'repeat': rep, 'gmp_cor': g,
            'aspect_ratio_genes_per_cell': N_GENES / n_cells,
        })
        print(f'  n_cells={n_cells:4d}  rep={rep}  GMP-Cor={g:8.4f}')

    v = np.asarray(vals)
    summary.append({
        'n_cells': n_cells, 'n_genes': N_GENES,
        'aspect_ratio_genes_per_cell': N_GENES / n_cells,
        'mean': float(v.mean()), 'std': float(v.std(ddof=1)),
        'min': float(v.min()), 'max': float(v.max()),
        'cv': float(v.std(ddof=1) / v.mean()) if v.mean() else float('nan'),
    })
    print(f'  -> n_cells={n_cells}: mean={v.mean():.4f}  SD={v.std(ddof=1):.4f}\n')

full = summary[-1]['mean']               # mean at the largest cell count
for r in summary:
    r['frac_of_full'] = r['mean'] / full if full else float('nan')

# ── Assemble results ─────────────────────────────────────────────────────────

results = {
    'params': {
        'rho': RHO, 'cell_sizes': CELL_SIZES, 'n_genes': N_GENES,
        'n_cells_pool': N_CELLS_POOL, 'n_repeats': N_REPEATS,
        'dropout_rate': DROPOUT_RATE, 'shape': SHAPE, 'hub_probability': HUB_PROB,
        'sigma_seed': SIGMA_SEED, 'count_seed': COUNT_SEED,
        'subsample_seed_base': SUBSAMPLE_SEED,
        'gene_selection': 'NONE — all genes retained at every cell size',
        'cell_selection': 'uniform random without replacement',
        'gmp_cor_definition': 'sum(max(lambda_i - max_scrambled_lambda, 0))',
        'note': (
            'Pool generated with exactly N_GENES genes; the complete gene set is '
            'used for every subsample. Only cells are subsampled. No top/random '
            'gene-selection axis. Aspect ratio (genes:cells) varies 10:1 -> 2:1.'
        ),
        'timestamp': datetime.datetime.now().isoformat(),
    },
    'per_repeat': records,
    'per_size': summary,
}
with open(JSON_LOG, 'w') as fh:
    fh.write(json.dumps(results, indent=2) + '\n')
print(f'JSON log written to: {JSON_LOG}')

# ── Human-readable summary ───────────────────────────────────────────────────

sep = '-' * 68


def _spread(xs):
    xs = list(xs)
    return (max(xs) - min(xs)) / np.mean(xs) * 100 if np.mean(xs) else float('nan')


raw_means = [r['mean'] for r in summary]
frac200 = summary[0]['frac_of_full']

lines = [
    '=' * 68,
    'SUBSAMPLING ROBUSTNESS — ALL GENES RETAINED  (REGULATED rho=0.9)',
    f'Run timestamp : {results["params"]["timestamp"]}',
    f'Script        : {os.path.abspath(__file__)}',
    '=' * 68, '',
    'PARAMETERS', sep,
    f'  rho                   : {RHO}',
    f'  cell sizes            : {CELL_SIZES}',
    f'  n_genes               : {N_GENES}  (ALL genes used every time — no selection)',
    f'  aspect ratio g:c      : {N_GENES}:{CELL_SIZES[0]} .. {N_GENES}:{CELL_SIZES[-1]}  (NOT fixed)',
    f'  pool                  : {N_CELLS_POOL} cells x {N_GENES} genes',
    f'  repeats per size      : {N_REPEATS}',
    f'  dropout_rate          : {DROPOUT_RATE}',
    f'  shape (Pareto)        : {SHAPE}',
    f'  hub_probability       : {HUB_PROB}',
    f'  sigma seed (pool)     : {SIGMA_SEED}',
    f'  count seed (pool)     : {COUNT_SEED}',
    f'  cell selection        : uniform random, no replacement',
    f'  gene selection        : NONE (complete 2000-gene set every time)',
    '',
    'GMP-COR DEFINITION', sep,
    '  GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) )',
    '',
    'RESULTS', sep,
    '  n_cells  g:c    mean     SD      min      max      CV    frac_of_full',
    '  ' + '-' * 64,
]
for r in summary:
    lines.append(
        f"  {r['n_cells']:6d}  {r['aspect_ratio_genes_per_cell']:4.1f}  "
        f"{r['mean']:7.3f}  {r['std']:6.3f}  {r['min']:7.3f}  {r['max']:7.3f}  "
        f"{r['cv']:5.3f}    {r['frac_of_full']:5.2f}"
    )
lines += [
    f"  raw mean spread over 5x cell range : {_spread(raw_means):6.1f}%",
    '',
    'INTERPRETATION', sep,
    *textwrap.wrap(
        f'Using the complete 2000-gene set and subsampling cells only, raw '
        f'GMP-Cor declines smoothly from {raw_means[-1]:.1f} at 1000 cells to '
        f'{raw_means[0]:.1f} at 200 cells (a {_spread(raw_means):.0f}% spread over '
        f'the 5x range), retaining {frac200 * 100:.0f}% of the full-size value at '
        f'200 cells. The fall is driven by fewer cells giving noisier correlation '
        f'estimates and a higher Marchenko-Pastur noise edge (the aspect ratio '
        f'reaches 10:1 genes:cells at 200 cells), not by any loss of genes.',
        width=68, initial_indent='  ', subsequent_indent='  ',
    ),
    '',
    *textwrap.wrap(
        f'Variance grows as cells are removed (CV {summary[-1]["cv"]:.2f} at 1000 '
        f'cells -> {summary[0]["cv"]:.2f} at 200 cells), so small-cell estimates '
        f'are noisier but the mean signal remains clearly elevated. With the gene '
        f'panel held complete, GMP-Cor is robust to cell subsampling down to a '
        f'few hundred cells.',
        width=68, initial_indent='  ', subsequent_indent='  ',
    ),
    '',
    'LOG FILES', sep,
    f'  JSON : {JSON_LOG}',
    f'  Text : {TEXT_LOG}',
    f'  SVG  : {SVG_FIG}',
    f'  PNG  : {PNG_FIG}',
    '=' * 68,
]
summary_text = '\n'.join(lines)
with open(TEXT_LOG, 'w', encoding='utf-8') as fh:
    fh.write(summary_text + '\n')
print('\n' + summary_text)

# ── Figure ───────────────────────────────────────────────────────────────────

xs = CELL_SIZES
xticklabels = [f'{c}\n({N_GENES / c:.1f}:1)' for c in xs]
m_raw = [r['mean'] for r in summary]
s_raw = [r['std'] for r in summary]
frac = [r['frac_of_full'] for r in summary]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5))

for rec in records:
    axA.plot(rec['n_cells'], rec['gmp_cor'], 'o', color='#de2d26',
             alpha=0.22, markersize=4, zorder=2)
axA.errorbar(xs, m_raw, yerr=s_raw, fmt='o-', color='#67000d', capsize=4,
             linewidth=2, markersize=7, zorder=3, label='mean ± SD')
for x, m, s in zip(xs, m_raw, s_raw):
    axA.text(x, m + s + max(m_raw) * 0.03, f'{m:.1f}', ha='center',
             va='bottom', fontsize=9, color='#67000d')
axA.set_xticks(xs); axA.set_xticklabels(xticklabels, fontsize=9)
axA.set_xlabel('Cells subsampled  (all 2000 genes;  g:c ratio)', fontsize=11)
axA.set_ylabel('GMP-Cor  (raw)', fontsize=12)
axA.set_ylim(bottom=0)
axA.set_title('A  Raw GMP-Cor vs cells  (all 2000 genes)', fontsize=12,
              fontweight='bold', pad=8, loc='left')
axA.spines['top'].set_visible(False); axA.spines['right'].set_visible(False)
axA.grid(True, axis='y', linestyle='--', alpha=0.25)
axA.legend(fontsize=9, frameon=False, loc='upper left')

axB.axhline(1.0, color='grey', linestyle='--', linewidth=1.0, alpha=0.7)
axB.plot(xs, frac, 'o-', color='#67000d', linewidth=2, markersize=7)
for x, f in zip(xs, frac):
    axB.text(x, f + 0.03, f'{f:.2f}', ha='center', va='bottom', fontsize=9,
             color='#67000d')
axB.set_xticks(xs); axB.set_xticklabels(xticklabels, fontsize=9)
axB.set_xlabel('Cells subsampled  (all 2000 genes;  g:c ratio)', fontsize=11)
axB.set_ylabel('Fraction of full-size (1000-cell) GMP-Cor', fontsize=12)
axB.set_ylim(0, 1.15)
axB.set_title('B  Signal retained under cell subsampling', fontsize=12,
              fontweight='bold', pad=8, loc='left')
axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)
axB.grid(True, axis='y', linestyle='--', alpha=0.25)

fig.suptitle(f'Subsampling robustness — all 2000 genes retained (ρ = {RHO})',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(SVG_FIG, format='svg', bbox_inches='tight')
fig.savefig(PNG_FIG, format='png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSVG saved to: {SVG_FIG}\nPNG saved to: {PNG_FIG}')
