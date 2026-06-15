"""
Subsampling robustness — fixed gene count (Reviewer #1, comment 2.3).

Companion to subsampling_genecmp_rho09_run.py. There the cell:gene ratio was
held fixed (genes scaled with cells). HERE the gene count is FIXED at
FIXED_N_GENES = 2000 and only the number of CELLS is subsampled
(200, 400, 600, 800, 1000). This is the more literal "profile fewer cells"
scenario: every dataset keeps the full 2000-gene panel.

Consequence: the cell:gene aspect ratio is NOT constant. It runs from
2000/200 = 10:1 (genes:cells) at 200 cells down to 2000/1000 = 2:1 at 1000
cells, so the Marchenko-Pastur noise edge shifts across sizes. Because n_genes
is constant, the per-gene index is just a constant rescale of the raw metric;
the informative robustness view is therefore the FRACTION of the full-size
(1000-cell) GMP-Cor that survives subsampling, shown in panel B.

Two gene-selection rules are compared (2000 genes drawn from the 2400-gene pool):
  - 'top'    : highest total expression first
  - 'random' : uniform random genes

Cell draws are seeded identically to the earlier runs (paired across methods).
At 1000 cells the 'top' values coincide with the ratio-fixed run (same 2000-gene
top panel, same cell draw).

GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) ).

Outputs (results/simulation_results/):
  - logs/subsampling_fixedgenes_rho09_<timestamp>.json
  - raw/subsampling_fixedgenes_rho09_<timestamp>.txt
  - figures/subsampling_fixedgenes_rho09_<timestamp>.svg / .png
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
FIXED_N_GENES  = 2000          # ← genes held constant; only cells subsampled
GENE_SIZES     = [FIXED_N_GENES] * len(CELL_SIZES)

N_CELLS_POOL   = 1200
N_GENES_POOL   = 2400

N_REPEATS      = 5
DROPOUT_RATE   = 1.0
SHAPE          = 1.5
HUB_PROB       = 0.2
SIGMA_SEED     = 31
COUNT_SEED     = 0
SUBSAMPLE_SEED = 1000          # base for CELL draws (matches earlier runs)
GENE_RNG_BASE  = 2_000_000     # base for RANDOM gene draws

METHODS = ['top', 'random']
METHOD_LABEL = {'top': 'highest expression first', 'random': 'uniform random genes'}

# ── Output paths ─────────────────────────────────────────────────────────────

_SIM_RESULTS = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
_FIG_DIR = os.path.join(_SIM_RESULTS, 'figures')
_RAW_DIR = os.path.join(_SIM_RESULTS, 'raw')
_LOG_DIR = os.path.join(_SIM_RESULTS, 'logs')
for _d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
    os.makedirs(_d, exist_ok=True)

_ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
JSON_LOG = os.path.join(_LOG_DIR, f'subsampling_fixedgenes_rho09_{_ts}.json')
TEXT_LOG = os.path.join(_RAW_DIR, f'subsampling_fixedgenes_rho09_{_ts}.txt')
SVG_FIG  = os.path.join(_FIG_DIR, f'subsampling_fixedgenes_rho09_{_ts}.svg')
PNG_FIG  = os.path.join(_FIG_DIR, f'subsampling_fixedgenes_rho09_{_ts}.png')


def gmp_cor(observed):
    pcs, pcs1, _ = get_eig_dist(observed, norm=True, log=False, norm_sum=100)
    return float(np.sum(np.maximum(pcs - pcs1.max(), 0)))


# ── Generate the pool (same seeds → identical pool to earlier runs) ───────────

print(f'Generating pool: {N_CELLS_POOL} cells x {N_GENES_POOL} genes at rho={RHO} ...')
sigma_pool = generate_gram_hub_matrix(N_GENES_POOL, RHO, SHAPE, HUB_PROB, seed=SIGMA_SEED)
_, pool = simulate_scRNA_data(
    n_cells=N_CELLS_POOL, n_genes=N_GENES_POOL, sigma=sigma_pool,
    dropout_rate=DROPOUT_RATE, seed=COUNT_SEED,
)
print(f'Pool generated: shape={pool.shape}  (genes fixed at {FIXED_N_GENES})\n')

# ── Subsampling sweep (cells only; genes fixed) ──────────────────────────────

records = []
summary = {m: [] for m in METHODS}

for n_cells in CELL_SIZES:
    n_genes = FIXED_N_GENES
    vals = {m: [] for m in METHODS}
    for rep in range(N_REPEATS):
        cell_rng = np.random.default_rng(SUBSAMPLE_SEED + 1000 * n_cells + rep)
        cell_idx = cell_rng.choice(pool.shape[0], size=n_cells, replace=False)
        sub_cells = pool[cell_idx, :]
        gene_totals = sub_cells.sum(axis=0)

        for method in METHODS:
            if method == 'top':
                gene_idx = np.argsort(gene_totals)[::-1][:n_genes]
            else:  # random
                gene_rng = np.random.default_rng(GENE_RNG_BASE + 1000 * n_cells + rep)
                gene_idx = gene_rng.choice(pool.shape[1], size=n_genes, replace=False)

            sub = sub_cells[:, gene_idx]
            g = gmp_cor(sub)
            records.append({
                'method': method, 'n_cells': n_cells, 'n_genes': n_genes,
                'repeat': rep, 'gmp_cor': g, 'gmp_cor_per_gene': g / n_genes,
                'aspect_ratio_genes_per_cell': n_genes / n_cells,
            })
            vals[method].append(g)
            print(f'  [{method:6s}] n_cells={n_cells:4d} n_genes={n_genes:4d} '
                  f'rep={rep} GMP-Cor={g:8.4f}')

    for method in METHODS:
        v = np.asarray(vals[method])
        summary[method].append({
            'n_cells': n_cells, 'n_genes': n_genes,
            'aspect_ratio_genes_per_cell': n_genes / n_cells,
            'mean': float(v.mean()), 'std': float(v.std(ddof=1)),
            'min': float(v.min()), 'max': float(v.max()),
            'cv': float(v.std(ddof=1) / v.mean()) if v.mean() else float('nan'),
        })
    print()

# Fraction of full-size (largest cell count) GMP-Cor retained, per method.
for method in METHODS:
    full = summary[method][-1]['mean']  # mean at the largest cell count
    for r in summary[method]:
        r['frac_of_full'] = r['mean'] / full if full else float('nan')

# ── Assemble results ─────────────────────────────────────────────────────────

results = {
    'params': {
        'rho': RHO, 'cell_sizes': CELL_SIZES, 'fixed_n_genes': FIXED_N_GENES,
        'n_cells_pool': N_CELLS_POOL, 'n_genes_pool': N_GENES_POOL,
        'n_repeats': N_REPEATS, 'dropout_rate': DROPOUT_RATE, 'shape': SHAPE,
        'hub_probability': HUB_PROB, 'sigma_seed': SIGMA_SEED,
        'count_seed': COUNT_SEED, 'subsample_seed_base': SUBSAMPLE_SEED,
        'gene_rng_base': GENE_RNG_BASE, 'methods': METHODS,
        'cell_selection': 'uniform random without replacement (paired across methods)',
        'gene_selection': METHOD_LABEL,
        'gmp_cor_definition': 'sum(max(lambda_i - max_scrambled_lambda, 0))',
        'note': (
            'Gene count fixed at FIXED_N_GENES; only cells subsampled. The '
            'cell:gene aspect ratio is therefore NOT held fixed (genes:cells runs '
            '10:1 -> 2:1). Companion to the ratio-fixed run '
            'subsampling_genecmp_rho09_run.py.'
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

sep = '-' * 72


def _spread(xs):
    xs = list(xs)
    return (max(xs) - min(xs)) / np.mean(xs) * 100 if np.mean(xs) else float('nan')


lines = [
    '=' * 72,
    'SUBSAMPLING ROBUSTNESS — FIXED GENE COUNT  (REGULATED rho=0.9)',
    f'Run timestamp : {results["params"]["timestamp"]}',
    f'Script        : {os.path.abspath(__file__)}',
    '=' * 72, '',
    'PARAMETERS', sep,
    f'  rho                   : {RHO}',
    f'  cell sizes            : {CELL_SIZES}',
    f'  n_genes (FIXED)       : {FIXED_N_GENES}  (selected from {N_GENES_POOL}-gene pool)',
    f'  aspect ratio g:c      : {FIXED_N_GENES}:{CELL_SIZES[0]} .. {FIXED_N_GENES}:{CELL_SIZES[-1]}  (NOT fixed)',
    f'  pool                  : {N_CELLS_POOL} cells x {N_GENES_POOL} genes',
    f'  repeats per size      : {N_REPEATS}',
    f'  dropout_rate          : {DROPOUT_RATE}',
    f'  shape (Pareto)        : {SHAPE}',
    f'  hub_probability       : {HUB_PROB}',
    f'  sigma seed (pool)     : {SIGMA_SEED}',
    f'  count seed (pool)     : {COUNT_SEED}',
    f'  cell selection        : uniform random, no replacement (paired)',
    f'  gene selection        : top = highest expression; random = uniform',
    '',
    'GMP-COR DEFINITION', sep,
    '  GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) )',
    '',
    'RESULTS', sep,
]

for method in METHODS:
    rows = summary[method]
    lines += [
        '',
        f'  GENE SELECTION: {method.upper()}  ({METHOD_LABEL[method]})',
        '    n_cells  g:c    mean     SD      min      max      CV    frac_of_full',
        '    ' + '-' * 64,
    ]
    for r in rows:
        lines.append(
            f"    {r['n_cells']:6d}  {r['aspect_ratio_genes_per_cell']:4.1f}  "
            f"{r['mean']:7.3f}  {r['std']:6.3f}  {r['min']:7.3f}  {r['max']:7.3f}  "
            f"{r['cv']:5.3f}    {r['frac_of_full']:5.2f}"
        )
    raw_means = [r['mean'] for r in rows]
    lines.append(f"    raw mean spread over 5x cell range : {_spread(raw_means):6.1f}%")

# Side-by-side comparison
lines += ['', 'GMP-COR — top vs random at each cell count (genes fixed = %d)' % FIXED_N_GENES,
          sep, '    n_cells    top      random    random/top']
for rt, rr in zip(summary['top'], summary['random']):
    ratio = rr['mean'] / rt['mean'] if rt['mean'] else float('nan')
    lines.append(f"    {rt['n_cells']:6d}   {rt['mean']:7.3f}  {rr['mean']:7.3f}    {ratio:5.2f}x")

top_m = [r['mean'] for r in summary['top']]
rnd_m = [r['mean'] for r in summary['random']]
top_frac = summary['top'][0]['frac_of_full']
rnd_frac = summary['random'][0]['frac_of_full']
lines += [
    '',
    'INTERPRETATION', sep,
    *textwrap.wrap(
        f'With the full {FIXED_N_GENES}-gene panel retained and only cells '
        f'subsampled, raw GMP-Cor still falls as cells are removed under both '
        f'gene-selection rules (top: {top_m[0]:.1f} -> {top_m[-1]:.1f}; random: '
        f'{rnd_m[0]:.1f} -> {rnd_m[-1]:.1f} from 200 to 1000 cells). Fewer cells '
        f'means noisier correlation estimates and a higher Marchenko-Pastur noise '
        f'edge (aspect ratio rises to 10:1 genes:cells at 200 cells), so less '
        f'eigenvalue mass clears the threshold.',
        width=72, initial_indent='  ', subsequent_indent='  ',
    ),
    '',
    *textwrap.wrap(
        f'At 200 cells the top-expressed panel still retains {top_frac * 100:.0f}% '
        f'of its full-size (1000-cell) GMP-Cor, whereas the random panel retains '
        f'only {rnd_frac * 100:.0f}%. Keeping all 2000 genes does not rescue the '
        f'random selection at small cell numbers: informative (highly-expressed) '
        f'genes are what carry the detectable correlation structure. The signal '
        f'is recovered far more stably across cell numbers when genes are chosen '
        f'by expression — consistent with restricting analysis to detected / '
        f'variable genes as done for the experimental data.',
        width=72, initial_indent='  ', subsequent_indent='  ',
    ),
    '',
    'LOG FILES', sep,
    f'  JSON : {JSON_LOG}',
    f'  Text : {TEXT_LOG}',
    f'  SVG  : {SVG_FIG}',
    f'  PNG  : {PNG_FIG}',
    '=' * 72,
]
summary_text = '\n'.join(lines)
with open(TEXT_LOG, 'w', encoding='utf-8') as fh:
    fh.write(summary_text + '\n')
print('\n' + summary_text)

# ── Figure ───────────────────────────────────────────────────────────────────

STYLE = {
    'top':    dict(color='#67000d', marker='o', label='top expression'),
    'random': dict(color='#08519c', marker='s', label='random genes'),
}
PT = {'top': '#de2d26', 'random': '#3182bd'}

xs = CELL_SIZES
xticklabels = [f'{c}\n({FIXED_N_GENES / c:.1f}:1)' for c in xs]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5))

for method in METHODS:
    rows = summary[method]
    m_raw = [r['mean'] for r in rows]
    s_raw = [r['std'] for r in rows]
    frac = [r['frac_of_full'] for r in rows]
    st = STYLE[method]

    for rec in records:
        if rec['method'] != method:
            continue
        axA.plot(rec['n_cells'], rec['gmp_cor'], st['marker'], color=PT[method],
                 alpha=0.2, markersize=4, zorder=2)

    axA.errorbar(xs, m_raw, yerr=s_raw, fmt=st['marker'] + '-', color=st['color'],
                 capsize=4, linewidth=2, markersize=7, zorder=3, label=st['label'])
    axB.plot(xs, frac, st['marker'] + '-', color=st['color'], linewidth=2,
             markersize=7, label=st['label'])

axA.set_xticks(xs); axA.set_xticklabels(xticklabels, fontsize=9)
axA.set_xlabel('Cells subsampled  (genes fixed = 2000;  g:c ratio)', fontsize=11)
axA.set_ylabel('GMP-Cor  (raw)', fontsize=12)
axA.set_ylim(bottom=0)
axA.set_title('A  Raw GMP-Cor vs cells  (2000 genes fixed)', fontsize=12,
              fontweight='bold', pad=8, loc='left')
axA.spines['top'].set_visible(False); axA.spines['right'].set_visible(False)
axA.grid(True, axis='y', linestyle='--', alpha=0.25)
axA.legend(fontsize=9, frameon=False, loc='upper left')

axB.axhline(1.0, color='grey', linestyle='--', linewidth=1.0, alpha=0.7)
axB.set_xticks(xs); axB.set_xticklabels(xticklabels, fontsize=9)
axB.set_xlabel('Cells subsampled  (genes fixed = 2000;  g:c ratio)', fontsize=11)
axB.set_ylabel('Fraction of full-size (1000-cell) GMP-Cor', fontsize=12)
axB.set_ylim(0, 1.15)
axB.set_title('B  Signal retained under cell subsampling', fontsize=12,
              fontweight='bold', pad=8, loc='left')
axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)
axB.grid(True, axis='y', linestyle='--', alpha=0.25)
axB.legend(fontsize=9, frameon=False, loc='lower right')

fig.suptitle(f'Subsampling robustness — fixed 2000-gene panel (ρ = {RHO})',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(SVG_FIG, format='svg', bbox_inches='tight')
fig.savefig(PNG_FIG, format='png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSVG saved to: {SVG_FIG}\nPNG saved to: {PNG_FIG}')
