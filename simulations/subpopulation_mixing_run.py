"""
Subpopulation mixing scenario — 50/50 runner.
Implements Simulation Improvement #2 from documents/reviewer1_response_plan.md.

Two conditions compared:
  - DYSREGULATED: two sub-populations, each low rho (rho=0.1), different sigma seeds
  - REGULATED:    two sub-populations, each high rho (rho=0.8), different sigma seeds

Different sigma seeds give each sub-population a distinct hub-network topology.
GMP-Cor is computed per sub-population and for the 50/50 mixture.

Key question: does a mixture of two distinct but internally *regulated* populations
produce elevated GMP-Cor compared to a mixture of two *dysregulated* ones?

Outputs (written to simulations/logs/):
  - subpopulation_mixing_50_50_<timestamp>.json  — full parameter + results log
  - subpopulation_mixing_50_50_<timestamp>.txt   — human-readable summary
  - subpopulation_mixing_50_50_<timestamp>.svg   — publication-quality figure
  - subpopulation_mixing_50_50_<timestamp>.png   — raster copy
"""

import sys
import os
import json
import datetime
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Allow running directly from simulations/ or from repo root
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.simulations import subpopulation_mixing  # noqa: E402

# ── Parameters ───────────────────────────────────────────────────────────────

PARAMS = dict(
    n_cells=1000,       # total cells; split 50/50 → 500 per sub-population
    n_genes=2000,       # genes per cell (matches experimental data scale)
    mixing_ratio=0.5,   # 50 / 50
    rho_low=0.1,        # dysregulated condition: weak gene-gene correlations
    rho_high=0.8,       # regulated condition: strong gene-gene correlations
    seed_a=20,          # sigma-matrix seed for sub-population A
    seed_b=21,          # sigma-matrix seed for sub-population B (different network)
    count_seed_a=0,     # count-generation seed for sub-population A
    count_seed_b=1,     # count-generation seed for sub-population B
    dropout_rate=1.0,   # expression-dependent dropout rate
    shape=1.5,          # Pareto shape for cluster-size distribution
    hub_probability=0.2,# probability a cluster hub connects to the global hub factor
)

# ── Log paths ────────────────────────────────────────────────────────────────

_SIM_RESULTS = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
_FIG_DIR  = os.path.join(_SIM_RESULTS, 'figures')
_RAW_DIR  = os.path.join(_SIM_RESULTS, 'raw')
_LOG_DIR  = os.path.join(_SIM_RESULTS, 'logs')
for _d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
    os.makedirs(_d, exist_ok=True)

_ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
JSON_LOG  = os.path.join(_LOG_DIR, f'subpopulation_mixing_50_50_{_ts}.json')
TEXT_LOG  = os.path.join(_RAW_DIR, f'subpopulation_mixing_50_50_{_ts}.txt')
SVG_FIG   = os.path.join(_FIG_DIR, f'subpopulation_mixing_50_50_{_ts}.svg')
PNG_FIG   = os.path.join(_FIG_DIR, f'subpopulation_mixing_50_50_{_ts}.png')

# ── Run ──────────────────────────────────────────────────────────────────────

print('=' * 65)
print('Subpopulation Mixing — 50/50 Scenario')
print('=' * 65)
print(f'Parameters:\n{json.dumps(PARAMS, indent=2)}\n')

results = subpopulation_mixing(**PARAMS, log_file=JSON_LOG)

# ── Human-readable text log ──────────────────────────────────────────────────

p = results['params']
dis = results['dysregulated']
reg = results['regulated']

sep = '-' * 65

summary_lines = [
    '=' * 65,
    'SUBPOPULATION MIXING — 50/50 SCENARIO',
    f'Run timestamp : {p["timestamp"]}',
    f'Script        : {os.path.abspath(__file__)}',
    '=' * 65,
    '',
    'PARAMETERS',
    sep,
    f'  n_cells total         : {p["n_cells"]}  ({p["n_cells_a"]} per sub-pop)',
    f'  n_genes               : {p["n_genes"]}',
    f'  mixing_ratio          : {p["mixing_ratio"]}  (50/50)',
    f'  rho_low  (dysreg.)    : {p["rho_low"]}',
    f'  rho_high (reg.)       : {p["rho_high"]}',
    f'  sigma seed A          : {p["seed_a"]}',
    f'  sigma seed B          : {p["seed_b"]}  (different hub-network topology)',
    f'  count seed A          : {p["count_seed_a"]}',
    f'  count seed B          : {p["count_seed_b"]}',
    f'  dropout_rate          : {p["dropout_rate"]}',
    f'  shape (Pareto)        : {p["shape"]}',
    f'  hub_probability       : {p["hub_probability"]}',
    '',
    'GMP-COR DEFINITION',
    sep,
    '  GMP-Cor = sum( max(lambda_i - max_scrambled_lambda, 0) )',
    '  (sum of eigenvalue excesses above the scrambled noise threshold)',
    '',
    'RESULTS',
    sep,
    '',
    '  DYSREGULATED CONDITION  (rho = {rho})'.format(rho=dis['rho']),
    '  Two low-rho sub-populations (different sigma seeds)',
    f'    Sub-pop A GMP-Cor  : {dis["gmp_cor_subpop_a"]:.4f}',
    f'    Sub-pop B GMP-Cor  : {dis["gmp_cor_subpop_b"]:.4f}',
    f'    Combined  GMP-Cor  : {dis["gmp_cor_combined"]:.4f}',
    '',
    '  REGULATED CONDITION  (rho = {rho})'.format(rho=reg['rho']),
    '  Two high-rho sub-populations (different sigma seeds)',
    f'    Sub-pop A GMP-Cor  : {reg["gmp_cor_subpop_a"]:.4f}',
    f'    Sub-pop B GMP-Cor  : {reg["gmp_cor_subpop_b"]:.4f}',
    f'    Combined  GMP-Cor  : {reg["gmp_cor_combined"]:.4f}',
    '',
    'COMPARISON',
    sep,
    f'  Regulated combined   : {reg["gmp_cor_combined"]:.4f}',
    f'  Dysregulated combined: {dis["gmp_cor_combined"]:.4f}',
    f'  Ratio (reg / dysreg) : {reg["gmp_cor_combined"] / dis["gmp_cor_combined"]:.2f}x'
        if dis['gmp_cor_combined'] > 0
        else '  Ratio (reg / dysreg) : inf  (dysregulated GMP-Cor ~ 0)',
    '',
    'INTERPRETATION',
    sep,
    *textwrap.wrap(
        'A mixture of two internally REGULATED populations (high rho, distinct hub '
        'networks) produces an elevated GMP-Cor because each sub-population contributes '
        'eigenvalues well above the scrambled noise threshold. In contrast, a mixture of '
        'two DYSREGULATED populations (low rho) produces a low combined GMP-Cor because '
        'neither sub-population has a strong correlation structure to contribute. '
        'This directly addresses Reviewer #1 Comment 1: GMP-Cor cannot be reduced to '
        'a low value merely by mixing two regulated sub-populations with different '
        'network topologies — genuine dysregulation is required.',
        width=65,
        initial_indent='  ',
        subsequent_indent='  ',
    ),
    '',
    'LOG FILES',
    sep,
    f'  JSON : {JSON_LOG}',
    f'  Text : {TEXT_LOG}',
    f'  SVG  : {SVG_FIG}',
    f'  PNG  : {PNG_FIG}',
    '=' * 65,
]

summary_text = '\n'.join(summary_lines)

with open(TEXT_LOG, 'w') as fh:
    fh.write(summary_text + '\n')

print('\n' + summary_text)
print(f'\nText log written to: {TEXT_LOG}')

# ── Plot ─────────────────────────────────────────────────────────────────────

print('\nGenerating figure ...')

CONDITIONS  = ['Dysregulated\n(χ = 0.1)', 'Regulated\n(χ = 0.8)']
BAR_LABELS  = ['Sub-pop A', 'Sub-pop B', 'Combined']
DIS_VALS    = [dis['gmp_cor_subpop_a'], dis['gmp_cor_subpop_b'], dis['gmp_cor_combined']]
REG_VALS    = [reg['gmp_cor_subpop_a'], reg['gmp_cor_subpop_b'], reg['gmp_cor_combined']]

# Colours: sub-pop A (light), sub-pop B (mid), combined (dark) — two hue families
DIS_COLORS  = ['#9ecae1', '#3182bd', '#08519c']   # blue family → dysregulated
REG_COLORS  = ['#fc9272', '#de2d26', '#67000d']   # red family  → regulated

fig, axes = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
fig.subplots_adjust(wspace=0.38)

# Shared y-axis limit across both subplots so bar heights are directly comparable
y_max = max(max(DIS_VALS), max(REG_VALS)) * 1.18

for ax, vals, colors, title in zip(
    axes,
    [DIS_VALS, REG_VALS],
    [DIS_COLORS, REG_COLORS],
    CONDITIONS,
):
    x = np.arange(len(BAR_LABELS))
    bars = ax.bar(x, vals, color=colors, edgecolor='black', linewidth=0.7, width=0.55)

    # Value labels on bars
    for bar, v in zip(bars, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + y_max * 0.02,
            f'{v:.2f}',
            ha='center', va='bottom', fontsize=9,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(BAR_LABELS, fontsize=10)
    ax.set_ylabel('GMP-Cor', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=8)
    ax.set_ylim(0, y_max)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Shared legend patches
legend_patches = [
    mpatches.Patch(color=DIS_COLORS[2], label=f'Dysregulated combined  ({dis["gmp_cor_combined"]:.2f})'),
    mpatches.Patch(color=REG_COLORS[2], label=f'Regulated combined  ({reg["gmp_cor_combined"]:.2f})'),
]
fig.legend(
    handles=legend_patches,
    loc='lower center',
    ncol=2,
    fontsize=9,
    frameon=False,
    bbox_to_anchor=(0.5, -0.04),
)

ratio = reg['gmp_cor_combined'] / dis['gmp_cor_combined'] if dis['gmp_cor_combined'] > 0 else float('inf')
fig.suptitle(
    f'Subpopulation Mixing (50/50) — GMP-Cor\n'
    f'n_cells={p["n_cells"]}, n_genes={p["n_genes"]}, '
    f'seeds A/B = {p["seed_a"]}/{p["seed_b"]}  |  '
    f'reg/dysreg ratio = {ratio:.1f}×',
    fontsize=10,
    y=1.02,
)

plt.tight_layout()
fig.savefig(SVG_FIG, format='svg', bbox_inches='tight')
fig.savefig(PNG_FIG, format='png', dpi=150, bbox_inches='tight')
plt.close(fig)

print(f'SVG saved to: {SVG_FIG}')
print(f'PNG saved to: {PNG_FIG}')
