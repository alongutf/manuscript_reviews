"""
Eigenvector block-recovery test — headline condition (alpha = 0.75).

Question: can the eigenvector analysis used on the experimental data
(`scripts/eigenvector_analysis.py`) reconstruct the gene-gene correlation blocks
when we KNOW they are there?

The synthetic generator builds its correlation matrix from contiguous gene clusters
plus a global hub factor, so the block membership of every gene is ground truth
(`generate_gram_hub_matrix(..., return_structure=True)`). We take the raw loadings
of the TOP 5 MODES only -- exactly what the experimental script inspects, no
coordination score -- and score them against those blocks.

Four stages of the same simulation are scored, which localizes which noise layer
(if any) destroys recovery:

  R_oracle    eigenvectors of the true correlation matrix      (5-mode ceiling)
  latent      the correlated Gaussian layer, no counts         (sampling noise only)
  true_counts negative-binomial counts + library-size scaling  (count noise)
  observed    after expression-dependent dropout               (what real data looks like)

Outputs (results/simulation_results/):
  logs/eigenvector_block_recovery_<ts>.json  full parameters + all per-mode/per-block rows
  raw/eigenvector_block_recovery_<ts>.txt    human-readable summary
  figures/eigenvector_block_recovery_<ts>.svg / .png
"""

import sys
import os
import json
import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.simulations import generate_gram_hub_matrix, simulate_scRNA_data  # noqa: E402
import src.eigenvector_validation as ev  # noqa: E402

# ── Parameters ───────────────────────────────────────────────────────────────

PARAMS = dict(
    n_genes=2000,          # matches the experimental panel scale
    n_cells=1000,
    alpha=0.75,            # coupling strength (shared-variance fraction)
    shape=1.5,             # Pareto shape for cluster sizes
    hub_probability=0.2,   # probability a cluster hub joins the global hub factor
    sigma_seed=31,         # correlation-matrix seed (fixes the ground-truth blocks)
    count_seed=0,          # count-sampling seed for the scored replicate
    count_seed_rep=1,      # second count seed, used only for the stability check
    dropout_rate=1.0,      # P(dropout) = exp(-rate * count)
    inv_gamma_shape=1.5,   # per-gene NB mean prior
    inv_gamma_scale=0.03,   # calibrated so the observed zero fraction (~96.8%)
                            # matches the real-data median (96.6%, range 87.5-99.0)
    n_top=5,               # modes inspected -- same as scripts/eigenvector_analysis.py
    n_top_genes=30,        # top-loading genes per mode (N_GENES in that script)
    min_block_size=5,      # blocks smaller than this cannot be recovered in principle
)

STAGES = ['R_oracle', 'latent', 'true_counts', 'observed']
STAGE_LABEL = {
    'R_oracle': 'true R\n(5-mode ceiling)',
    'latent': 'latent Gaussian\n(sampling only)',
    'true_counts': 'NB counts\n(+ library size)',
    'observed': 'observed\n(+ dropout)',
}

# ── Output paths ─────────────────────────────────────────────────────────────

_SIM_RESULTS = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
_FIG_DIR = os.path.join(_SIM_RESULTS, 'figures')
_RAW_DIR = os.path.join(_SIM_RESULTS, 'raw')
_LOG_DIR = os.path.join(_SIM_RESULTS, 'logs')
for _d in (_FIG_DIR, _RAW_DIR, _LOG_DIR):
    os.makedirs(_d, exist_ok=True)

_ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
_stem = f'eigenvector_block_recovery_{_ts}'
JSON_LOG = os.path.join(_LOG_DIR, f'{_stem}.json')
TEXT_LOG = os.path.join(_RAW_DIR, f'{_stem}.txt')
SVG_FIG = os.path.join(_FIG_DIR, f'{_stem}.svg')
PNG_FIG = os.path.join(_FIG_DIR, f'{_stem}.png')

# ── Simulate ─────────────────────────────────────────────────────────────────

print('=' * 70)
print('Eigenvector block recovery — headline condition')
print('=' * 70)
print(f'Parameters:\n{json.dumps(PARAMS, indent=2)}\n')

p = PARAMS
R, structure = generate_gram_hub_matrix(
    p['n_genes'], p['alpha'], p['shape'], p['hub_probability'],
    seed=p['sigma_seed'], return_structure=True)

sizes, strengths = structure['sizes'], structure['strengths']
print(f'Ground truth: {len(sizes)} blocks, {(sizes == 1).sum()} singletons, '
      f'{(sizes >= p["min_block_size"]).sum()} of size >= {p["min_block_size"]}, '
      f'largest = {sizes.max()} genes')

# latent Gaussian layer, drawn with the same seed the count simulation uses so the
# ladder stages describe one coherent dataset
np.random.seed(p['count_seed'])
latent = np.random.multivariate_normal(np.zeros(p['n_genes']), R, size=p['n_cells'])

print('Simulating counts (this is the slow step)...')
true_counts, observed_counts = simulate_scRNA_data(
    n_cells=p['n_cells'], n_genes=p['n_genes'], sigma=R,
    dropout_rate=p['dropout_rate'], seed=p['count_seed'],
    inv_gamma_shape=p['inv_gamma_shape'], inv_gamma_scale=p['inv_gamma_scale'])

# ── Score every ladder stage ─────────────────────────────────────────────────

stages = ev.ladder_stages(R, latent, true_counts, observed_counts, n_top=p['n_top'])
results = {}
for name in STAGES:
    eigvals, eigvecs, threshold, kept_cols = stages[name]
    scored = ev.score_stage(eigvecs, kept_cols, structure,
                            n_top_genes=p['n_top_genes'],
                            min_block_size=p['min_block_size'])
    scored['eigenvalues'] = [float(x) for x in eigvals]
    scored['threshold_scrambled'] = float(threshold)
    scored['n_modes_above_threshold'] = int(np.sum(np.asarray(eigvals) > threshold)) \
        if np.isfinite(threshold) else None
    results[name] = scored
    print(f'  scored {name}: {scored["n_genes_kept"]} genes kept, '
          f'{len(scored["blocks"])} scorable blocks')

# ── GMP-Cor and sparsity of the same matrices ────────────────────────────────
# The scalar the paper reports, computed on the very matrices scored above, so the
# eigenvector read-out and GMP-Cor can be compared on one dataset.

from src.simulations import gmp_cor  # noqa: E402


def _sparsity(m):
    """Zero fraction raw and after the all-zero gene/cell filter get_eig_dist applies."""
    m = np.asarray(m, dtype=float)
    f = m[:, (m > 0).sum(axis=0) >= 1]
    f = f[(f > 0).sum(axis=1) >= 1, :]
    return dict(zero_fraction_raw=float(1 - (m > 0).mean()),
                zero_fraction_filtered=float(1 - (f > 0).mean()),
                median_counts_per_cell=float(np.median(m.sum(axis=1))),
                median_genes_per_cell=float(np.median((m > 0).sum(axis=1))))


spectrum = {}
for name, m in (('true_counts', true_counts), ('observed', observed_counts)):
    spectrum[name] = dict(gmp_cor(m, norm=True, norm_sum=50), **_sparsity(m))
    r = spectrum[name]
    print(f'  {name}: GMP-Cor={r["gmp_cor"]:.1f} (per gene {r["gmp_cor_per_gene"]:.4f}), '
          f'lambda_max={r["lambda_max"]:.2f}, zeros={r["zero_fraction_filtered"]:.4f}')
# the true correlation matrix has no cells to scramble, so GMP-Cor is undefined there;
# its leading eigenvalues are the noiseless ceiling
spectrum['R_oracle'] = dict(lambda_max=float(results['R_oracle']['eigenvalues'][0]),
                            gmp_cor=None)

# ── Stability: same R, independent count-sampling replicate ──────────────────

print('Simulating replicate for the stability check...')
_, observed_rep = simulate_scRNA_data(
    n_cells=p['n_cells'], n_genes=p['n_genes'], sigma=R,
    dropout_rate=p['dropout_rate'], seed=p['count_seed_rep'],
    inv_gamma_shape=p['inv_gamma_shape'], inv_gamma_scale=p['inv_gamma_scale'])

from src.analysis_functions import get_eig_vectors  # noqa: E402
vecs_a, kept_a = stages['observed'][1], stages['observed'][3]
_, vecs_b, _, kept_b = get_eig_vectors(observed_rep, n_top=p['n_top'],
                                       norm_method='sum', norm_sum=50)
common = np.asarray(kept_a, dtype=bool) & np.asarray(kept_b, dtype=bool)
# restrict both replicates to the genes kept in both, then re-normalize the rows
idx_a = common[np.asarray(kept_a, dtype=bool)]
idx_b = common[np.asarray(kept_b, dtype=bool)]
va, vb = vecs_a[:, idx_a], vecs_b[:, idx_b]
va = va / np.linalg.norm(va, axis=1, keepdims=True)
vb = vb / np.linalg.norm(vb, axis=1, keepdims=True)
stability = dict(
    n_genes_common=int(common.sum()),
    per_mode_abs_r=ev.mode_agreement(va, vb),
    subspace_overlap=ev.subspace_overlap(va, vb),
)
print(f'  stability: subspace overlap = {stability["subspace_overlap"]:.3f}, '
      f'per-mode |r| = {[round(x, 3) for x in stability["per_mode_abs_r"]]}')

# ── JSON log ─────────────────────────────────────────────────────────────────

log = dict(
    experiment='eigenvector_block_recovery',
    timestamp=_ts,
    description=('Can the top-5 eigenvector loadings reconstruct the ground-truth '
                 'correlation blocks of the synthetic generator?'),
    params=p,
    ground_truth=dict(
        n_blocks=int(len(sizes)),
        n_singletons=int((sizes == 1).sum()),
        n_blocks_scorable=int((sizes >= p['min_block_size']).sum()),
        largest_block=int(sizes.max()),
        block_sizes=[int(x) for x in sizes],
        block_strengths=[float(x) for x in strengths],
        n_hubs_connected=int((structure['hub_weights'] != 0).sum()),
    ),
    stages=results,
    spectrum=spectrum,
    stability=stability,
)
with open(JSON_LOG, 'w', encoding='utf-8') as fh:
    json.dump(log, fh, indent=2)

# ── Text summary ─────────────────────────────────────────────────────────────

lines = [
    '=' * 78,
    'EIGENVECTOR BLOCK RECOVERY — headline condition',
    '=' * 78,
    f'timestamp        : {_ts}',
    f'coupling alpha   : {p["alpha"]}',
    f'matrix           : {p["n_cells"]} cells x {p["n_genes"]} genes, dropout '
    f'rate {p["dropout_rate"]}',
    f'modes inspected  : top {p["n_top"]} (raw loadings, no coordination score)',
    '',
    'GROUND TRUTH',
    f'  blocks                 : {len(sizes)} ({(sizes == 1).sum()} singletons)',
    f'  blocks >= {p["min_block_size"]} genes      : {(sizes >= p["min_block_size"]).sum()}',
    f'  largest block          : {sizes.max()} genes '
    f'(strength {strengths[np.argmax(sizes)]:.2f})',
    f'  hubs on global factor  : {(structure["hub_weights"] != 0).sum()}',
    '',
]

for name in STAGES:
    r = results[name]
    lines.append('-' * 78)
    thr = r['threshold_scrambled']
    thr_s = 'n/a' if not np.isfinite(thr) else f'{thr:.3f}'
    lines.append(f'STAGE: {name}   ({r["n_genes_kept"]} genes kept, '
                 f'lambda_max^scr = {thr_s}, '
                 f'modes above threshold = {r["n_modes_above_threshold"]})')
    lines.append('  per mode — is the mode a single block?')
    lines.append('    mode  eig     best block (size/strength)   mass   [null]   '
                 'eff.blocks  top30 purity  p_hyper')
    for m in r['modes']:
        lines.append(
            f'    {m["mode"]:>4}  {r["eigenvalues"][m["mode"]-1]:>6.2f}  '
            f'#{m["best_block"]:<5} ({m["best_block_size_kept"]:>4}g / '
            f'{m["best_block_strength"]:.2f})   {m["best_block_mass"]:.3f}  '
            f'[{m["best_block_mass_null"]:.4f}]   '
            f'{m["effective_n_blocks"]:>8.1f}   {m["top_gene_purity"]:>8.2f}     '
            f'{m["top_gene_p_hypergeom"]:.1e}')
    blocks = r['blocks']
    if blocks:
        proj = np.array([b['projection'] for b in blocks])
        ratio = np.array([b['mass_ratio'] for b in blocks])
        bsz = np.array([b['size_kept'] for b in blocks])
        order = np.argsort(bsz)[::-1][:8]
        lines.append(f'  per block — {len(blocks)} blocks of >= {p["min_block_size"]} genes; '
                     f'projection onto the {p["n_top"]}-mode subspace: '
                     f'median {np.median(proj):.3f}, max {proj.max():.3f}, '
                     f'{np.sum(proj > 0.5)} of {len(proj)} above 0.5')
        lines.append('    largest blocks:  size  strength  projection  mass/null  best mode')
        for i in order:
            b = blocks[i]
            lines.append(f'                    {b["size_kept"]:>5}  {b["strength"]:>8.2f}  '
                         f'{b["projection"]:>10.3f}  {b["mass_ratio"]:>9.1f}  '
                         f'{b["best_mode"]:>9}')
        lines.append(f'    mass/null over all scorable blocks: median {np.median(ratio):.1f}')
    lines.append('')

lines += [
    '-' * 78,
    'GMP-Cor AND SPARSITY (same matrices as scored above)',
]
for name in ('true_counts', 'observed'):
    r = spectrum[name]
    lines += [
        f'  {name}',
        f'    GMP-Cor            : {r["gmp_cor"]:.1f}  '
        f'(per gene {r["gmp_cor_per_gene"]:.4f}, {r["n_modes_above_threshold"]} modes above)',
        f'    lambda_max / scr   : {r["lambda_max"]:.2f} / {r["lambda_max_scrambled"]:.2f}',
        f'    zero fraction      : {r["zero_fraction_raw"]:.4f} raw, '
        f'{r["zero_fraction_filtered"]:.4f} after filtering',
        f'    median per cell    : {r["median_counts_per_cell"]:.0f} counts, '
        f'{r["median_genes_per_cell"]:.0f} genes',
    ]
lines += [
    f'  true R lambda_max     : {spectrum["R_oracle"]["lambda_max"]:.2f} '
    f'(GMP-Cor undefined: no cells to scramble)',
    '',
    '-' * 78,
    'STABILITY (independent count replicate, same ground-truth R)',
    f'  genes in common      : {stability["n_genes_common"]}',
    f'  per-mode |r|         : '
    f'{[round(x, 3) for x in stability["per_mode_abs_r"]]}',
    f'  {p["n_top"]}-mode subspace overlap: {stability["subspace_overlap"]:.3f}',
    '',
    '=' * 78,
    'Reading the table: best_block_mass is the fraction of a mode\'s squared loading',
    'that falls in its best-matching true block; [null] is what a random gene set of',
    'that block\'s size would carry. eff.blocks is exp(entropy) over blocks — 1 means',
    'the mode IS one block, hundreds means it is smeared over the whole panel.',
    'projection is how much of a block\'s coherent direction lies inside the 5-mode',
    'subspace (1 = fully captured, 0 = invisible).',
    '=' * 78,
]
summary = '\n'.join(lines)
with open(TEXT_LOG, 'w', encoding='utf-8') as fh:
    fh.write(summary + '\n')
print('\n' + summary)

# ── Figure ───────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(STAGES)))

# (A) per-mode best-block mass across the ladder
ax = axes[0, 0]
w = 0.2
modes = np.arange(1, p['n_top'] + 1)
for j, name in enumerate(STAGES):
    vals = [m['best_block_mass'] for m in results[name]['modes']]
    ax.bar(modes + (j - 1.5) * w, vals, width=w, color=colors[j],
           label=STAGE_LABEL[name].replace('\n', ' '))
nulls = [m['best_block_mass_null'] for m in results['observed']['modes']]
ax.axhline(np.mean(nulls), ls='--', c='0.4', lw=1,
           label='random-gene-set null')
ax.set_xticks(modes)
ax.set_xlabel('mode')
ax.set_ylabel('squared loading mass in best block')
ax.set_title('A  Is each mode one block?', loc='left', fontweight='bold')
ax.legend(fontsize=7, frameon=False)

# (B) effective number of blocks per mode
ax = axes[0, 1]
for j, name in enumerate(STAGES):
    vals = [m['effective_n_blocks'] for m in results[name]['modes']]
    ax.plot(modes, vals, 'o-', color=colors[j], label=STAGE_LABEL[name].replace('\n', ' '))
ax.axhline(1, ls='--', c='0.4', lw=1)
ax.set_yscale('log')
ax.set_xticks(modes)
ax.set_xlabel('mode')
ax.set_ylabel('effective number of blocks, exp(H)')
ax.set_title('B  How many blocks does a mode mix?', loc='left', fontweight='bold')
ax.legend(fontsize=7, frameon=False)

# (C) per-block projection vs block size
ax = axes[1, 0]
for j, name in enumerate(STAGES):
    b = results[name]['blocks']
    if not b:
        continue
    ax.scatter([x['size_kept'] for x in b], [x['projection'] for x in b],
               s=18, alpha=0.7, color=colors[j],
               label=STAGE_LABEL[name].replace('\n', ' '))
ax.set_xscale('log')
ax.set_xlabel(f'true block size (genes, >= {p["min_block_size"]})')
ax.set_ylabel(f'projection onto {p["n_top"]}-mode subspace')
ax.set_ylim(-0.02, 1.02)
ax.set_title('C  Are the true blocks in the subspace?', loc='left', fontweight='bold')
ax.legend(fontsize=7, frameon=False)

# (D) the five strongest blocks, tracked individually across the ladder
ax = axes[1, 1]
# rank blocks by size * strength -- the proxy for how large an eigenvalue a block
# generates, i.e. which blocks a 5-mode analysis could plausibly reach at all
rank = np.argsort(sizes * strengths)[::-1][:p['n_top']]
xs = np.arange(len(STAGES))
cmap = plt.cm.plasma(np.linspace(0.05, 0.75, len(rank)))
for c, b_id in zip(cmap, rank):
    ys = []
    for name in STAGES:
        row = next((x for x in results[name]['blocks'] if x['block'] == b_id), None)
        ys.append(row['projection'] if row else np.nan)
    ax.plot(xs, ys, 'o-', color=c,
            label=f'block #{b_id}: {sizes[b_id]}g, w={strengths[b_id]:.2f}')
ax.set_xticks(xs)
ax.set_xticklabels([STAGE_LABEL[n] for n in STAGES], fontsize=7)
ax.set_ylabel('projection onto 5-mode subspace')
ax.set_ylim(-0.02, 1.02)
ax.set_title('D  The 5 strongest blocks down the ladder', loc='left', fontweight='bold')
ax.legend(fontsize=7, frameon=False)

fig.suptitle(f'Block recovery from the top-{p["n_top"]} eigenvector loadings '
             f'(alpha = {p["alpha"]}, {p["n_cells"]} cells x {p["n_genes"]} genes)',
             fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(SVG_FIG, format='svg', bbox_inches='tight')
fig.savefig(PNG_FIG, dpi=200, bbox_inches='tight')
plt.close(fig)

print(f'\nJSON log : {JSON_LOG}')
print(f'Text log : {TEXT_LOG}')
print(f'Figure   : {SVG_FIG}')
print(f'           {PNG_FIG}')
