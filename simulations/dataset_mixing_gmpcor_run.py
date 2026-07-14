"""
Dataset-mixing GMP-Cor runner.

Concatenates two experimental cell-gene matrices, subsamples to a fixed
n_cells x n_genes matrix, computes GMP-Cor, and plots the CCDF eigenvalue
spectrum in the project's figure-3 style.

Subsampling procedure:
  1. Load both matrices. Optionally drop non-gene columns (``Unnamed*`` /
     ``INTR_*``) and case-fold gene names so e.g. ``accA`` and ``acca`` align
     across datasets (duplicate columns are collapsed by summing).
  2. Build the gene space as the UNION or the SHARED (intersection) set of the
     two datasets' genes. Report how many genes each contributes.
  3. Cells: take the ``n_cells_per_dataset`` highest-total-count cells from EACH
     dataset, reindexed onto the gene space (missing genes -> 0), then stacked.
  4. Genes: keep the ``n_genes`` most highly variable by Fano factor
     (variance / mean), computed on the selected cells. If the gene space is
     already <= n_genes, all genes are kept.
  5. GMP-Cor = sum( lambda_i - lambda_max^scr ) over eigenvalues above the
     scrambled noise threshold, using src.analysis_functions.get_eig_dist.

Usage (from repo root or simulations/):
    python simulations/dataset_mixing_gmpcor_run.py \
        --file1 EXP_biorep_t0A_filtered.csv \
        --file2 VAPC_biorep_t2A_filtered.csv \
        --data-dir data_for_umap --gene-space shared --n-genes 2000

Outputs (results/simulation_results/):
  - figures/<base>.svg / .png   — CCDF eigenvalue plot
  - raw/<base>.txt              — human-readable summary
  - logs/<base>.json            — full parameter + results log
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

# Allow running directly from simulations/ or from repo root
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

import src.analysis_functions as af  # noqa: E402  (pure numpy, no umap dependency)

# Reporter / plasmid-marker genes excluded from every panel (case-insensitive)
EXCLUDE_GENES = {'gfp', 'laci', 'ampr'}


# ── Parameters (CLI overridable) ─────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--file1', default='EXP_biorep_t0A_filtered.csv',
                   help='first dataset CSV (cells x genes)')
    p.add_argument('--file2', default='VAPC_biorep_t2A_filtered.csv',
                   help='second dataset CSV (cells x genes)')
    p.add_argument('--data-dir', default='data_for_umap',
                   help='directory (relative to repo root, or absolute) holding the CSVs')
    p.add_argument('--gene-space', choices=['shared', 'union'], default='shared',
                   help="build gene space from the intersection ('shared') or 'union'")
    p.add_argument('--n-cells-per-dataset', type=int, default=500,
                   help='number of highest-total-count cells to take from each dataset')
    p.add_argument('--n-genes', type=int, default=2000,
                   help='number of genes to keep by Fano-factor cutoff')
    p.add_argument('--seed', type=int, default=0, help='RNG seed (scrambling)')
    p.add_argument('--no-case-fold', dest='case_fold', action='store_false',
                   help='do NOT case-fold gene names before matching')
    p.add_argument('--keep-nongene', dest='drop_nongene', action='store_false',
                   help='do NOT drop Unnamed*/INTR_* columns')
    return p.parse_args()


# ── Data loading ─────────────────────────────────────────────────────────────

def load_matrix(path, fname, case_fold=True, drop_nongene=True):
    d = pd.read_csv(path, index_col=0).fillna(0.0)
    if drop_nongene:
        drop = [c for c in d.columns
                if str(c).lower().startswith('unnamed') or str(c).startswith('INTR_')]
        if drop:
            d = d.drop(columns=drop)
            print(f'  {fname}: dropped {len(drop)} non-gene columns: {drop}')
    if case_fold:
        d.columns = [str(c).casefold() for c in d.columns]
        if d.columns.duplicated().any():
            # collapse duplicate gene columns by summing counts
            d = d.T.groupby(level=0).sum().T
    return d  # rows = cells, cols = genes


def top_cells(d, n):
    n = min(n, d.shape[0])
    keep = d.sum(axis=1).sort_values(ascending=False).index[:n]
    return d.loc[keep]


# ── CCDF plot (figure-3 style) ───────────────────────────────────────────────

def plot_ccdf(ax, data1, data2, title, fsize=12, signal_color='steelblue'):
    data1 = data1[data1 > 0]
    data2 = data2[data2 > 0]
    x2 = float(np.max(data2))
    d1s = np.sort(data1)
    d2s = np.sort(data2)
    p1 = len(d1s)
    ccdf1 = 1 - np.arange(1, p1 + 1) / p1 + 1 / p1
    p2 = len(d2s)
    ccdf2 = 1 - np.arange(1, p2 + 1) / p2 + 1 / p2
    noise = d1s < x2
    ax.loglog(d1s[noise], ccdf1[noise], '.', linestyle='-',
              color='darkgray', alpha=0.7, label='spurious', markersize=3)
    ax.loglog(d1s[~noise], ccdf1[~noise], '.', linestyle='-',
              color=signal_color, label='signal', markersize=3)
    ax.loglog(d2s, ccdf2, '.', linestyle='-',
              color='black', alpha=0.5, label='scrambled', markersize=3)
    ax.set_xlim([0.1, 40])
    ax.axvline(x2, color='k', linestyle='--', alpha=0.6)
    ax.text(x2 * 1.1, 0.8, r'$\lambda_\mathrm{max}^\mathrm{scr}$',
            fontsize=fsize - 2, va='center', ha='left', color='k', alpha=0.7,
            transform=ax.get_xaxis_transform())
    ax.set_xlabel(r'$\lambda$', fontsize=fsize)
    ax.set_ylabel('CCDF', fontsize=fsize)
    ax.set_title(title, fontsize=fsize)
    ax.legend(fontsize=fsize - 2)
    ax.tick_params(labelsize=fsize - 2)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    np.random.seed(args.seed)

    data_dir = args.data_dir if os.path.isabs(args.data_dir) \
        else os.path.join(_REPO_ROOT, args.data_dir)
    shared_only = args.gene_space == 'shared'

    print('=' * 65)
    print('Dataset-mixing GMP-Cor')
    print('=' * 65)

    df1 = load_matrix(os.path.join(data_dir, args.file1), args.file1,
                      args.case_fold, args.drop_nongene)
    df2 = load_matrix(os.path.join(data_dir, args.file2), args.file2,
                      args.case_fold, args.drop_nongene)
    print(f'{args.file1}: {df1.shape[0]} cells x {df1.shape[1]} genes')
    print(f'{args.file2}: {df2.shape[0]} cells x {df2.shape[1]} genes')

    genes1, genes2 = set(df1.columns), set(df2.columns)
    shared = genes1 & genes2
    only1, only2 = genes1 - genes2, genes2 - genes1
    gene_space = sorted(shared) if shared_only else sorted(genes1 | genes2)
    gene_mode = 'SHARED-ONLY (intersection)' if shared_only else 'UNION'

    # Drop reporter / plasmid-marker genes (case-insensitive)
    excluded = [g for g in gene_space if str(g).casefold() in EXCLUDE_GENES]
    if excluded:
        gene_space = [g for g in gene_space if str(g).casefold() not in EXCLUDE_GENES]
        print(f'excluded reporter genes: {excluded}')

    print(f'\n--- gene space: {gene_mode} ---')
    print(f'shared genes      : {len(shared)}')
    print(f'only in dataset 1 : {len(only1)}')
    print(f'only in dataset 2 : {len(only2)}')
    print(f'gene space size   : {len(gene_space)}')

    # Cell selection: top total-count cells from each dataset, on the gene space
    sel1 = top_cells(df1, args.n_cells_per_dataset).reindex(columns=gene_space, fill_value=0.0)
    sel2 = top_cells(df2, args.n_cells_per_dataset).reindex(columns=gene_space, fill_value=0.0)
    combined = pd.concat([sel1, sel2], axis=0)
    print(f'\nselected {sel1.shape[0]} + {sel2.shape[0]} cells; '
          f'combined {combined.shape[0]} cells x {combined.shape[1]} genes')

    # Gene selection: top n_genes by Fano factor (var / mean)
    mat = combined.to_numpy(dtype=float)
    mean = mat.mean(axis=0)
    var = mat.var(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        fano = np.where(mean > 0, var / mean, 0.0)
    n_genes = min(args.n_genes, mat.shape[1])
    top_idx = np.sort(np.argsort(fano)[::-1][:n_genes])
    final = mat[:, top_idx]
    print(f'final matrix: {final.shape[0]} cells x {final.shape[1]} genes')

    # GMP-Cor
    pcs, pcs1, frac_nonzero = af.get_eig_dist(
        final, norm=True, log=False, norm_method='sum', norm_sum=1)
    max_ev = float(np.max(pcs))
    max_ev_scr = float(np.max(pcs1))
    signal = pcs > max_ev_scr
    gmp_cor = float(np.sum(pcs[signal] - max_ev_scr))
    n_signal = int(np.sum(signal))
    print('\n--- GMP-Cor ---')
    print(f'fraction non-zero : {frac_nonzero:.4f}')
    print(f'lambda_max        : {max_ev:.4f}')
    print(f'lambda_max^scr    : {max_ev_scr:.4f}')
    print(f'# signal eigenvals: {n_signal}')
    print(f'GMP-Cor           : {gmp_cor:.4f}')

    # ── Output paths ─────────────────────────────────────────────────────────
    sim_results = os.path.join(_REPO_ROOT, 'results', 'simulation_results')
    fig_dir = os.path.join(sim_results, 'figures')
    raw_dir = os.path.join(sim_results, 'raw')
    log_dir = os.path.join(sim_results, 'logs')
    for d in (fig_dir, raw_dir, log_dir):
        os.makedirs(d, exist_ok=True)

    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    tag1 = os.path.splitext(args.file1)[0]
    tag2 = os.path.splitext(args.file2)[0]
    mode_tag = 'shared' if shared_only else 'union'
    base = f'mixing_{tag1}__{tag2}_{mode_tag}_{final.shape[0]}x{final.shape[1]}_{ts}'

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    title = (f'{tag1} + {tag2}\n({mode_tag} genes)  GMP-Cor = {gmp_cor:.2f}')
    plot_ccdf(ax, pcs.copy(), pcs1.copy(), title)
    fig.tight_layout()
    svg_path = os.path.join(fig_dir, base + '.svg')
    png_path = os.path.join(fig_dir, base + '.png')
    fig.savefig(svg_path, transparent=True)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    # ── Text summary ─────────────────────────────────────────────────────────
    summary = f"""Dataset-mixing GMP-Cor
================================
run timestamp       : {ts}
script              : {os.path.abspath(__file__)}
datasets            : {args.file1} + {args.file2}
data_dir            : {data_dir}
dataset1 shape      : {df1.shape[0]} cells x {df1.shape[1]} genes
dataset2 shape      : {df2.shape[0]} cells x {df2.shape[1]} genes

preprocessing
  case_fold         : {args.case_fold}
  drop_nongene      : {args.drop_nongene}
  excluded genes    : {excluded if excluded else 'none'}

gene space          : {gene_mode}
  shared genes      : {len(shared)}
  only in dataset1  : {len(only1)}
  only in dataset2  : {len(only2)}
  gene space size   : {len(gene_space)}

subsampling
  cells/dataset     : {args.n_cells_per_dataset} (highest total counts)
  gene selection    : top {n_genes} by Fano factor (var/mean)
  final matrix      : {final.shape[0]} cells x {final.shape[1]} genes

GMP-Cor
  fraction non-zero : {frac_nonzero:.4f}
  lambda_max        : {max_ev:.4f}
  lambda_max^scr    : {max_ev_scr:.4f}
  # signal eigenvals: {n_signal}
  GMP-Cor           : {gmp_cor:.4f}

outputs
  SVG  : {svg_path}
  PNG  : {png_path}
"""
    txt_path = os.path.join(raw_dir, base + '.txt')
    with open(txt_path, 'w') as f:
        f.write(summary)

    # ── JSON log ─────────────────────────────────────────────────────────────
    log = {
        'experiment': 'dataset_mixing_gmpcor',
        'timestamp': ts,
        'seed': args.seed,
        'params': {
            'file1': args.file1, 'file2': args.file2, 'data_dir': data_dir,
            'gene_space': args.gene_space,
            'n_cells_per_dataset': args.n_cells_per_dataset,
            'n_genes': args.n_genes, 'case_fold': args.case_fold,
            'drop_nongene': args.drop_nongene,
            'excluded_genes': excluded,
            'cell_selection': 'top total counts',
            'gene_selection': 'top fano factor',
            'norm': True, 'log': False, 'norm_method': 'sum', 'norm_sum': 1,
            'get_eig_dist_reps': 10,
        },
        'dataset1_shape': list(df1.shape),
        'dataset2_shape': list(df2.shape),
        'gene_space_info': {
            'mode': gene_mode, 'shared': len(shared),
            'only1': len(only1), 'only2': len(only2), 'size': len(gene_space),
        },
        'final_shape': list(final.shape),
        'results': {
            'fraction_non_zero': frac_nonzero,
            'lambda_max': max_ev,
            'lambda_max_scrambled': max_ev_scr,
            'n_signal_eigenvalues': n_signal,
            'gmp_cor': gmp_cor,
        },
    }
    json_path = os.path.join(log_dir, base + '.json')
    with open(json_path, 'w') as f:
        json.dump(log, f, indent=2)

    print(f'\nsaved figure : {svg_path}')
    print(f'saved figure : {png_path}')
    print(f'saved summary: {txt_path}')
    print(f'saved log    : {json_path}')


if __name__ == '__main__':
    main()
