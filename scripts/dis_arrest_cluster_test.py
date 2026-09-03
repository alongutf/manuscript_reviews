"""Are the two dis-arrest clusters biology, or do they just track expression level?

Builds a UMAP from the two dis-arrest samples alone (sample_13a + sample_15a,
data_for_paper, which share a single 2042-gene panel), tunes the leiden
resolution until exactly two clusters are returned, then asks two questions:

  (1) Do the two clusters track sequencing depth / detection rather than a
      transcriptional program?  Measured by the AUC of total counts and of
      genes-detected between the clusters, plus a logistic model on log-depth
      alone.
  (2) Is GMP-Cor within each cluster meaningful?  Computed per cluster at a
      matched cell number and on the identical gene panel, alongside three
      controls -- a median split on depth, a random split, and the batch split.

dGMP = 1 - GMP-Cor(group-centered) / GMP-Cor(raw): the fraction of GMP-Cor that
is carried by the between-group mean difference rather than within-group
covariance.

Outputs (results/dis_arrest_clusters/):
  figures/  <stem>_<timestamp>.svg / .png
  raw/      <stem>_<timestamp>.txt   human-readable summary
  logs/     <stem>_<timestamp>.json  every parameter and every number
            <stem>_<timestamp>_cells.csv  per-cell embedding + labels + depth
"""
import os
import sys
import json
import argparse
import datetime
import contextlib
import io

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import src.analysis_functions as af

OUTDIR = os.path.join(ROOT, 'results', 'dis_arrest_clusters')
PAPER = os.path.join(ROOT, 'data_for_paper')

SAMPLES = {'dis1': 'sample_13a_filtered.csv', 'dis2': 'sample_15a_filtered.csv'}


# ---------------------------------------------------------------- metrics ---
def gmp_cor(m, norm_sum=1.0, rep=10, seed=None):
    """GMP-Cor = sum_i max(lambda_i - lambda_max^scrambled, 0) for a cells x genes matrix.

    Same definition and preprocessing as af.get_eig_dist (rep=10 scrambles, threshold =
    max of the rep-averaged scrambled spectrum), evaluated through the smaller Gram
    matrix so that the many repeats below stay affordable.
    """
    z = _z(m, norm=True, norm_sum=norm_sum)
    ev = _spectrum(z)
    rng = np.random.default_rng(seed)
    acc = np.zeros_like(ev)
    for _ in range(rep):
        acc += _spectrum(_scramble(z, rng))
    thr = float((acc / rep)[0])
    excess = ev[ev > thr] - thr
    return {
        'gmp_cor': float(excess.sum()),
        'lambda_max': float(ev[0]),
        'lambda_max_scrambled': thr,
        'n_modes_above_threshold': int(excess.size),
    }


def _z(m, norm=True, norm_sum=1.0):
    """The exact preprocessing get_eig_dist applies, returned as a matrix."""
    m = np.asarray(m, dtype=float)
    m = m[:, (m > 0).sum(axis=0) >= 1]
    m = m[(m > 0).sum(axis=1) >= 1, :]
    if norm:
        m = af.normalize(m, method='sum', target_sum=norm_sum)
    return af.z_transform(m)


def _spectrum(z, vectors=False):
    """Non-zero eigenvalues of the gene-gene covariance, via the smaller Gram matrix."""
    n = z.shape[0]
    g = z @ z.T if n <= z.shape[1] else z.T @ z
    if vectors:
        w, v = np.linalg.eigh(g)
        return (w / n)[::-1], v[:, ::-1]
    return (np.linalg.eigvalsh(g) / n)[::-1]


def _scramble(z, rng):
    """Independently permute the rows within every column (af.scramble, vectorized)."""
    idx = rng.random(z.shape).argsort(axis=0)
    return np.take_along_axis(z, idx, axis=0)


def leading_mode_concentration(m, norm_sum=1.0):
    """How many cells carry the leading mode.

    Only meaningful when the Gram matrix is the smaller one (n <= p), which holds for
    every group here; then the leading eigenvector lives in cell space directly.
    Returns the participation ratio (an effective number of contributing cells), the
    same as a fraction of n, and the share of the mode's weight held by the top 1% of
    cells. A mode driven by a handful of outlier cells has a participation ratio of
    order 1-10 rather than of order n.
    """
    z = _z(m, norm=True, norm_sum=norm_sum)
    if z.shape[0] > z.shape[1]:
        return None
    _, vecs = _spectrum(z, vectors=True)
    v = vecs[:, 0]
    w = v ** 2
    w = w / w.sum()
    n = w.size
    k = max(1, int(round(0.01 * n)))
    return {
        'n_cells': int(n),
        'participation_ratio': float(1.0 / np.sum(w ** 2)),
        'participation_fraction': float(1.0 / np.sum(w ** 2) / n),
        'top_1pct_cells_weight_share': float(np.sort(w)[::-1][:k].sum()),
    }


def permutation_null(m, n_perm=200, norm_sum=1.0, seed=0):
    """The paper's permutation test (add_permutation_metrics.py) applied to one matrix.

    Returns the observed spectrum, the null distribution of lambda_1 over n_perm
    column-wise scrambles, the empirical p-value, and GMP-Cor referenced to the
    null mean with its sqrt(N)*sigma uncertainty.
    """
    z = _z(m, norm=True, norm_sum=norm_sum)
    ev = _spectrum(z)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for b in range(n_perm):
        null[b] = _spectrum(_scramble(z, rng))[0]
    thr = float(null.mean())
    sd = float(null.std(ddof=1))
    n_above = int((ev > thr).sum())
    return {
        'n_perm': int(n_perm),
        'lambda_max': float(ev[0]),
        'null_mean': thr, 'null_sd': sd,
        'p_empirical': float((1 + int((null >= ev[0]).sum())) / (n_perm + 1)),
        'n_above_threshold': n_above,
        'gmp_cor': float(np.sum(np.maximum(ev - thr, 0.0))),
        'gmp_cor_ci': float(np.sqrt(n_above) * sd),
    }


def group_center(m, labels, norm_sum=1.0):
    """Row-normalize, then subtract each group's own gene means."""
    x = af.normalize(np.asarray(m, dtype=float), method='sum', target_sum=norm_sum)
    out = x.copy()
    for g in np.unique(labels):
        sel = labels == g
        out[sel] -= x[sel].mean(axis=0)
    return out


def gmp_cor_centered(m, labels, norm_sum=1.0):
    """GMP-Cor after removing the between-group mean difference (no re-normalization)."""
    x = group_center(m, labels, norm_sum=norm_sum)
    with contextlib.redirect_stdout(io.StringIO()):
        pcs, pcs1, _ = af.get_eig_dist(x, norm=False, log=False)
    thr = float(pcs1.max())
    return float(np.sum(np.maximum(pcs - thr, 0.0)))


def auc(x, labels):
    """Rank AUC of continuous x separating labels==1 from labels==0."""
    from scipy.stats import rankdata
    x = np.asarray(x, dtype=float)
    labels = np.asarray(labels)
    pos = labels == 1
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return float('nan')
    r = rankdata(x)
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def cohens_d(x, labels):
    """Standardized mean difference (label==1 minus label==0), pooled SD (ddof=1)."""
    x = np.asarray(x, dtype=float)
    labels = np.asarray(labels)
    a, b = x[labels == 1], x[labels == 0]
    s = np.sqrt(((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1))
                / max(a.size + b.size - 2, 1))
    return float((a.mean() - b.mean()) / s) if s > 0 else float('nan')


# ------------------------------------------------------------- embedding ---
def build_embedding(counts, batch, args):
    """Published preprocessing (scanpy_analysis.ipynb) restricted to the dis pair."""
    import scanpy as sc
    import anndata as ad

    adata = ad.AnnData(X=counts.values.astype(np.float32),
                       obs=pd.DataFrame({'batch': batch}, index=counts.index),
                       var=pd.DataFrame(index=counts.columns))
    adata.obs['total_counts'] = counts.values.sum(axis=1)
    adata.obs['n_detected'] = (counts.values > 0).sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    n_hvg = min(args.n_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=True)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=args.n_comps)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.umap(adata, min_dist=args.min_dist, random_state=args.umap_seed)

    # bisect the leiden resolution until exactly two clusters come back
    trace, res_used, lo, hi = [], None, args.res_lo, args.res_hi
    for _ in range(args.max_res_steps):
        res = 0.5 * (lo + hi)
        sc.tl.leiden(adata, resolution=res, key_added='leiden', random_state=args.leiden_seed)
        k = int(adata.obs['leiden'].nunique())
        trace.append({'resolution': float(res), 'n_clusters': k})
        if k == 2:
            res_used = float(res)
            break
        if k < 2:
            lo = res
        else:
            hi = res
    if res_used is None:
        raise RuntimeError('no resolution in [%g, %g] gave exactly 2 clusters; trace=%s'
                           % (args.res_lo, args.res_hi, trace))
    return adata, res_used, trace, n_hvg


# ------------------------------------------------------------------ main ---
def main():
    """Run the full dis-arrest two-cluster test end to end.

    Loads the two dis-arrest sample matrices on their shared gene panel, builds the
    published-style embedding and bisects leiden resolution to exactly two clusters
    (build_embedding), then answers question 1 (does the split track sequencing
    depth?) and question 2 (is GMP-Cor meaningful within each cluster, at matched n,
    against a permutation null and against depth/random/batch control splits), and
    finally writes the JSON log, per-cell CSV, text summary and figure described in
    the module docstring.
    """
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stem', default='dis_arrest_cluster_test')
    ap.add_argument('--target-sum', type=float, default=3507,
                    help='per-cell normalization for the embedding (published value for dis)')
    ap.add_argument('--min-cells', type=int, default=3)
    ap.add_argument('--n-hvg', type=int, default=2000)
    ap.add_argument('--n-comps', type=int, default=50)
    ap.add_argument('--n-neighbors', type=int, default=40)
    ap.add_argument('--n-pcs', type=int, default=40)
    ap.add_argument('--min-dist', type=float, default=0.3)
    ap.add_argument('--umap-seed', type=int, default=0)
    ap.add_argument('--leiden-seed', type=int, default=0)
    ap.add_argument('--res-lo', type=float, default=0.01)
    ap.add_argument('--res-hi', type=float, default=2.0)
    ap.add_argument('--max-res-steps', type=int, default=30)
    ap.add_argument('--permutations', type=int, default=200,
                    help='B for the per-group permutation test on lambda_1 (0 to skip)')
    ap.add_argument('--repeats', type=int, default=5,
                    help='seeded cell draws per group for the matched-n GMP-Cor')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--norm-sum', type=float, default=1.0,
                    help='row normalization for GMP-Cor; z-transform makes it scale-free')
    args = ap.parse_args()

    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    for sub in ('figures', 'raw', 'logs'):
        os.makedirs(os.path.join(OUTDIR, sub), exist_ok=True)

    # ---- load: the two dis samples share one gene panel by construction ----
    mats, batches = [], []
    for key, fname in SAMPLES.items():
        m = pd.read_csv(os.path.join(PAPER, fname), index_col=0)
        m.index = [str(b) + '-' + key for b in m.index]
        mats.append(m)
        batches += [key] * m.shape[0]
    shared = [g for g in mats[0].columns if g in set(mats[1].columns)]
    if len(shared) != mats[0].shape[1] or len(shared) != mats[1].shape[1]:
        print('note: panels differ; using the %d shared genes' % len(shared))
    counts = pd.concat([m[shared] for m in mats], axis=0)
    batch = np.array(batches)

    adata, res_used, res_trace, n_hvg = build_embedding(counts, batch, args)

    cells = pd.DataFrame({
        'batch': adata.obs['batch'].values,
        'leiden': adata.obs['leiden'].values.astype(str),
        'total_counts': adata.obs['total_counts'].values,
        'n_detected': adata.obs['n_detected'].values,
        'UMAP_1': adata.obsm['X_umap'][:, 0],
        'UMAP_2': adata.obsm['X_umap'][:, 1],
        'PC_1': adata.obsm['X_pca'][:, 0],
        'PC_2': adata.obsm['X_pca'][:, 1],
    }, index=adata.obs_names)

    # orient: the cluster labelled 1 is the deeper one, so AUCs come out positive
    med = cells.groupby('leiden')['total_counts'].median()
    deep = med.idxmax()
    y = (cells['leiden'] == deep).values.astype(int)
    cells['cluster'] = np.where(y == 1, 'high-depth', 'low-depth')

    X = counts.loc[cells.index].values.astype(float)

    # ------------------------------- question 1: do the clusters = depth? ---
    from scipy.stats import mannwhitneyu, spearmanr
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_predict

    logd = np.log10(cells['total_counts'].values)
    depth_tests = {}
    for name, v in (('total_counts', cells['total_counts'].values.astype(float)),
                    ('n_detected', cells['n_detected'].values.astype(float))):
        depth_tests[name] = {
            'median_high_depth_cluster': float(np.median(v[y == 1])),
            'median_low_depth_cluster': float(np.median(v[y == 0])),
            'fold_change': float(np.median(v[y == 1]) / np.median(v[y == 0])),
            'auc': auc(v, y),
            'cohens_d': cohens_d(v, y),
            'mannwhitney_p': float(mannwhitneyu(v[y == 1], v[y == 0]).pvalue),
        }
    lr = LogisticRegression()
    pred = cross_val_predict(lr, logd.reshape(-1, 1), y, cv=5, method='predict_proba')[:, 1]
    depth_only = {
        'cv_auc_logdepth_predicts_cluster': auc(pred, y),
        'cv_accuracy': float(((pred > 0.5).astype(int) == y).mean()),
    }
    embedding_vs_depth = {
        'spearman_' + c + '_vs_log_total_counts':
            float(spearmanr(cells[c].values, logd).statistic)
        for c in ('UMAP_1', 'UMAP_2', 'PC_1', 'PC_2')
    }
    batch_table = pd.crosstab(cells['batch'], cells['cluster'])
    batch_auc = auc(y, (cells['batch'] == 'dis2').values.astype(int))

    # ------------------------- question 2: GMP-Cor within each cluster -----
    rng_master = np.random.default_rng(args.seed)
    n_match = int(min((y == 1).sum(), (y == 0).sum()))

    def matched(idx_pool, label, seed_offset):
        """GMP-Cor over `repeats` seeded draws of n_match cells from idx_pool."""
        vals, rows = [], []
        for rep in range(args.repeats):
            rng = np.random.default_rng(args.seed + 1000 * seed_offset + rep)
            take = rng.choice(idx_pool, size=n_match, replace=False)
            r = gmp_cor(X[take], norm_sum=args.norm_sum, seed=args.seed + rep)
            r['rep'] = rep
            rows.append(r)
            vals.append(r['gmp_cor'])
        vals = np.asarray(vals)
        # the rep-0 cell set is the canonical one; the permutation test runs on it
        rng0 = np.random.default_rng(args.seed + 1000 * seed_offset)
        take0 = rng0.choice(idx_pool, size=n_match, replace=False)
        perm = (permutation_null(X[take0], n_perm=args.permutations,
                                 norm_sum=args.norm_sum, seed=args.seed + seed_offset)
                if args.permutations > 0 else None)
        return {
            'label': label, 'n_cells': n_match, 'n_pool': int(len(idx_pool)),
            'permutation_test': perm,
            'leading_mode': leading_mode_concentration(X[take0], norm_sum=args.norm_sum),
            'gmp_cor_draws': [float(v) for v in vals],
            'gmp_cor_mean': float(np.mean(vals)),
            'gmp_cor_std': float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
            'gmp_cor_min': float(vals.min()), 'gmp_cor_max': float(vals.max()),
            'gmp_cor_fold_range': (float(vals.max() / vals.min()) if vals.min() > 0
                                   else float('inf')),
            'lambda_max_mean': float(np.mean([r['lambda_max'] for r in rows])),
            'lambda_max_scrambled_mean': float(np.mean([r['lambda_max_scrambled'] for r in rows])),
            'detected_per_cell': float(cells['n_detected'].values[idx_pool].mean()),
            'total_counts_per_cell': float(cells['total_counts'].values[idx_pool].mean()),
            'reps': rows,
        }

    all_idx = np.arange(X.shape[0])
    depth_split = (cells['total_counts'].values
                   > np.median(cells['total_counts'].values)).astype(int)
    rand_split = rng_master.permutation(
        np.repeat([0, 1], [X.shape[0] // 2, X.shape[0] - X.shape[0] // 2]))

    groups = {
        'all cells': all_idx,
        'high-depth cluster': all_idx[y == 1],
        'low-depth cluster': all_idx[y == 0],
        'depth split: upper half': all_idx[depth_split == 1],
        'depth split: lower half': all_idx[depth_split == 0],
        'random half A': all_idx[rand_split == 1],
        'random half B': all_idx[rand_split == 0],
        'dis1': all_idx[cells['batch'].values == 'dis1'],
        'dis2': all_idx[cells['batch'].values == 'dis2'],
    }
    per_group = {k: matched(v, k, i) for i, (k, v) in enumerate(groups.items())}

    # dGMP on all cells for each candidate partition
    raw_all = gmp_cor(X, norm_sum=args.norm_sum)['gmp_cor']
    # normalized (but not z-transformed) matrix, shared by every partition below: it
    # gives the group mean-difference axis and the AUC of projecting cells onto it
    xn = af.normalize(X, method='sum', target_sum=args.norm_sum)
    partitions = {
        'leiden (2 clusters)': y,
        'depth median split': depth_split,
        'batch (dis1 vs dis2)': (cells['batch'].values == 'dis2').astype(int),
        'random split': rand_split,
    }
    dgmp = {}
    for name, lab in partitions.items():
        cen = gmp_cor_centered(X, lab, norm_sum=args.norm_sum)
        # direction of the group-1 minus group-0 mean difference in gene space;
        # projecting cells onto it and taking the AUC checks how well that single
        # axis alone separates the two groups (a sanity check on dGMP's premise)
        axis = xn[lab == 1].mean(0) - xn[lab == 0].mean(0)
        dgmp[name] = {
            'gmp_cor_raw': raw_all, 'gmp_cor_group_centered': cen,
            'dGMP': float(1.0 - cen / raw_all) if raw_all > 0 else float('nan'),
            'group_axis_auc': auc(xn @ axis, lab),
        }

    # markers: is the split just "more genes detected"?
    from scipy.stats import ranksums
    stat = np.array([ranksums(xn[y == 1, j], xn[y == 0, j]).statistic
                     for j in range(xn.shape[1])])
    order = np.argsort(-np.abs(stat))[:50]
    markers = {
        'n_top': 50,
        'n_up_in_high_depth': int((stat[order] > 0).sum()),
        'n_up_in_low_depth': int((stat[order] < 0).sum()),
        'top_genes': [{'gene': str(counts.columns[j]), 'z': float(stat[j])} for j in order[:15]],
    }

    # ------------------------------------------------------------- output --
    log = {
        'generated': stamp, 'script': os.path.basename(__file__),
        'parameters': vars(args),
        'inputs': {k: os.path.join('data_for_paper', v) for k, v in SAMPLES.items()},
        'matrix': {'n_cells': int(X.shape[0]), 'n_genes_panel': int(X.shape[1]),
                   'n_genes_hvg_embedding': int(n_hvg)},
        'leiden': {'resolution': res_used, 'trace': res_trace,
                   'cluster_sizes': cells['cluster'].value_counts().to_dict()},
        'q1_depth': {'tests': depth_tests, 'depth_only_classifier': depth_only,
                     'embedding_vs_depth': embedding_vs_depth,
                     'batch_by_cluster': batch_table.to_dict(),
                     'auc_cluster_predicts_batch': batch_auc},
        'q2_gmp_cor': {'matched_n': n_match, 'per_group': per_group, 'dGMP': dgmp},
        'markers': markers,
    }
    jpath = os.path.join(OUTDIR, 'logs', args.stem + '_' + stamp + '.json')
    with open(jpath, 'w') as f:
        json.dump(log, f, indent=2, default=float)
    cpath = os.path.join(OUTDIR, 'logs', args.stem + '_' + stamp + '_cells.csv')
    cells.to_csv(cpath)

    # ---- summary -----------------------------------------------------------
    L = []
    L.append('Dis-arrest two-cluster test')
    L.append('generated ' + stamp)
    L.append('')
    L.append('Input : data_for_paper sample_13a + sample_15a, %d cells x %d genes (shared panel)'
             % (X.shape[0], X.shape[1]))
    L.append('Leiden: resolution %.4f gives exactly 2 clusters (%s)'
             % (res_used, cells['cluster'].value_counts().to_dict()))
    L.append('')
    L.append('Q1 -- do the clusters track expression level?')
    for name, d in depth_tests.items():
        L.append('  %-14s median %8.1f (low) vs %8.1f (high)  fold %5.2f  AUC %.3f  d %5.2f  p %.2e'
                 % (name, d['median_low_depth_cluster'], d['median_high_depth_cluster'],
                    d['fold_change'], d['auc'], d['cohens_d'], d['mannwhitney_p']))
    L.append('  log10(total counts) alone predicts the cluster label: CV AUC %.3f, accuracy %.3f'
             % (depth_only['cv_auc_logdepth_predicts_cluster'], depth_only['cv_accuracy']))
    for k, v in embedding_vs_depth.items():
        L.append('  %-38s %+.3f' % (k, v))
    L.append('  cluster vs batch:')
    for line in batch_table.to_string().split('\n'):
        L.append('    ' + line)
    L.append('  AUC(cluster separates dis1 from dis2) = %.3f' % batch_auc)
    L.append('  top-50 markers: %d up in the high-depth cluster, %d up in the low-depth cluster'
             % (markers['n_up_in_high_depth'], markers['n_up_in_low_depth']))
    L.append('')
    L.append('Q2 -- GMP-Cor per group, all at n = %d cells, p = %d' % (n_match, X.shape[1]))
    L.append('  Significance: the paper\'s permutation test on lambda_1 (B = %d column scrambles), '
             'run on one fixed cell set per group.' % args.permutations)
    L.append('  Sampling spread: mean +/- SD of GMP-Cor over %d independent cell draws.'
             % args.repeats)
    L.append('')
    L.append('  %-26s %5s %8s %8s %7s %7s %7s %17s %17s %14s %5s'
             % ('group', 'pool', 'det/cell', 'cnt/cell', 'lam_1', 'null_mu', 'perm_p',
                'GMP-Cor fixed set', 'over cell draws', 'draw range', 'PR'))
    for k, d in per_group.items():
        p = d['permutation_test']
        lm = d['leading_mode']
        pr = ('%5.1f' % lm['participation_ratio']) if lm else '    -'
        rng_s = '%5.2f - %-6.2f' % (d['gmp_cor_min'], d['gmp_cor_max'])
        if p is None:
            L.append('  %-26s %5d %8.1f %8.1f %7s %7s %7s %17s %8.2f +/- %-6.2f %14s %s'
                     % (k, d['n_pool'], d['detected_per_cell'], d['total_counts_per_cell'],
                        '-', '-', '-', '-', d['gmp_cor_mean'], d['gmp_cor_std'], rng_s, pr))
            continue
        L.append('  %-26s %5d %8.1f %8.1f %7.3f %7.3f %7.4f %8.2f +/- %-6.2f %8.2f +/- %-6.2f '
                 '%14s %s'
                 % (k, d['n_pool'], d['detected_per_cell'], d['total_counts_per_cell'],
                    p['lambda_max'], p['null_mean'], p['p_empirical'],
                    p['gmp_cor'], p['gmp_cor_ci'], d['gmp_cor_mean'], d['gmp_cor_std'],
                    rng_s, pr))
    L.append('  PR = participation ratio of the leading mode in cell space: the effective number '
             'of cells carrying it, out of n = %d.' % n_match)
    L.append('  perm_p is censored at 1/(B+1) = %.4f.' % (1.0 / (args.permutations + 1))
             if args.permutations > 0 else '')
    L.append('')
    L.append('  dGMP = fraction of GMP-Cor carried by the between-group mean difference')
    L.append('  (all %d cells, raw GMP-Cor = %.3f)' % (X.shape[0], raw_all))
    L.append('  %-26s %10s %8s %9s' % ('partition', 'centered', 'dGMP', 'axis AUC'))
    for k, d in dgmp.items():
        L.append('  %-26s %10.3f %8.3f %9.3f'
                 % (k, d['gmp_cor_group_centered'], d['dGMP'], d['group_axis_auc']))
    L.append('')
    a_cnt = depth_tests['total_counts']['auc']
    a_det = depth_tests['n_detected']['auc']
    d_leiden = dgmp['leiden (2 clusters)']['dGMP']
    d_depth = dgmp['depth median split']['dGMP']
    d_rand = dgmp['random split']['dGMP']
    hi = per_group['high-depth cluster']['gmp_cor_mean']
    lo = per_group['low-depth cluster']['gmp_cor_mean']
    L.append('Reading:')
    if max(a_cnt, a_det) > 0.9:
        L.append('  The split is essentially a depth split: AUC %.3f on total counts and %.3f on '
                 'genes detected.' % (a_cnt, a_det))
    elif max(a_cnt, a_det) > 0.75:
        L.append('  The split is strongly confounded with depth (AUC %.3f / %.3f) but not fully '
                 'explained by it.' % (a_cnt, a_det))
    else:
        L.append('  The split is not explained by depth (AUC %.3f on total counts, %.3f on genes '
                 'detected).' % (a_cnt, a_det))
    # Q2 verdict, driven by the permutation test rather than by the point estimates
    alpha = 0.05
    sig = {k: (d['permutation_test'] is not None
               and d['permutation_test']['p_empirical'] < alpha)
           for k, d in per_group.items()}
    named = ['all cells', 'high-depth cluster', 'low-depth cluster']
    if args.permutations > 0:
        yes = [k for k in named if sig[k]]
        no = [k for k in named if not sig[k]]
        if not yes:
            L.append('  No dis-arrest group carries a leading eigenvalue distinguishable from its '
                     'own permutation null at p < %.2f -- neither cluster, nor the pooled set. '
                     'GMP-Cor within these clusters is not a meaningful quantity to interpret.'
                     % alpha)
        elif no:
            L.append('  Significant against the permutation null: %s. Not significant: %s.'
                     % (', '.join(yes), ', '.join(no)))
        else:
            L.append('  All of %s exceed their permutation null at p < %.2f.'
                     % (', '.join(yes), alpha))
        ra, rb = per_group['random half A'], per_group['random half B']
        L.append('  Scale of the noise: two random halves of the same cells -- which differ by '
                 'construction only by which cells landed in them -- give GMP-Cor %.2f and %.2f '
                 '(perm_p %.3f and %.3f).'
                 % (ra['permutation_test']['gmp_cor'], rb['permutation_test']['gmp_cor'],
                    ra['permutation_test']['p_empirical'], rb['permutation_test']['p_empirical']))
        ac = per_group['all cells']
        L.append('  Over %d independent %d-cell draws from the same pool, GMP-Cor ranges %.2f to '
                 '%.2f (%.1f-fold). A cluster-to-cluster difference smaller than that is not '
                 'evidence of anything.'
                 % (args.repeats, n_match, ac['gmp_cor_min'], ac['gmp_cor_max'],
                    ac['gmp_cor_fold_range']))
        lm = per_group['all cells']['leading_mode']
        if lm is not None:
            L.append('  Why it is unstable: the leading mode is carried by an effective %.1f of '
                     '%d cells (top 1%% of cells hold %.0f%% of its weight), so which few cells a '
                     'draw happens to include sets the value.'
                     % (lm['participation_ratio'], lm['n_cells'],
                        100 * lm['top_1pct_cells_weight_share']))
    L.append('  Within-cluster GMP-Cor at matched n, averaged over %d cell draws: %.2f +/- %.2f '
             '(high-depth) vs %.2f +/- %.2f (low-depth); all cells %.2f +/- %.2f.'
             % (args.repeats, hi, per_group['high-depth cluster']['gmp_cor_std'],
                lo, per_group['low-depth cluster']['gmp_cor_std'],
                per_group['all cells']['gmp_cor_mean'], per_group['all cells']['gmp_cor_std']))
    if raw_all < 1.0:
        L.append('  dGMP (leiden %.3f, depth split %.3f, random split %.3f) is NOT interpretable '
                 'here: its denominator is the pooled GMP-Cor, %.3f, which is at the noise floor, '
                 'so the ratio is dominated by scramble noise and even the random split -- which '
                 'has no between-group difference by construction -- returns %.3f.'
                 % (d_leiden, d_depth, d_rand, raw_all, d_rand))
    else:
        L.append('  dGMP: leiden %.3f, depth median split %.3f, random split %.3f (the random '
                 'split is the null).' % (d_leiden, d_depth, d_rand))
    L.append('  GMP-Cor is extensive in p and depends on n through the Marchenko-Pastur edge, so '
             'only the matched-n, matched-p rows above are comparable to one another.')
    txt = '\n'.join(L)
    tpath = os.path.join(OUTDIR, 'raw', args.stem + '_' + stamp + '.txt')
    with open(tpath, 'w') as f:
        f.write(txt + '\n')
    print(txt)

    make_figure(cells, per_group, dgmp, depth_tests, res_used, X.shape,
                os.path.join(OUTDIR, 'figures', args.stem + '_' + stamp))
    print('\nwrote:\n  %s\n  %s\n  %s' % (jpath, cpath, tpath))


def make_figure(cells, per_group, dgmp, depth_tests, res_used, shape, stem):
    """Six-panel summary figure, written to <stem>.svg and <stem>.png (dpi=200).

    A: UMAP colored by leiden cluster.  B: UMAP colored by sample of origin.
    C: UMAP colored by log10 total counts (the depth confound).
    D: histogram of total counts per cluster, log-x.
    E: GMP-Cor (vs. permutation null when available, else draw mean +/- SD) per
       group, with the empirical p-value annotated and starred at p < 0.05.
    F: per-group GMP-Cor across the independent matched-n cell draws, one row per
       group, points jittered vertically only for readability.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(13.5, 8))
    cmap = {'high-depth': '#c0392b', 'low-depth': '#2980b9'}

    a = ax[0, 0]
    for k, sub in cells.groupby('cluster'):
        a.scatter(sub.UMAP_1, sub.UMAP_2, s=4, c=cmap[k], label='%s (n=%d)' % (k, len(sub)),
                  alpha=.6, linewidths=0)
    a.legend(fontsize=7, frameon=False)
    a.set_title('A  leiden, resolution %.3f' % res_used, fontsize=9, loc='left')

    a = ax[0, 1]
    for k, sub in cells.groupby('batch'):
        a.scatter(sub.UMAP_1, sub.UMAP_2, s=4, label='%s (n=%d)' % (k, len(sub)),
                  alpha=.6, linewidths=0)
    a.legend(fontsize=7, frameon=False)
    a.set_title('B  sample of origin', fontsize=9, loc='left')

    a = ax[0, 2]
    sc = a.scatter(cells.UMAP_1, cells.UMAP_2, s=4, c=np.log10(cells.total_counts),
                   cmap='viridis', alpha=.8, linewidths=0)
    plt.colorbar(sc, ax=a, label='log10 total counts', fraction=.046)
    a.set_title('C  sequencing depth', fontsize=9, loc='left')
    for a in ax[0]:
        a.set_xlabel('UMAP 1', fontsize=8)
        a.set_ylabel('UMAP 2', fontsize=8)
        a.tick_params(labelsize=7)

    a = ax[1, 0]
    bins = np.logspace(np.log10(cells.total_counts.min()),
                       np.log10(cells.total_counts.max()), 40)
    for k, sub in cells.groupby('cluster'):
        a.hist(sub.total_counts, bins=bins, alpha=.6, color=cmap[k], label=k)
    a.set_xscale('log')
    a.set_xlabel('total counts per cell', fontsize=8)
    a.set_ylabel('cells', fontsize=8)
    a.legend(fontsize=7, frameon=False)
    a.set_title('D  depth by cluster (AUC %.3f)' % depth_tests['total_counts']['auc'],
                fontsize=9, loc='left')
    a.tick_params(labelsize=7)

    a = ax[1, 1]
    keys = list(per_group.keys())
    has_perm = per_group[keys[0]]['permutation_test'] is not None
    vals = [(per_group[k]['permutation_test']['gmp_cor'] if has_perm
             else per_group[k]['gmp_cor_mean']) for k in keys]
    errs = [(per_group[k]['permutation_test']['gmp_cor_ci'] if has_perm
             else per_group[k]['gmp_cor_std']) for k in keys]
    cols = ['#555555' if k == 'all cells'
            else cmap.get(k.replace(' cluster', ''), '#7f8c8d') for k in keys]
    a.barh(range(len(keys)), vals, xerr=errs, color=cols, height=.7)
    if has_perm:
        for i, k in enumerate(keys):     # star the groups that beat their own null
            p = per_group[k]['permutation_test']['p_empirical']
            a.text(vals[i] + errs[i] + .3, i, 'p=%.3f%s' % (p, ' *' if p < .05 else ''),
                   va='center', fontsize=6, color='k')
    a.set_yticks(range(len(keys)))
    a.set_yticklabels(keys, fontsize=7)
    a.invert_yaxis()
    a.margins(x=.25)
    a.set_xlabel('GMP-Cor vs permutation null (error = sqrt(N)*sigma)', fontsize=8)
    a.set_title('E  GMP-Cor at matched n = %d, p = %d'
                % (per_group['all cells']['n_cells'], shape[1]), fontsize=9, loc='left')
    a.tick_params(labelsize=7)

    a = ax[1, 2]
    keys = list(per_group.keys())
    cols = ['#555555' if k == 'all cells'
            else cmap.get(k.replace(' cluster', ''), '#7f8c8d') for k in keys]
    for i, k in enumerate(keys):
        d = np.asarray(per_group[k]['gmp_cor_draws'])
        a.plot([d.min(), d.max()], [i, i], color=cols[i], lw=1, zorder=1)
        a.scatter(d, np.full_like(d, i) + np.random.default_rng(i).normal(0, .07, d.size),
                  s=9, color=cols[i], alpha=.75, linewidths=0, zorder=2)
    a.set_yticks(range(len(keys)))
    a.set_yticklabels(keys, fontsize=7)
    a.invert_yaxis()
    a.set_xlabel('GMP-Cor, one point per independent cell draw', fontsize=8)
    a.set_title('F  sampling spread at the same n (%d draws)'
                % len(per_group[keys[0]]['gmp_cor_draws']), fontsize=9, loc='left')
    a.tick_params(labelsize=7)

    fig.tight_layout()
    fig.savefig(stem + '.svg')
    fig.savefig(stem + '.png', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    main()
