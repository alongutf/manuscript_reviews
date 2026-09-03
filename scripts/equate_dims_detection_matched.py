"""equate_dims matched on genes-detected, within replicate pairs only.

Matching total counts across all four samples left reg detecting ~2x more genes
than dis at equal depth -- and detection, not raw depth, is what the eigenvalue
spectrum responds to. So here:

  * matching variable: genes detected per cell on the sample's final panel
  * matching scope:    within each replicate pair only (dis1<->dis2,
                       reg1<->reg2), not across conditions

Everything else follows the paper pipeline: cell pool from
filter_by_umi_count(400, 20000) on totals, gene panel from
filter_by_gene_dispersion(min_dispersion=1) intersected within the pair,
target ~1000 cells.

Depth (total counts) is reported but NOT matched, so the residual depth
difference after detection-matching is visible.

Eligible pools are cached so reruns skip the full CSV passes.
"""
import os
import sys
import json
import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import src.analysis_functions as af

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# NOTE: absolute path to a sibling repo on this machine holding the raw
# unfiltered probe-count CSVs; must be updated to run on another machine.
OTHER = r'C:\Users\owner\Documents\Projects\rnaseq_correlations'
DATA = os.path.join(OTHER, 'data')
PAPER = os.path.join(ROOT, 'data_for_paper')
OUTDIR = os.path.join(ROOT, 'results', 'cluster_gmp_cor')
POOLDIR = os.path.join(OUTDIR, 'eligible_pools')

UNFILT = {
    'dis1': 'sample_13a_unfiltered.csv',
    'dis2': 'sample_15a_unfiltered.csv',
    'reg1': 'sample_13b_unfiltered.csv',
    'reg2': 'sample_15b_unfiltered.csv',
}
PUBLISHED = {
    'dis1': 'sample_13a_filtered.csv',
    'dis2': 'sample_15a_filtered.csv',
    'reg1': 'sample_13b_filtered.csv',
    'reg2': 'sample_15b_filtered.csv',
}
PAIRS = [('dis1', 'dis2'), ('reg1', 'reg2')]
ORDER = ['dis1', 'dis2', 'reg1', 'reg2']

UMI_MIN, UMI_MAX = 400, 20000  # total-count window for the eligible cell pool (paper pipeline)
MIN_DISPERSION = 1.0           # AnnMat.filter_by_gene_dispersion threshold, same as the paper
# Dropped from the final panel, matching the published matrices. Cell selection
# still uses TOTAL counts (rRNA included) -- that is what reproduces the paper's
# barcode set -- but these genes are excluded from the matrix itself.
# 16s_mature alone is ~92% of counts, so leaving it in dominates the spectrum.
DROP_GENES = ['16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR', 'mCherry']
TARGET_CELLS = 1000  # cells requested per matched draw (actual count may be lower, see matched_draw)
N_BINS = 40           # quantile bins used to match the genes-detected distribution
REPS = 5              # independent matched draws per replicate pair, for a mean +/- spread
CHUNK = 20000         # rows per pandas read_csv chunk when streaming the large unfiltered CSVs
NORM_SUM = 50         # row-sum target passed to af.get_eig_dist's normalize step
SEED = 0              # base RNG seed; rep k uses SEED + rep, so reruns are reproducible per rep


def load_eligible(sample):
    # Build (or load from a cached .npz) the pool of cells passing the UMI-count
    # window for one sample: read the raw unfiltered probe-count CSV in chunks
    # (it is too large to load at once), keep only rows with UMI_MIN < total <
    # UMI_MAX, and stack them into one dense matrix. Caching lets reruns of this
    # script (e.g. with different REPS or matching parameters) skip the slow
    # full-CSV pass entirely.
    cache = os.path.join(POOLDIR, f'{sample}.npz')
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True)
        print(f'{sample}: cache hit ({d["M"].shape})')
        return d['M'].astype(float), d['bc'], d['genes']
    path = os.path.join(DATA, UNFILT[sample])
    header = pd.read_csv(path, nrows=0)
    # int32 counts keep memory down; barcode column stays str
    dtypes = {c: np.int32 for c in header.columns[1:]}
    dtypes[header.columns[0]] = str
    genes = np.asarray(header.columns[1:])
    rows, bcs = [], []
    for ch in pd.read_csv(path, index_col=0, chunksize=CHUNK, dtype=dtypes):
        V = ch.values
        tot = V.sum(1)
        sel = (tot > UMI_MIN) & (tot < UMI_MAX)
        if sel.any():
            rows.append(V[sel])
            # strip the GEM well suffix ("-1" etc.) so barcodes match across files
            bcs.append(np.array([str(b).split('-')[0] for b in ch.index])[sel])
    M = np.vstack(rows)
    bc = np.concatenate(bcs)
    os.makedirs(POOLDIR, exist_ok=True)
    np.savez_compressed(cache, M=M, bc=bc, genes=genes)
    print(f'{sample}: cached {M.shape}')
    return M.astype(float), bc, genes


def dispersion_genes(M, genes):
    # Genes passing AnnMat.filter_by_gene_dispersion's default threshold
    # (variance/mean > MIN_DISPERSION), computed over the eligible cell pool M.
    mean = M.mean(axis=0)
    var = M.var(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        disp = np.where(mean > 0, var / mean, 0.0)
    return set(genes[np.nan_to_num(disp) > MIN_DISPERSION])


def matched_draw(vals, target, n_bins, rng):
    """Draw `target` indices per sample sharing one histogram of `vals`.

    `vals` maps sample name -> 1-D array (here, genes-detected per cell). All
    samples' values are pooled to define shared quantile bin edges, so every
    sample is binned on the same scale. Within each bin, the number of cells
    drawn is capped by whichever sample has the fewest cells in that bin
    (`avail`), so the returned per-sample subsets share one histogram shape --
    this is what "matching" the genes-detected distribution means here.
    If the total available cells across bins is below `target`, all of them are
    taken (`take = avail`) and the actual yield is reported as `n_got` by the
    caller; otherwise bins are proportionally down-sampled to hit `target` with
    a greedy remainder pass so the sum matches exactly.

    Returns (picks, n_drawn) where picks maps sample -> array of row indices
    into that sample's own matrix, and n_drawn is the total cells taken.
    """
    allv = np.concatenate(list(vals.values()))
    edges = np.unique(np.quantile(allv, np.linspace(0, 1, n_bins + 1)))
    binned = {s: np.digitize(v, edges[1:-1]) for s, v in vals.items()}
    nb = len(edges) - 1
    avail = np.array([min(int((binned[s] == b).sum()) for s in vals)
                      for b in range(nb)])
    total = int(avail.sum())
    if total == 0:
        raise RuntimeError('no overlap')
    # proportionally allocate `target` draws across bins, then top up bin-by-bin
    # (largest remaining headroom first) until the total matches exactly
    take = avail if total <= target else np.floor(avail * target / total).astype(int)
    while take.sum() < min(target, total):
        take[int(np.argmax(avail - take))] += 1
    picks = {}
    for s in vals:
        out = []
        for b in range(nb):
            idx = np.where(binned[s] == b)[0]
            k = int(take[b])
            if k > 0:
                out.append(rng.choice(idx, k, replace=False))
        picks[s] = np.concatenate(out) if out else np.array([], dtype=int)
    return picks, int(take.sum())


def gmp_cor(m):
    # GMP-Cor = sum of eigenvalue excess above the scrambled-null maximum
    # (lambda_i - lambda_max_scrambled, clipped at 0), summed over all genes.
    pcs, pcs1, _ = af.get_eig_dist(m, norm=True, log=False,
                                   norm_method='sum', norm_sum=NORM_SUM)
    thr = float(pcs1.max())
    return (float(np.sum(np.maximum(pcs - thr, 0))), float(pcs.max()), thr)


def main():
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    # ---- load eligible cell pools (UMI-count filtered) for every sample ----
    pool = {}
    for s in ORDER:
        M, bc, genes = load_eligible(s)
        pool[s] = {'M': M, 'bc': bc, 'genes': genes,
                   'disp': dispersion_genes(M, genes)}

    # ---- gene panel: dispersion-filtered genes shared within each replicate
    # pair (dis1&dis2, reg1&reg2), NOT across the two pairs. The dis panel and
    # the reg panel can therefore end up with different gene counts (see
    # n_genes printed per sample below) -- GMP-Cor is dimension-dependent
    # (it scales with p, the number of genes/kept modes), so an unequal panel
    # size between the two conditions is a potential confound for any dis-vs-
    # reg comparison, even though depth/detection are matched within each pair.
    panel = {}
    for a, b in PAIRS:
        common = sorted((pool[a]['disp'] & pool[b]['disp']) - set(DROP_GENES))
        panel[a] = panel[b] = common
        print(f'panel {a}+{b}: {len(common)} genes (after dropping {DROP_GENES})')

    # ---- subset each sample's pool to its panel, and record its detected-genes
    # and total-depth distributions (used by matched_draw and for QC printing)
    sub, det, dep = {}, {}, {}
    for s in ORDER:
        ix = pd.Index(pool[s]['genes']).get_indexer(panel[s])
        sub[s] = pool[s]['M'][:, ix]
        det[s] = (sub[s] > 0).sum(1)
        dep[s] = sub[s].sum(1)
        print(f'{s}: pool n={len(det[s])}  detected median {np.median(det[s]):6.1f}  '
              f'depth median {np.median(dep[s]):8.1f}')

    # barcodes of the published (paper) filtered matrices, for the "how much of
    # the originally-published cell set does this matched draw retain" check
    published = {}
    for s in ORDER:
        idx = pd.read_csv(os.path.join(PAPER, PUBLISHED[s]),
                          index_col=0, usecols=[0]).index
        published[s] = set(str(i).split('-')[0] for i in idx)

    # ---- repeated matched draws + GMP-Cor -----------------------------------
    records, overlaps = [], []
    for rep in range(REPS):
        rng = np.random.default_rng(SEED + rep)
        for a, b in PAIRS:
            picks, n_got = matched_draw({a: det[a], b: det[b]},
                                        TARGET_CELLS, N_BINS, rng)
            if rep == 0:
                print(f'\npair {a}+{b}: matched draw yields {n_got} cells each')
            for s in (a, b):
                ix = picks[s]
                m = sub[s][ix]
                bcs = set(pool[s]['bc'][ix])
                g, lam, thr = gmp_cor(m)
                records.append({
                    'sample': s, 'pair': f'{a}+{b}', 'rep': rep,
                    'n_cells': int(m.shape[0]), 'n_genes': len(panel[s]),
                    'mean_genes_detected': float((m > 0).sum(1).mean()),
                    'median_genes_detected': float(np.median((m > 0).sum(1))),
                    'mean_depth': float(m.sum(1).mean()),
                    'median_depth': float(np.median(m.sum(1))),
                    'lambda_max': lam, 'lambda_max_scrambled': thr, 'gmp_cor': g,
                })
                overlaps.append({
                    'sample': s, 'rep': rep, 'n_new': len(bcs),
                    'shared_with_published': len(bcs & published[s]),
                    'pct_of_published_retained':
                        round(100 * len(bcs & published[s]) / len(published[s]), 1),
                })
                if rep == 0:
                    r = records[-1]
                    print(f'  {s}: detected {r["mean_genes_detected"]:.1f}  '
                          f'depth {r["mean_depth"]:.1f}  GMP-Cor {g:.3f}  '
                          f'| {overlaps[-1]["pct_of_published_retained"]}% of '
                          f'published barcodes')

    df = pd.DataFrame(records)
    ov = pd.DataFrame(overlaps)

    # ---- summary printouts ---------------------------------------------------
    print('\n=== detection (matched) vs depth (not matched) ===')
    summ = df.groupby('sample')[['n_cells', 'mean_genes_detected',
                                 'mean_depth']].agg(['mean', 'std']).loc[ORDER]
    print(summ.to_string(float_format=lambda v: '%.2f' % v))

    print('\n=== GMP-Cor ===')
    g = df.groupby('sample')['gmp_cor'].agg(['mean', 'std', 'min', 'max']).loc[ORDER]
    print(g.to_string(float_format=lambda v: '%.3f' % v))

    print('\n=== barcode change vs data_for_paper ===')
    print(ov.groupby('sample')[['shared_with_published',
                                'pct_of_published_retained']]
            .mean().loc[ORDER].to_string(float_format=lambda v: '%.1f' % v))

    d1, d2 = g.loc['dis1', 'mean'], g.loc['dis2', 'mean']
    r1, r2 = g.loc['reg1', 'mean'], g.loc['reg2', 'mean']
    print(f'\n  dis mean {np.mean([d1,d2]):.3f} | reg mean {np.mean([r1,r2]):.3f}')
    print(f'  within-dis gap {abs(d1-d2):.3f} | within-reg gap {abs(r1-r2):.3f} '
          f'| between-condition gap {abs(np.mean([d1,d2])-np.mean([r1,r2])):.3f}')

    # ---- write outputs ---------------------------------------------------
    df.to_csv(os.path.join(OUTDIR, f'equate_detection_matched_{stamp}.csv'), index=False)
    ov.to_csv(os.path.join(OUTDIR, f'equate_detection_overlap_{stamp}.csv'), index=False)
    with open(os.path.join(OUTDIR, f'equate_detection_matched_{stamp}.json'), 'w') as fh:
        json.dump({'params': {'umi_min': UMI_MIN, 'umi_max': UMI_MAX,
                              'min_dispersion': MIN_DISPERSION,
                              'target_cells': TARGET_CELLS, 'n_bins': N_BINS,
                              'reps': REPS, 'norm_sum': NORM_SUM, 'seed': SEED,
                              'match_variable': 'genes_detected',
                              'match_scope': 'within replicate pair'},
                   'results': records, 'overlaps': overlaps}, fh, indent=2)
    print(f'\nwrote results/cluster_gmp_cor/equate_detection_matched_{stamp}.*')


if __name__ == '__main__':
    main()
