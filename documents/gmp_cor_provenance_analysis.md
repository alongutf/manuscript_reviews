# GMP-Cor: data provenance and depth-confound analysis

**Date:** 2026-07-19 · **Status:** in progress, to be continued
**Scripts:** `scripts/` · **Outputs:** `results/cluster_gmp_cor/`

---

## 1. What this started as

Original request: get leiden clusters from `data/scanpy_shx.h5ad`, use them to
filter the five SHX samples in `data_for_paper/`, and compute GMP-Cor per
sample × cluster.

That surfaced a discrepancy — the barcodes in `data_for_paper/` and the h5ad
barely overlapped — which turned into an investigation of how the published
matrices were generated and whether the published GMP-Cor comparisons are
confounded. That investigation is the substance of this document.

**GMP-Cor definition used throughout:** `sum_i max(lambda_i - lambda_max_scrambled, 0)`
via `af.get_eig_dist(m, norm=True, log=False, norm_method='sum', norm_sum=50)`,
with the scrambled threshold averaged over 10 permutations (`rep=10`, hardcoded
in `analysis_functions.get_eig_dist`).

---

## 2. Sample naming

| h5ad batch | source matrix | condition |
|---|---|---|
| `exp`  | `sample_2b`  | exponential |
| `dis1` | `sample_13a` | dysregulated |
| `dis2` | `sample_15a` | dysregulated |
| `reg1` | `sample_13b` | regulated |
| `reg2` | `sample_15b` | regulated |

Established from `scripts/scanpy_analysis.ipynb` cell 0 (concat order
`[exp, shx1, shx2, casp1, casp2]` with keys `['exp','dis1','dis2','reg1','reg2']`)
and confirmed by barcode overlap.

---

## 3. Provenance of the published matrices — solved

### 3.1 `data_for_umap/` (and therefore the h5ad)

Reproduced at **96%**: top-1000 barcodes by **total** count (rRNA included),
i.e. `filter_by_umi_count` with `target_cells=1000` on the unfiltered matrix,
no biotype filter first. Since `16s_mature` is ~92% of counts in these samples,
this is effectively *"the 1000 highest-rRNA cells."*

### 3.2 `data_for_paper/` — reproduced exactly

`data_for_paper/sample_15a_filtered.csv` was reproduced **1000/1000 barcodes**:

```
filter_by_umi_count(400, 20000)  on TOTAL counts   -> 15,645 cells for 15a
                                                      1,600 cells for 13a
filter_by_gene_dispersion(min_dispersion=1)
equate_dims(amat_13a, amat_15a, 1000):
    np.random.seed(0)
    draw 1 -> 1000 of  1,600  (13a)
    draw 2 -> 1000 of 15,645  (15a)   <- matches published file exactly
```

This is `analysis_notebook.ipynb` **cell 3**, not the cells 4→6 path. Chance
overlap would be ~64 barcodes; we got 1000/1000.

Gene panel: the pairwise **intersection** of the two samples' dispersion-passing
genes, minus `['16s_mature','16s_unprocessed','LELOBEKK','kanR','mCherry']`.
Reproduced 2042/2042 genes for the dis pair — exact.

This explains the paired gene counts: 13a/15a both 2042, 13b/15b both 2023,
2b 2071 (paired differently).

### 3.3 The defect

`equate_dims` draws **uniformly** from each pool. The pools differ enormously
(13a: 1,600 cells; 15a: 15,645), so the same 1000-cell draw reaches far deeper
down 15a's depth distribution. Result: **the deeper sample ends up with the
shallower matrix** — published mean total 90 (13a) vs 50 (15a), inverting the
raw-data relationship.

Since GMP-Cor is strongly depth-dependent, any 13a-vs-15a comparison built from
these files carries a depth confound running *opposite* to the biology.

---

## 4. Why GMP-Cor values differ so much between sources

Three mechanisms, all verified:

1. **GMP-Cor is extensive in `p`.** After z-transforming columns, trace = p, so
   the total spectral mass *is* the gene count. Verified: trace = 3874/3859/3895
   /3852/3999 for h5ad matrices, 2071/2037/2004/2023/2023 for paper matrices —
   exactly p in every case. Doubling genes doubles GMP-Cor. **It is not
   comparable across gene sets of different size.** Report `GMP-Cor / p`, or
   always state p.

2. **The noise threshold is pure matrix shape.** `lambda_max_scrambled` is the
   Marchenko–Pastur edge `(sqrt(gamma)+1)^2`. h5ad gamma≈3.9 → predicted 8.8,
   observed 9.4–10.0. Paper gamma≈2.0 → predicted 5.9, observed 6.1–6.5. Carries
   no biological information.

3. **Detection sparsity destroys signal.** Ordering the paper matrices by genes
   detected per cell: exp 101 → 33.5, reg2 87 → 29.1, dis1 63 → 2.9, reg1 59 →
   22.5, dis2 39 → **0.00**. For `dis2:paper`, lambda_max (6.095) fell *below*
   its own scrambled edge (6.184) — no mode above noise at all.

---

## 5. Raw data is much thinner than the published files imply

Depth/cell tradeoff on mRNA counts (rRNA+tRNA removed), 100,000 barcodes/sample:

| T (mRNA) | exp | dis1 | dis2 | reg1 | reg2 | min | usable genes |
|---|---|---|---|---|---|---|---|
| 50  | 1241 | 743 | 2723 | 769 | 957 | 743 | 148 |
| 75  |  759 | 453 | 1449 | 386 | 567 | **386** | **279** |
| 100 |  532 | 310 |  960 | 229 | 383 | 229 | 392 |
| 150 |  321 | 193 |  536 | 115 | 207 | 115 | 630 |

Median mRNA per barcode is **1–12 counts**; most of the 100k barcodes are empty
droplets. The published matrices contain 1000 cells each, which exceeds the
number of real cells at usable depth in reg1 (229 at T=100) and dis1 (310).

`dis2` is a quality outlier: median mRNA 12 vs 1–3 elsewhere, 4.5× more cells at
every threshold. Same nominal condition as dis1. **Unexplained — worth checking
whether this is a library-prep or sequencing-depth difference.**

---

## 6. Conclusions on the dis-vs-reg question

The answer flipped twice as confounds were removed. Chronology matters because
earlier runs are superseded:

| run | control | result | verdict |
|---|---|---|---|
| published matrices | none | reg >> dis | confounded (gene sets, depth) |
| `thinned_match_all` | depth+detection matched **across** conditions, 279 genes, n=386 | dis 3.71 vs reg 3.02 — no separation | **likely over-controlled**; collapsing to 279 genes plausibly destroyed signal |
| `equate_dims_detection_matched` (rRNA present — **discard**) | detection matched within pairs | replicate gaps ≥ condition gap | invalid: rRNA not removed |
| `equate_dims_detection_matched` (rRNA removed) | detection matched within pairs | dis 2.64 vs reg 39.44 (15×) | reg 3× deeper — confounded |
| **`equate_reg_n1000` (current best)** | reg `umi_min` lowered to 200 → detection matched across conditions | **dis 2.70 vs reg 25.15 (9.3×)** | see below |

### Current best result (`reg_n1000_20260719_185739`)

| sample | n | genes | detected | depth | GMP-Cor |
|---|---|---|---|---|---|
| dis1 | 1000 | 2042 | 62.6 ± 0.4 | 88.6 ± 1.0 | 3.74 ± 0.54 |
| dis2 | 1000 | 2042 | 62.0 ± 0.6 | 81.1 ± 1.5 | 1.66 ± 1.12 |
| reg1 |  888 | 2100 | 66.5 ± 0.0 | 138.3 ± 0.2 | 25.50 ± 2.76 |
| reg2 |  888 | 2100 | 67.7 ± 1.5 | 127.7 ± 4.3 | 24.81 ± 1.89 |

```
dis mean 2.70 | reg mean 25.15          -> 9.3x
within-dis gap 2.08 | within-reg gap 0.69 | between-condition gap 22.45
```

**The key diagnostic:** lowering reg's `umi_min` pulled detection 113 → 66.5
(onto dis's 62) and depth 262 → 138, but GMP-Cor only fell 42 → 25.5. A 1.7×
detection reduction gave a 1.65× GMP-Cor reduction while the reg/dis ratio
barely moved (15× → 9.3×). **If detection drove the condition difference,
matching it should have collapsed the gap. It did not.** The condition gap is
10× the largest replicate gap, and reg replicates agree to 2.7%.

Direction (reg > dis) matches the manuscript's prediction.

**Incidental finding:** `umi_min=200` is probably what the paper actually used
for the reg pair — barcode retention jumped to reg1 88.4% / reg2 77.9% (from
28%/25%) and the panel moved to 2100 genes vs the published 2023.

---

## 7. Exact parameters — current best result

Script: `scripts/equate_reg_n1000.py`
Outputs: `results/cluster_gmp_cor/reg_n1000_20260719_185739.{csv,json}`,
`reg_n1000_overlap_20260719_185739.csv`

```json
{
  "umi_floor": 50,
  "reg_umi_min": 200,
  "dis_umi_min": 400,
  "umi_max": 20000,
  "min_dispersion": 1.0,
  "target_cells": 1000,
  "reps": 5,
  "drop_genes": ["16s_mature", "16s_unprocessed", "LELOBEKK", "kanR", "mCherry"],
  "norm_sum": 50,
  "seed": 0
}
```

Additional settings fixed in the script:

```
N_BINS        = 40      # quantile bins for the detection-matched draw
matching var  = genes detected per cell, on the sample's final panel
matching scope= within replicate pair only (dis1<->dis2, reg1<->reg2)
cell pool     = UMI_FLOOR < total < UMI_MAX, total over ALL 4184 genes (rRNA in)
gene panel    = pairwise intersection of dispersion-passing genes, minus DROP_GENES
GMP-Cor       = af.get_eig_dist(norm=True, log=False, norm_method='sum', norm_sum=50)
scramble reps = 10 (hardcoded in analysis_functions.get_eig_dist)
RNG           = np.random.default_rng(SEED + rep), rep in 0..4
```

Cached inputs (avoid re-scanning the 840 MB CSVs):
- `results/cluster_gmp_cor/eligible_pools/{dis1,dis2}.npz` — pools at total>400
- `results/cluster_gmp_cor/eligible_pools/{reg1,reg2}_low.npz` — pools at total>50
- `results/cluster_gmp_cor/barcode_stats/{sample}.npz` — per-barcode mRNA depth
- Unfiltered source CSVs: `C:\Users\owner\Documents\Projects\rnaseq_correlations\data\`

---

## 8. Open issue — must resolve before publishing this result

**Depth is still measured and thresholded inconsistently across conditions.**

Two different quantities are in play:
- **selection depth** = row sum over all 4184 genes, rRNA included (what
  `umi_min` thresholds)
- **reported depth** = row sum over the final panel, rRNA excluded (the
  `mean_depth` column)

They differ by a factor that is **not constant across conditions**, because
rRNA fraction differs:

| sample | total (all genes) | panel sum | rRNA % | ratio |
|---|---|---|---|---|
| dis1 | 1218.3 |  87.4 | 91.6% | 13.9× |
| dis2 |  760.3 |  50.4 | 92.0% | 15.1× |
| reg1 |  525.3 | 129.1 | 73.1% |  4.1× |
| reg2 |  664.1 | 176.3 | 69.6% |  3.8× |

So `umi_min=400` for dis and `200` for reg was **not** a neutral adjustment: a
dis cell at total=400 carries ~29 mRNA counts, a reg cell at total=200 carries
~50. The two conditions were selected on scales differing ~1.7× in the quantity
that matters. This is the same class of error as the original `equate_dims` bug
— thresholding on a number dominated by rRNA whose fraction varies between the
samples being compared.

Also open:
- Residual **1.6× depth gap** (reg 138/128 vs dis 89/81) in panel units.
- dis and reg panels are **different gene sets** (2042 vs 2100), so row sums are
  not over identical genes — fine within a pair, loose across conditions.
- `n` is 1000 (dis) vs 888 (reg).

### Next step

Re-run selecting cells on **panel/mRNA counts directly**, with the *same*
threshold for all four samples, so the criterion means the same thing regardless
of a sample's rRNA content. Then equalise `n` and thin reg to dis's exact panel
depth. If the ~9× gap survives, the result is solid.

---

## 9. Corrections made during this analysis

Recorded so they are not re-introduced:

1. **`equate_dims` is not dead code.** It generated the published matrices. It
   is simply not called from `figure5.py`, which only reads its output.
2. **The published p values are not "matched by accident"** — they are matched
   deliberately, pairwise, by `equate_dims`.
3. **rRNA removal was initially omitted** from the `equate_dims_*` reconstruction
   scripts, making depths ~14× too high and GMP-Cor ~100× too high. All results
   from those runs before the fix are invalid.
4. **Selecting on mRNA vs total was not the cause** of the paper/umap barcode
   divergence (that hypothesis was tested and refuted — the mRNA-based
   reconstruction overlapped the *umap* file at 66%). The actual cause is
   top-N-by-total vs uniform-random-from-a-large-pool.
5. **The "no dis-vs-reg difference" conclusion was premature** — it rested on
   rRNA-contaminated runs and on an over-controlled 279-gene comparison.

---

## 10. Script index

| script | purpose |
|---|---|
| `cluster_gmp_cor.py` | GMP-Cor per sample × leiden cluster, cells from `data_for_paper` |
| `cluster_gmp_cor_h5ad.py` | same, sourced entirely from the h5ad |
| `sample_geneset_gmp_cor.py` | 2000-gene selections (max/mean/Fano) per sample |
| `source_comparison_diagnostics.py` | h5ad vs data_for_paper structural + spectral stats |
| `matched_downsampling.py` | first matched-p/n/depth attempt, both sources |
| `reproduce_filter_paths.py` | reproduces the two notebook filter paths |
| `locate_paper_cells.py` | profiles published barcodes inside the unfiltered matrix |
| `depth_tradeoff_curve.py` | cells-vs-depth-vs-genes feasibility curve |
| `capped_match_13a_15a.py` | cap 15a at 13a's max, top-1000 below |
| `thinned_match_13a_15a.py` | exact depth match via multinomial thinning, one pair |
| `thinned_match_all.py` | same across all five samples |
| `equate_dims_depth_matched.py` | equate_dims with depth-histogram matching |
| `equate_dims_detection_matched.py` | equate_dims with detection matching, within pairs |
| `equate_reg_n1000.py` | **current best** — reg raised to n≈1000 via umi_min=200 |
