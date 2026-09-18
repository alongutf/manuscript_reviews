# GMP-Cor: data provenance and depth-confound analysis

**Date:** 2026-07-19, updated 2026-09-18 · **Status:** in progress
**Latest:** section 9 — the rRNA-fraction confound was tested six ways and rejected
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
| `equate_dims_depth_matched` (rRNA removed) | depth matched (~1110 all four), n=287 | dis 190 vs reg 347 (1.8× raw, 2.1× per gene) | weaker control: detection left **2.3× apart** (61 vs 114–141); between-condition gap 157 < within-condition gaps 153/188; dis1 > reg1 |
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

## 8. Open issue — largely resolved, see section 9

**Status (2026-09-18):** the confound described here is real but runs in the
*opposite* direction to the result — see section 9. The scale mismatch below is
also understated: it is ~3×, not ~1.7× (section 9.1).

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

### Next step — superseded

The original plan was to re-select on mRNA counts with one common threshold.
Section 9 shows why that is not the right fix: an mRNA floor conditions directly
on the quantity GMP-Cor measures, trading this bias for a worse one. The
remaining work is instead to fix the **uniform draw** (section 9.7) — select on
a cell-size proxy not derived from the mRNA signal, then equalise `n` and match
depth post hoc by multinomial thinning, at a gene count large enough to avoid
`thinned_match_all`'s 279-gene collapse.

---

## 9. rRNA-fraction selection bias — tested and rejected

**Date:** 2026-09-18 · **Outputs:** `results/cluster_gmp_cor/`

**Question.** Cell selection thresholds `total` counts, which are 70–92% rRNA
depending on the sample (section 8). Does that systematically elevate reg's
mRNA content and depth relative to dis, and so manufacture the reg > dis
GMP-Cor gap?

**Answer: no.** Six tests below. The bias from using total counts runs
*against* reg, not for it. The one mechanism that does inflate reg is the
`equate_dims` uniform draw (section 3.3), which is a pool-size artefact and has
nothing to do with rRNA content.

### 9.1 Composition, all five samples

Full pass over the unfiltered matrices, 100,000 barcodes each. `mRNA` =
protein_coding + pseudogene + ncRNA per the K-12 biotype map; `panel` = every
column except `16s_mature`, `16s_unprocessed`, `mCherry`, `kanR`.

| sample | pool | mean total | mean mRNA | rRNA % | total/mRNA |
|---|---|---|---|---|---|
| exp (2b) | all | 116.4 | 5.0 | **91.9%** | 23.4× |
| exp (2b) | total>400 | 1304.1 | 69.2 | 91.8% | 18.8× |
| dis1 | total>400 | 1330.5 | 81.3 | 91.7% | 16.4× |
| dis2 | total>400 | 856.7 | 44.3 | 92.3% | 19.4× |
| reg1 | total>200 | 563.0 | 78.4 | 73.2% | 7.2× |
| reg2 | total>200 | 682.1 | 110.3 | 70.4% | 6.2× |

exp sits with the dis group on rRNA fraction (91.9%), not with reg. This makes
it a natural control — see 9.6.

**Correction to section 8.** That section estimated the dis/reg selection-scale
mismatch at ~1.7×. Measured directly by asking what total threshold yields the
same cell count as a given mRNA floor, it is closer to **3×**:

| mRNA floor | exp | dis1 | dis2 | reg1 | reg2 |
|---|---|---|---|---|---|
| ≥25 | tot~400 | tot~393 | tot~479 | tot~134 | tot~138 |
| ≥50 | tot~719 | tot~722 | tot~797 | tot~239 | tot~224 |
| ≥100 | tot~1300 | tot~1451 | tot~1474 | tot~487 | tot~446 |

So `umi_min` 400 (dis) / 200 (reg) was a 2× adjustment where ~3× was needed —
i.e. reg was still held to a *stricter* effective mRNA bar than dis.

### 9.2 Total and mRNA counts are not interchangeable

Correlation between total and mRNA per barcode, `total>400` pool:

| sample | Pearson raw | Pearson log1p | Spearman | r² (log) |
|---|---|---|---|---|
| exp | 0.838 | 0.689 | 0.602 | 0.48 |
| dis1 | 0.873 | 0.689 | 0.583 | 0.47 |
| dis2 | 0.846 | 0.781 | 0.616 | 0.61 |
| reg1 | 0.951 | 0.639 | 0.542 | 0.41 |
| reg2 | 0.678 | 0.627 | 0.560 | 0.39 |

The raw Pearson is inflated by a few very deep barcodes spanning orders of
magnitude. On rank, **Spearman is only 0.54–0.62**: total counts explain less
than half the variance in mRNA content. Directly — barcodes with total in
[380, 420], mRNA quantiles:

```
            n     p10   p25   p50   p75   p90
exp       448      5    10    16    24    39
dis1      203      8    12    17    35    48
dis2     5097     12    14    17    21    27
reg1       55      9    29    51    96   179
reg2       64     17    39    60   120   193
```

At an identical total of ~400 a reg1 cell carries 3× the mRNA of an exp cell,
and within reg1 alone the p10–p90 span is 20-fold.

### 9.3 Test 1 — swapping the selection variable changes ~40% of the cells

Top-N by total vs top-N by mRNA, matched n:

| sample | scheme | overlap | overlap % | Jaccard |
|---|---|---|---|---|
| exp | top-1000 | 595 | 59.5% | 0.42 |
| dis1 | top-1000 | 638 | 63.8% | 0.47 |
| dis2 | top-1000 | 665 | 66.5% | 0.50 |
| reg1 | top-1000 | 487 | **48.7%** | 0.32 |
| reg2 | top-1000 | 544 | 54.4% | 0.37 |

35–51% of selected cells change identity. Worst in reg, whose per-cell rRNA
fraction is lowest and most variable.

### 9.4 Test 2 — the bias runs *against* reg

mRNA recovery = mean mRNA of the top-n by total ÷ mean mRNA of the top-n by
mRNA (<1 = a penalty); enrichment = same numerator ÷ a random draw from the
`total>50` pool. Averaged over n = 500/1000/2000/5000:

| condition | mRNA recovery | detection recovery | enrichment vs random |
|---|---|---|---|
| dis | **0.912** | 0.900 | 9.9× |
| exp | 0.843 | 0.836 | 17.6× |
| reg | **0.833** | 0.835 | 5.3× |

Total-based selection costs reg ~17% of available mRNA but dis only ~9% — a 2×
difference in the penalty, in the direction that *suppresses* reg. It is still
4–20× better than a random draw, so it works as a cell-caller; it is just an
inefficient and condition-dependent one.

Mechanism, isolated by holding mRNA fixed (barcodes with mRNA in [40, 60],
split at that band's median total):

| sample | kept mRNA | dropped mRNA | kept rRNA frac | dropped rRNA frac |
|---|---|---|---|---|
| exp | 48.8 | 47.8 | 0.903 | 0.707 |
| dis1 | 49.8 | 47.4 | 0.912 | 0.786 |
| dis2 | 48.8 | 46.8 | 0.916 | 0.860 |
| reg1 | 49.8 | 47.3 | 0.608 | **0.305** |
| reg2 | 49.5 | 48.7 | 0.637 | **0.353** |

At equal mRNA content the criterion does not discriminate on mRNA at all — it
picks the high-rRNA half, and in reg that is a 2× split.

### 9.5 Test 3 — one uniform rule across all five samples (decisive)

Top-1000 by total, nothing else varied:

| sample | cond | rRNA frac | panel depth | mRNA depth | tmRNA share of panel |
|---|---|---|---|---|---|
| exp | exp | 0.922 | 208.4 | 138.4 | 0.204 |
| dis1 | dis | 0.918 | 143.0 | 111.7 | 0.092 |
| dis2 | dis | 0.928 | 329.8 | 250.9 | 0.106 |
| reg1 | reg | 0.733 | 149.1 | 78.5 | 0.367 |
| reg2 | reg | 0.705 | 217.5 | 120.4 | 0.326 |

**Under an identical rule reg is not elevated.** On mRNA depth reg is the
*lowest* of the five (78.5 / 120.4 vs dis 111.7 / 250.9); on panel depth
reg1 ≈ dis1 and dis2 > reg2. Total-count selection does not hand reg a depth
advantage over dis.

### 9.6 Test 4 — exp as a natural control

Across the five samples, rRNA fraction and post-filter panel depth are
unrelated: **Spearman ρ = 0.30, p = 0.62**. exp has dis-like rRNA (92.2%) and
reg-like depth (208, between reg1's 149 and reg2's 218). Final depth is set by
how deeply the library was sequenced, not by its rRNA fraction. With n = 5 the
ρ alone is weak; exp is the load-bearing part of the argument.

### 9.7 Test 5 — what actually inflates reg: the uniform draw

Published-style thresholds, then n = 1000 drawn uniformly from the pool (what
`equate_dims` does) vs taking the top n:

| sample | umi_min | pool | panel depth, top-n | uniform draw | dilution |
|---|---|---|---|---|---|
| exp | 400 | 2862 | 172.5 | 93.1 | 0.54 |
| dis1 | 400 | 1600 | 130.1 | 96.4 | 0.74 |
| dis2 | 400 | 15645 | 258.6 | 56.4 | **0.22** |
| reg1 | 200 | 1004 | 139.6 | 139.5 | **1.00** |
| reg2 | 200 | 1142 | 217.1 | 198.4 | 0.91 |

reg's pool is barely larger than n, so its draw is nearly the whole pool and
loses nothing; dis2's pool is 15.6× n, so the draw reaches far down its depth
distribution and retains 22% of top-n depth. **This asymmetry — not the use of
total counts — is what puts reg above dis in the published comparison.** It is
the section 3.3 defect, and it follows from reg's *small* eligible pool, itself
a consequence of reg's low rRNA content pushing fewer cells past any total
threshold.

### 9.8 Test 6 — tmRNA is not a confound

`DROP_GENES` removes only the 16s probes and the markers, so tmRNA stays in the
panel, and its share differs ~3.5× between conditions (reg1 0.367, reg2 0.326
vs dis1 0.092, dis2 0.106; exp 0.204). The panel-to-mRNA depth ratio is 1.90
for reg1 vs 1.28 for dis1. This was checked separately: **GMP-Cor does not drop
significantly when tmRNA is dropped**, so its share of the panel is not doing
the work. tmRNA is treated as valid biology and retained.

### 9.9 Conclusion

The rRNA-fraction confound hypothesis is rejected on three independent grounds:
under a uniform rule reg is not elevated (9.5); rRNA fraction does not predict
post-filter depth, with exp as the control (9.6); and the criterion's bias
penalises reg roughly twice as hard as dis (9.4). Section 8's open issue is
therefore **largely resolved** — the inconsistent thresholding is real and
worth fixing, but it cannot explain the reg > dis gap because it operates in
the opposite direction.

What remains is the uniform draw (9.7), which is a genuine defect independent
of rRNA. The best control to date, `equate_reg_n1000`, matches detection across
conditions (62.0–67.7 in all four) at comparable n (1000 vs 888) and leaves a
9.3× gap (9.1× per gene) with a between-condition gap 10.8× the largest
within-condition gap, every dis replicate below every reg replicate, and reg
replicates agreeing to 2.7%. Its residual 1.6× depth gap is not enough to close
that: the section 6 diagnostic showed a 1.7× detection cut moved GMP-Cor only
42 → 25.5.

For comparison, `equate_dims_depth_matched` matches depth (~1110 in all four)
but leaves detection **2.3× apart** (dis 61 vs reg 114–141) at n = 287, and
gives only 1.8× raw / 2.1× per gene, with a between-condition gap (157) smaller
than either within-condition gap (153, 188) and dis1 (267) above reg1 (253). It
is the weaker control of the two and should not be quoted as the primary result.

---

## 10. Corrections made during this analysis

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
6. **rRNA fraction does not confound the reg > dis result.** Tested six ways
   (section 9); the bias from selecting on total counts penalises reg ~2× harder
   than dis, and under a uniform rule reg has the *lowest* mRNA depth of the five
   samples. Do not re-raise this without new evidence.
7. **Do not "fix" selection by thresholding mRNA counts.** That conditions on the
   quantity GMP-Cor measures. The defensible design selects on a cell-size proxy
   and matches depth post hoc by thinning.
8. **tmRNA stays in the panel.** Only `16s_mature`, `16s_unprocessed`, `mCherry`
   and `kanR` are dropped. tmRNA is ~33–37% of the reg panel vs ~9–11% of dis,
   but GMP-Cor does not drop significantly without it, and it is valid biology.

---

## 11. Script index

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
| `total_vs_mrna_stats.py` | per-barcode total / mRNA / rRNA, all five samples; rRNA fraction and total-vs-mRNA correlation (9.1, 9.2) |
| `filter_criterion_overlap.py` | selecting on total vs on mRNA — barcode overlap at matched n (9.3) |
| `total_selection_bias.py` | mRNA recovery / enrichment of total-based selection, and its bias at fixed mRNA (9.4) |
| `gene_composition_by_sample.py` | per-gene count composition by biotype; tmRNA share of the panel (9.8) |
| `panel_depth_by_criterion.py` | panel / mRNA depth under one uniform rule across all five samples (9.5) |
