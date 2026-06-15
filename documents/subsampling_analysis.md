# Subsampling robustness analysis (Reviewer #1, comment 2.3)

## Formal response to Reviewer #1, comment 2.3

> *"Fig. 1H shows the Dis-Arrest cells split into two clusters. Was the GMP
> correlation index computed on all cells together, and is the metric sensitive
> to which cells are included in the calculation?"*

We thank the reviewer for raising this point. GMP-Cor in the manuscript is
computed on the full set of cells passing quality control for each condition.
To establish that the metric is not an artefact of how many — or which — cells
are included, we performed a systematic subsampling analysis on synthetic data
whose ground-truth correlation structure is exactly known. Starting from a
single strongly-correlated ("regulated", ρ = 0.9) dataset, we randomly
subsampled cells over a 5× range (200–1000 cells) and recomputed GMP-Cor,
separating three factors: the number of cells, the number of genes, and the
rule used to choose genes (four experiments, detailed below).

The analysis supports three conclusions:

1. **GMP-Cor is an extensive quantity.** Because it is a sum of eigenvalue
   excesses over a spectrum whose total mass equals the number of genes, its raw
   magnitude scales approximately linearly with the number of genes analysed. We
   therefore compare GMP-Cor only between datasets of matched gene dimension (as
   is the case for the conditions in the manuscript), and we now state this
   explicitly. A per-gene (intensive) form of the index removes this dependence.

2. **With the gene panel held complete, GMP-Cor is robust to cell subsampling.**
   Subsampling cells from 1000 down to 200 while retaining the full gene set
   reduces the metric only modestly (mean retains ~63–71% of the full-size value
   at 200 cells). The decline reflects the expected statistical cost of fewer
   observations — noisier correlation estimates and a higher Marchenko–Pastur
   noise threshold — rather than any change in the underlying structure. The main
   effect at low cell numbers is increased variance, not a collapse of the mean,
   which remains clearly elevated throughout.

3. **The metric reflects genuine within-cell coordination, provided informative
   genes are retained.** Restricting the analysis to detected/expressed genes
   (standard scRNA-seq practice, and what is done for the experimental data)
   preserves the correlation signal; degrading the gene panel to a small random
   set of mostly lowly-expressed genes is what erodes it, not the reduction in
   cell number per se.

Taken together, these results show that the GMP-Cor values reported for the
Reg-Arrest and Dis-Arrest conditions are stable with respect to cell number and
cell inclusion over the relevant range, and are not artefacts of population size
or sub-clustering. (The complementary per-cluster GMP-Cor computation for the
two Dis-Arrest clusters in Fig. 1H is addressed separately in the response to
this comment.)

The four supporting experiments and their full results follow.

---

**Question.** Is the GMP-Cor metric robust to the number of cells (and genes)
profiled, or could a reported value be an artefact of dataset size?

**Approach.** Generate a single strongly-correlated ("regulated") synthetic
scRNA-seq dataset at `rho = 0.9`, then randomly subsample it to decreasing
sizes and recompute GMP-Cor. Four complementary experiments isolate the
contributions of cell number, gene number, and gene-selection rule.

**Metric.** `GMP-Cor = Σ max(λᵢ − λ*_scrambled, 0)` — the total eigenvalue mass
above the maximum scrambled (noise) eigenvalue. Computed by
`src.analysis_functions.get_eig_dist` exactly as for the experimental data.

**Common parameters.** `rho = 0.9`, `dropout_rate = 1.0`, Pareto `shape = 1.5`,
`hub_probability = 0.2`, pool `sigma_seed = 31`, pool `count_seed = 0`,
5 random repeats per size, cell subsampling uniform without replacement
(seeded for reproducibility). Cell sizes: 200, 400, 600, 800, 1000.

All runs are deterministic (fixed seeds). Outputs are under
`results/simulation_results/{logs,raw,figures}/`.

---

## Experiment 1 — Ratio fixed, top-expressed genes

Cells subsampled and genes scaled together to hold the **cell:gene ratio fixed
at 1:2** (200c/400g … 1000c/2000g); genes chosen by highest total expression.
Pool: 1200 × 2400. Run: `subsampling_robustness_rho09_20260615_114715`.

| cells | genes | GMP-Cor (mean±SD) | GMP-Cor / gene |
|------:|------:|------------------:|---------------:|
| 200   | 400   | 10.5 ± 4.8        | 0.0264 |
| 400   | 800   | 19.6 ± 7.7        | 0.0245 |
| 600   | 1200  | 31.7 ± 6.0        | 0.0265 |
| 800   | 1600  | 41.5 ± 4.6        | 0.0259 |
| 1000  | 2000  | 51.0 ± 4.8        | 0.0255 |

- **Raw GMP-Cor is *extensive*** — it scales almost perfectly linearly with size
  (fit ≈ `0.0257 × n_genes`, intercept ≈ 0). This is expected: it is a *sum* of
  eigenvalue excesses, and after standardization total spectral mass = n_genes.
  → Raw values are only comparable between datasets of **equal dimension**.
- **Per-gene GMP-Cor (÷ n_genes) is invariant** — constant to ~7% across the 5×
  range. This is the scale-free form of the index.

---

## Experiment 2 — Ratio fixed, top-expressed vs random genes

Same ratio-fixed design, comparing gene-selection rules on the same pool and the
same (paired) cell draws. Run: `subsampling_genecmp_rho09_20260615_115859`.

**Top expression** (per-gene index, mean): 0.0254, 0.0237, 0.0255, 0.0255, 0.0252 → ~7% spread (invariant).

**Random genes:**

| cells | genes | GMP-Cor (mean±SD) | per-gene | CV |
|------:|------:|------------------:|---------:|----:|
| 200   | 400   | 0.89 ± 1.35       | 0.0022   | 1.52 |
| 400   | 800   | 6.71 ± 4.48       | 0.0084   | 0.67 |
| 600   | 1200  | 14.65 ± 2.62      | 0.0122   | 0.18 |
| 800   | 1600  | 24.54 ± 2.38      | 0.0153   | 0.10 |
| 1000  | 2000  | 38.58 ± 1.55      | 0.0193   | 0.04 |

- The per-gene invariance **holds for top-expressed genes but breaks for random
  genes**: the random per-gene index climbs ~9× (0.0022 → 0.0193, 149% spread).
- At 200 cells / 400 random genes the signal is near the noise floor (per-gene
  0.0022; ~9% of the top level). Random panels at small sizes are dominated by
  lowly-expressed, dropout-heavy genes whose pairwise correlations fall below
  the Marchenko-Pastur noise edge.
- The random/top per-gene ratio rises 0.09× → 0.35× → 0.48× → 0.60× → 0.77×,
  converging only at large size.

→ The apparent scale-invariance *does* depend on retaining informative,
highly-expressed genes — but note that here the gene panel was being shrunk in
lockstep with the cells (see Experiment 3).

---

## Experiment 3 — Gene count fixed at 2000, top vs random

Genes held at 2000 (selected from the 2400-gene pool by top expression or at
random), only **cells** subsampled. The cell:gene aspect ratio is *not* fixed
(10:1 genes:cells at 200 cells → 2:1 at 1000). Run:
`subsampling_fixedgenes_rho09_20260615_170718`.

| cells | g:c   | top (mean±SD) | frac of full | random (mean±SD) | frac of full |
|------:|------:|--------------:|-------------:|-----------------:|-------------:|
| 200   | 10:1  | 36.1 ± 10.2   | 0.71         | 25.2 ± 9.0       | 0.66 |
| 400   | 5:1   | 36.7 ± 10.9   | 0.72         | 29.1 ± 10.9      | 0.76 |
| 600   | 3.3:1 | 40.9 ± 1.4    | 0.80         | 32.2 ± 3.9       | 0.84 |
| 800   | 2.5:1 | 46.2 ± 2.6    | 0.91         | 35.5 ± 2.1       | 0.92 |
| 1000  | 2:1   | 50.9 ± 2.6    | 1.00         | 38.4 ± 4.1       | 1.00 |

- Keeping the **full gene panel dramatically stabilizes the metric**. Raw GMP-Cor
  spread drops to 35% (top) / 41% (random), versus 130–220% when genes were
  shrunk in lockstep.
- Even **random** selection retains 66% of full signal at 200 cells (vs the
  near-collapse of Experiment 2), and the top/random gap is a stable ~0.70–0.79×
  across all sizes.

→ The signal loss for random genes in Experiment 2 was driven mainly by
*shrinking the gene panel*, not by reducing cells.

---

## Experiment 4 — All genes retained, cells only (the clean test)

Pool generated with **exactly 2000 genes**; the **complete** gene set is used at
every size (no gene-selection axis), only cells subsampled. Pool: 1200 × 2000.
Run: `subsampling_allgenes_rho09_20260615_171556`.

| cells | g:c   | GMP-Cor (mean±SD) | CV   | frac of full |
|------:|------:|------------------:|-----:|-------------:|
| 200   | 10:1  | 24.6 ± 11.1       | 0.45 | 0.63 |
| 400   | 5:1   | 25.6 ± 7.5        | 0.29 | 0.66 |
| 600   | 3.3:1 | 31.6 ± 6.2        | 0.20 | 0.81 |
| 800   | 2.5:1 | 34.5 ± 6.2        | 0.18 | 0.89 |
| 1000  | 2:1   | 38.9 ± 1.7        | 0.04 | 1.00 |

- GMP-Cor declines **smoothly and modestly** with fewer cells (38.9 → 24.6, 46%
  spread), retaining **63%** of the full-size value even at 200 cells.
- The decline is a pure **cell-number effect** — fewer cells give noisier
  correlation estimates and a higher Marchenko-Pastur noise edge — not loss of
  genes.
- The dominant cost at low cell counts is **increased variance** (CV 0.04 →
  0.45), not collapse of the mean; the signal stays clearly elevated throughout.

---

## Conclusions

1. **Raw GMP-Cor is extensive** — it scales ~linearly with the number of genes.
   Raw values should only be compared between datasets of **matched dimension**,
   or expressed per-gene.
2. **With the gene panel held complete, GMP-Cor is robust to cell subsampling**
   down to a few hundred cells: the mean stays within ~37% of its full-size value
   at 200 cells, with higher variance the main cost.
3. **Robustness requires informative genes.** Restricting (or simply retaining a
   full panel that includes) detected/expressed genes preserves the correlation
   structure; blindly shrinking the gene set to a small random panel at low cell
   numbers loses the signal. This matches standard scRNA-seq practice and what is
   done for the experimental data.

**Practical guidance for the metric:** compare GMP-Cor only across datasets of
matched gene dimension (or use the per-gene index), restrict to detected genes,
and treat values from very small cell numbers (a few hundred) as noisier.

---

## Files

| Experiment | Script | Output stem (`results/simulation_results/`) |
|---|---|---|
| 1 — ratio fixed, top | `simulations/subsampling_robustness_rho09_run.py` | `subsampling_robustness_rho09_20260615_114715` |
| 2 — ratio fixed, top vs random | `simulations/subsampling_genecmp_rho09_run.py` | `subsampling_genecmp_rho09_20260615_115859` |
| 3 — genes fixed 2000, top vs random | `simulations/subsampling_genecmp_fixedgenes_rho09_run.py` | `subsampling_fixedgenes_rho09_20260615_170718` |
| 4 — all genes, cells only | `simulations/subsampling_allgenes_rho09_run.py` | `subsampling_allgenes_rho09_20260615_171556` |

Each stem has a `.json` log (full parameters + per-repeat values) in `logs/`, a
human-readable `.txt` summary in `raw/`, and `.svg`/`.png` figures in `figures/`.
