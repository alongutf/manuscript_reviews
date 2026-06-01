# Response Plan: Reviewer #1 Comments

## Overview

Reviewer #1 has two broad concerns: (1) the biological interpretation and experimental design, and (2) the GMP-Cor metric validation. The metric-validation concerns (Comment 2) are the most tractable computationally and map directly onto the existing simulations. Comments 1, 3, 4, and 5 require either re-analysis of existing data or targeted experimental work.

**Updated GMP-Cor definition (replaces GMP model fit):** The GMP model fit to the Generalized Marchenko-Pastur distribution was found to be numerically unreliable — fits are not robust or reproducible across datasets. GMP-Cor is now defined as:

> **GMP-Cor = Σ max(λᵢ − λ\*\_scrambled, 0)** for all eigenvalues λᵢ
>
> where λ\*\_scrambled is the maximum eigenvalue of the scrambled (gene-permuted) correlation matrix. Only positive differences are summed — this counts the total excess of observed eigenvalues above the noise threshold set by scrambling.

Implementation: `get_eig_dist()` in `src/analysis_functions.py` already returns `pcs` (observed eigenvalues) and `pcs1` (scrambled eigenvalues). The new scalar is `sum(pcs[pcs > pcs1.max()] - pcs1.max())`. Figure 2 and its caption will need updating to reflect this change.

**Simulation approach:** All response simulations use the synthetic scRNA-seq pipeline from `simulations/simulated_data.ipynb` (`simulate_scRNA_data()` function). The dynamical ODE/SDE simulations (`main.py`, `competitive_binding.py`, `subpopulations.py`) are unreliable for this purpose — the correlation signal in their output is dominated by high data dispersion rather than true gene-gene couplings. The synthetic approach generates correlated Gaussian latent variables and maps them through Negative Binomial inverse-transform sampling to produce realistic count data with a known, controlled correlation structure. The `rho` parameter of `generate_gram_hub_matrix()` is the clean control knob: high rho → regulated, low rho → dysregulated.

---

## Comment 1 — Why scRNA-seq? Subpopulations vs. Uniform Dysregulation

**What the reviewer is asking:** Could the observed loss of gene-gene correlations arise from a mixture of two regulated subpopulations rather than uniform dysregulation? And why is scRNA-seq the right tool here?

**Proposed response steps:**

1. **Simulation argument (computational):** Use `simulate_scRNA_data()` to generate two contrasting scenarios:
   - **Scenario A (dysregulation):** Single population, low rho (e.g., rho=0.1) → genuinely low GMP-Cor
   - **Scenario B (subpopulations):** Two sub-populations, each with high rho (e.g., rho=0.8) but different sigma matrices (different hub network structures from `generate_gram_hub_matrix()` with different seeds), combined into one matrix — a mixed population of internally regulated but transcriptomically distinct cells
   - Compute GMP-Cor on the full mixed population in Scenario B and show it is *higher* than Scenario A, not lower. This directly addresses the reviewer's concern: a mixture of regulated subpopulations produces an elevated GMP-Cor, while genuine dysregulation produces a low one.

2. **Manuscript text:** Add a paragraph in the Results/Methods clarifying that scRNA-seq is essential because GMP-Cor requires per-cell gene expression vectors — bulk data gives one average expression vector per condition and cannot compute the correlation spectrum across cells.

---

## Comment 2.1 — Simulations for Expected GMP-Cor Ranges

**What the reviewer is asking:** What GMP-Cor values should be expected for fully regulated, partially regulated, and dysregulated states? The broad spread in Fig. 2D-F needs to be contextualized.

**Proposed response steps:**

1. **Sweep rho using `simulate_scRNA_data()`:** Run the synthetic data generator across a range of correlation strengths (e.g., rho = 0.1, 0.3, 0.5, 0.7, 0.9). For each rho, generate a cell × gene count matrix and compute GMP-Cor using the new definition (sum of eigenvalue excesses above max scrambled eigenvalue).

2. **Plot GMP-Cor vs. rho:** A single figure showing how GMP-Cor monotonically increases with rho. This defines a calibration curve: what correlation strength corresponds to each biological condition?

3. **Overlay experimental values:** Superimpose the experimentally observed GMP-Cor values (Exp, Reg-Arrest, Dis-Arrest) on this calibration curve to show they fall in interpretable positions.

4. **Noise vs. coupling disambiguation:** Run two additional series:
   - Fixed rho, varying dropout rate — show GMP-Cor changes little with dropout noise
   - Fixed dropout, varying rho — show GMP-Cor tracks correlation strength
   - This directly addresses whether the metric reflects dysregulation vs. enhanced technical variability.

---

## Comment 2.2 — Positive Control / Benchmark

**What the reviewer is asking:** Is there a positive control with a known phenotype that can anchor the GMP-Cor scale?

**Proposed response steps:**

1. **Use simulated extremes as anchors:** Generate synthetic data with rho → 1 (nearly fully correlated, "maximally regulated") and rho → 0 (uncorrelated, "maximally dysregulated") as theoretical upper and lower bounds for GMP-Cor. These simulated extremes bracket the calibration curve and give a concrete scale.

2. **Link to experimental conditions:** Show that exponential growth cells have GMP-Cor consistent with the high-rho regime, while Dis-Arrest cells fall near the low-rho end.

3. **Quantitative statement:** Provide a table of GMP-Cor values per condition (mean ± SD across replicates) alongside the equivalent rho from the calibration curve. This gives readers a concrete interpretation of the metric magnitude. Publicly available E. coli scRNA-seq datasets from GEO (exponential growth) may serve as an independent supplementary validation.

---

## Comment 2.3 — GMP Computed on Full Population; Cluster Sensitivity

**What the reviewer is asking:** Fig. 1H shows Dis-Arrest cells split into two clusters. Was GMP-Cor computed on all cells together? Is the metric sensitive to which cells are included?

**Proposed response steps:**

1. **Re-run GMP-Cor per cluster:** For the two Dis-Arrest clusters identified in Fig. 1H, compute GMP-Cor separately. If both clusters show similarly low GMP-Cor, this demonstrates the metric is not an artifact of mixing. Report this result explicitly.

2. **Subsampling robustness (computational):** Using `simulate_scRNA_data()` with a fixed rho, subsample cells at varying fractions (e.g., 10%, 25%, 50%, 75%, 100% of n_cells) and compute GMP-Cor at each level. Because the ground truth correlation structure is exactly known from the simulation, this is a clean test. Show that above some minimum threshold (likely ~50–100 cells) GMP-Cor stabilizes.

3. **Simulated mixed population:** Generate two synthetic sub-populations (high rho, different sigma matrices) and vary their mixing ratio (10/90, 25/75, 50/50, etc.). Show that GMP-Cor remains elevated regardless of mixing ratio. Compare to a single low-rho population at the same total cell count to demonstrate the subpopulation/dysregulation contrast.

4. **Marker gene analysis for clusters:** Run Scanpy marker gene analysis on the two Dis-Arrest clusters (already in `scanpy_analysis.ipynb`) and show they are not growing/non-growing but represent biologically similar states. Cite the GMP-Cor per cluster as evidence.

---

## Comment 2.4 — Generalizability Across Conditions and Data Types

**What the reviewer is asking:** The metric is only tested on two conditions. Can it generalize to other datasets, technologies, and bulk RNA-seq?

**Proposed response steps:**

1. **Apply GMP-Cor to publicly available E. coli scRNA-seq datasets:** Several E. coli scRNA-seq datasets exist in GEO. Pick 1–2 from well-defined growth states and compute GMP-Cor. Show that exponential growth → high GMP-Cor, stationary/stress → lower GMP-Cor, even in independent datasets with different protocols.

2. **Address the bulk RNA-seq question conceptually:** Explain in the manuscript that GMP-Cor is inherently a single-cell metric — it computes correlations *across cells*, which is not possible with bulk data (one sample = one point). Bulk data provides a different view via enrichment scores or differential expression.

3. **Cell count scaling analysis:** Show empirically (from the subsampling analysis in step 2.3 above) that the method works for both large and small cell counts, defining the practical lower bound.

---

## Comment 2.5 — Comparison to Entropy-Based Metrics

**What the reviewer is asking:** Prior work (Zhu et al. 2020, doi: 10.1038/s41467-020-18134-z) used transcriptional entropy on bulk RNA-seq to quantify dysregulation. GMP-Cor is conceptually related. Show advantages/limitations vs. entropy. Also note that the prior work was misrepresented in the Introduction.

**Proposed response steps:**

1. **Fix the Introduction:** Correct the description of the entropy paper's conclusions. Acknowledge that the entropy approach captures global dysregulation, and explicitly position GMP-Cor as complementary in the single-cell context.

2. **Explain why entropy is not applicable here (manuscript text):** Transcriptional entropy quantifies how broadly cells sample gene expression state space over time or across perturbations — it is a measure of dynamic disorder. In antibiotic-arrested conditions, cells are in a steady state: they are not traversing state space. Entropy of a steady-state snapshot is therefore not a meaningful readout of dysregulation. GMP-Cor, by contrast, directly quantifies whether the steady-state correlation structure is preserved. The two metrics address fundamentally different questions and operate on different data types (time-series/bulk trajectories vs. single-cell snapshots). No simulation comparison is needed — this is a conceptual distinction that should be stated clearly in the manuscript.

3. **Address the Jensen et al. 2017 paper** (flagged also by Reviewer #4): Add a brief discussion paragraph comparing the approach to Tn-seq entropy comparisons.

---

## Comment 3 — GO Term Analysis on scRNA-seq Data

**What the reviewer is asking:** GO enrichment was only done on bulk RNA-seq. Why not on scRNA-seq? And how correlated are the two datasets?

**Proposed response steps:**

1. **Run DE analysis on scRNA-seq clusters:** Use the existing Scanpy framework (`scanpy_analysis.ipynb`) to perform differential expression between Reg-Arrest and Dis-Arrest clusters. The marker gene output already exists (referenced in `bulk_functions.run_go_single_cell()`).

2. **Run GO enrichment on scRNA-seq DE genes:** `bulk_functions.run_go_single_cell()` is already implemented and uses the same goatools pipeline. Apply it to the scRNA-seq DE results. Compare GO terms enriched in scRNA-seq vs. bulk and report the overlap.

3. **Correlate scRNA-seq and bulk expression:** Compute Pearson/Spearman correlation between per-gene average expression from scRNA-seq and the bulk RNA-seq log-fold changes. A scatter plot showing strong correlation will support internal consistency.

---

## Comment 4 — Specificity of Fig. 4A; GO Term Parameter Sensitivity

**What the reviewer is asking:** Are the upregulated genes in Fig. 4A truly specific to chemotaxis/flg/fli? How sensitive are enrichment scores to GO term set-size parameters? Fix Fig. 4B x-axis labels.

**Proposed response steps:**

1. **GO term parameter sensitivity analysis:** Re-run the GSEA/GO enrichment from `bulk_analysis.ipynb` with varying minimum and maximum gene-set sizes (e.g., min 5, 10, 15; max 100, 200, 500). Show that the direction and significance of differences between Dis-Arrest and Reg-Arrest are stable across parameter choices.

2. **Annotate Fig. 4A genes to all GO terms:** For the genes highlighted in Fig. 4A, list all their GO term memberships (not just chemotaxis). Show whether the pattern is specific to motility/chemotaxis or reflects broader dysregulation.

3. **Fix Fig. 4B x-axis:** Replace GO term numeric IDs with "GO term name (n=X genes)" format. This is a straightforward matplotlib label update in the figure script.

---

## Comment 5 — VapC Discrepancy: Long Lag Times vs. Dis-Arrest-Like GMP-Cor

**What the reviewer is asking:** All VapC time points show prolonged lag, but only 24h shows Dis-Arrest-like GMP-Cor. Does the approach actually predict persistence? Direct persistence assay needed. Are lag times in Fig. 5H comparable to Figs. 1E and 4D?

**Proposed response steps:**

1. **Reframe the narrative:** Clarify in the text that the GMP-Cor metric does not predict *whether* cells have long lags, but rather *why* — distinguishing genuine transcriptional dysregulation (Dis-Arrest-like) from a regulated prolonged arrest. The 2h/5h VapC points have long lags but retain coordination; the 24h point has lost coordination. This is a finding, not a contradiction.

2. **Compare absolute lag times across figures:** Add a supplementary table with mean ± SD lag times in minutes for all conditions across Figs. 1E, 4D, and 5H to confirm whether the 24h VapC lag times are quantitatively similar to Dis-Arrest lags.

3. **Address the persistence assay request:** Either (a) perform a biphasic killing curve on VapC-induced cells at different time points and show that the 24h sample has higher persister frequency, or (b) clearly acknowledge in the Discussion that a direct persistence assay was not performed and that this is a limitation, while explaining why the correlation between long lag and persistence is supported by prior literature.

4. **Fix the Fig. 5H legend:** Add the missing fourth label to the legend. This is a trivial fix in the figure script.

---

## Minor Comments

- **Minor 1 (Fig. 5H legend):** Edit the figure script to add the fourth distribution label. Trivial fix.
- **Minor 2 (package versions):** Run `conda env export` or `pip freeze` and add a `session_info` call at the end of key notebooks to embed exact package versions. `session-info` is already in the environment.

---

## Simulation Improvement Ideas

All simulations should use `simulate_scRNA_data()` from `simulations/simulated_data.ipynb`. The dynamical simulations (`main.py`, `competitive_binding.py`, `subpopulations.py`) are not reliable for this purpose — the correlation signal is dominated by high data dispersion rather than true gene-gene couplings.

### 1. Add a rho sweep runner
Add a `rho_sweep()` function to `simulated_data.ipynb` (or a new `simulations/sweep.py`) that:
- Loops rho over e.g. [0.1, 0.3, 0.5, 0.7, 0.9]
- At each rho, calls `simulate_scRNA_data()` and computes GMP-Cor = `sum(pcs[pcs > pcs1.max()] - pcs1.max())`
- Returns a DataFrame of rho → GMP-Cor
- Produces the calibration curve for Comments 2.1 and 2.2

### 2. Add a subpopulation mixing scenario
Extend `simulated_data.ipynb` to simulate two sub-populations with different sigma matrices (different seeds for `generate_gram_hub_matrix()`), each with high rho, combined into one matrix. Vary the mixing ratio (10/90 → 50/50). Compare GMP-Cor of the mixture to a single-population low-rho case. Produces the core figure for Comments 1 and 2.3.

### 3. Add a cell-count subsampling loop
Fixed rho, vary n_cells from 50 → 1000. Compute GMP-Cor at each level. Plot to define the minimum cell count for reliable metric estimation. Directly addresses Comment 2.3 robustness question.