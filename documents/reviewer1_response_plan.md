# Response Plan: Reviewer #1 Comments

## Overview

Reviewer #1 has two broad concerns: (1) the biological interpretation and experimental design, and (2) the GMP-Cor metric validation. The metric-validation concerns (Comment 2) are the most tractable computationally and map directly onto the existing simulations. Comments 1, 3, 4, and 5 require either re-analysis of existing data or targeted experimental work.

---

## Comment 1 — Why scRNA-seq? Subpopulations vs. Uniform Dysregulation

**What the reviewer is asking:** Could the observed loss of gene-gene correlations arise from a mixture of two regulated subpopulations rather than uniform dysregulation? And why is scRNA-seq the right tool here?

**Proposed response steps:**

1. **Simulation argument (computational):** Use the existing `subpopulations.py` framework to generate two contrasting scenarios:
   - **Scenario A (dysregulation):** Single population, low coupling strength J → genuinely low GMP-Cor
   - **Scenario B (subpopulations):** Two subpopulations each with high J but different interaction networks (as in `subpopulations.py`) → mixed population
   - Compute GMP-Cor on the full mixed population in Scenario B and show it is *higher* than Scenario A, not lower. This directly addresses the reviewer's concern: subpopulations increase apparent correlation under GMP, while genuine dysregulation decreases it.

2. **Comparison to bulk entropy:** For the same two scenarios, compute bulk-level transcriptional entropy (from `utils.calculate_entropy()`) on "pseudo-bulk" samples aggregated from the simulated cells. Show that bulk entropy cannot distinguish Scenario A from Scenario B, but GMP-Cor can. This also addresses Comment 2.5 about entropy comparisons.

3. **Manuscript text:** Add a paragraph in the Results/Methods clarifying that scRNA-seq is essential because GMP-Cor requires per-cell gene expression vectors — bulk data gives one average expression vector per condition and cannot compute the correlation spectrum across cells.

---

## Comment 2.1 — Simulations for Expected GMP-Cor Ranges

**What the reviewer is asking:** What GMP-Cor values should be expected for fully regulated, partially regulated, and dysregulated states? The broad spread in Fig. 2D-F needs to be contextualized.

**Proposed response steps:**

1. **Sweep J in `main.py`:** Run the main simulation across a range of coupling strengths (e.g., J = 0, 0.5, 1, 2, 3, 5). For each J, take the last N_STEPS snapshot (the steady-state distribution across trajectories) and run `utils.get_eig_dist()` on it as if it were scRNA-seq data. Compute GMP-Cor from the eigenvalue distribution.

2. **Plot GMP-Cor vs. J:** A single figure showing how GMP-Cor monotonically maps to coupling strength. This defines the "ruler" for interpreting experimental values: what J corresponds to each biological condition?

3. **Overlay experimental values:** Superimpose the experimentally observed GMP-Cor values (Exp, Reg-Arrest, Dis-Arrest) on this calibration curve to show they fall in sensible positions.

4. **Noise vs. coupling disambiguation:** Run two additional series:
   - Fixed J, varying noise amplitude (from `main.py`: `noise_amp` parameter) — show GMP-Cor changes little
   - Fixed noise, varying J — show GMP-Cor tracks coupling
   - This directly addresses whether the metric reflects dysregulation vs. enhanced variability.

---

## Comment 2.2 — Positive Control / Benchmark

**What the reviewer is asking:** Is there a positive control with a known phenotype that can anchor the GMP-Cor scale?

**Proposed response steps:**

1. **Use exponential growth datasets as positive control:** Apply GMP-Cor to publicly available E. coli scRNA-seq data from exponentially growing cells (already mentioned in the draft rebuttal). Present as the "maximally correlated" reference.

2. **Link to simulation calibration:** After establishing the J-vs-GMP-Cor curve (step above), show that exponential growth cells have a GMP-Cor consistent with a high-J (well-coupled) regime. Dis-Arrest cells fall at the low-J end.

3. **Quantitative statement:** Provide a table of GMP-Cor values per condition (mean ± SD across replicates) alongside the simulated J-equivalent. This gives readers a concrete interpretation of the metric magnitude.

---

## Comment 2.3 — GMP Computed on Full Population; Cluster Sensitivity

**What the reviewer is asking:** Fig. 1H shows Dis-Arrest cells split into two clusters. Was GMP-Cor computed on all cells together? Is the metric sensitive to which cells are included?

**Proposed response steps:**

1. **Re-run GMP-Cor per cluster:** For the two Dis-Arrest clusters identified in Fig. 1H, compute GMP-Cor separately. If both clusters show similarly low GMP-Cor, this demonstrates the metric is not an artifact of mixing. Report this result explicitly.

2. **Subsampling robustness (computational):** Using any of the existing simulations (e.g., `main.py`), subsample cells at varying fractions (10%, 25%, 50%, 75%, 100% of trajectories) and compute GMP-Cor at each level. Show that above some minimum threshold (likely ~50–100 cells), GMP-Cor stabilizes. This directly answers the reviewer's question about robustness to population composition.

3. **Simulated mixed population:** Using `subpopulations.py`, vary the mixing ratio of two sub-populations (10/90, 25/75, 50/50, etc.) and show that GMP-Cor remains elevated regardless of mixing ratio (since both sub-populations are individually regulated). Compare to a single dysregulated population at the same total cell count.

4. **Marker gene analysis for clusters:** Run Scanpy marker gene analysis on the two Dis-Arrest clusters (already in `scanpy_analysis.ipynb`) and show they are not growing/non-growing but represent biologically similar states. Cite the GMP-Cor per cluster as evidence.

---

## Comment 2.4 — Generalizability Across Conditions and Data Types

**What the reviewer is asking:** The metric is only tested on two conditions. Can it generalize to other datasets, technologies, and bulk RNA-seq?

**Proposed response steps:**

1. **Apply GMP-Cor to publicly available E. coli scRNA-seq datasets:** Several E. coli scRNA-seq datasets exist in GEO. Pick 1–2 from well-defined growth states and compute GMP-Cor. Show that exponential growth → high GMP-Cor, stationary/stress → lower GMP-Cor, even in independent datasets with different protocols.

2. **Address the bulk RNA-seq question conceptually:** Explain in the manuscript that GMP-Cor is inherently a single-cell metric — it computes correlations *across cells*, which is not possible with bulk data (one sample = one point). What bulk data provides is a different (and as shown, less informative) view via entropy or enrichment scores.

3. **Cell count scaling analysis:** Show empirically (from the subsampling analysis in step 2.3 above) that the method works for both large and small cell counts, defining the practical lower bound.

---

## Comment 2.5 — Comparison to Entropy-Based Metrics

**What the reviewer is asking:** Prior work (Zhu et al. 2020, doi: 10.1038/s41467-020-18134-z) used transcriptional entropy on bulk RNA-seq to quantify dysregulation. GMP-Cor is conceptually related. Show advantages/limitations vs. entropy. Also note that the prior work was misrepresented in the Introduction.

**Proposed response steps:**

1. **Fix the Introduction:** Correct the description of the entropy paper's conclusions. Acknowledge that the entropy approach also captures global dysregulation, and explicitly position GMP-Cor as complementary or superior in the single-cell context.

2. **Head-to-head simulation comparison:** This is the key new analysis. Using the simulation framework:
   - Generate 3 datasets from simulations: (A) regulated, (B) dysregulated, (C) mixed regulated subpopulations
   - Compute both GMP-Cor and bulk transcriptional entropy on each
   - Show the key result: entropy conflates (B) and (C) as both "disordered," while GMP-Cor correctly distinguishes them
   - `utils.calculate_entropy()` already exists for this purpose (though it currently imports sklearn which may need enabling)

3. **Apply both metrics to real data:** Run both GMP-Cor and entropy on the existing experimental scRNA-seq datasets (Exp, Reg-Arrest, Dis-Arrest) and compare their values side by side. Entropy from bulk will need the bulk RNA-seq data from `bulk_analysis.ipynb`.

4. **Address the Jensen et al. 2017 paper** (flagged also by Reviewer #4): Add a brief discussion paragraph comparing the approach to Tn-seq entropy comparisons.

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

The current simulations (`simulations/`) produce raw time-series arrays but do not auto-compute GMP-Cor or produce comparison figures. Here are targeted improvements aligned with the rebuttal:

### 1. Add a GMP-Cor computation pipeline to each simulation
Currently, simulations save `.npy` files but don't analyze them. Add a post-processing step that:
- Takes the last `T` time steps of each trajectory as a "snapshot" (analogous to scRNA-seq)
- Passes the cell×gene matrix to `utils.get_eig_dist()`
- Computes a scalar GMP-Cor value (e.g., ratio of max eigenvalue to MP prediction, or the entropy measure)
- Returns that value for downstream comparison

### 2. Add a parameter sweep runner
Create a new script (e.g., `simulations/sweep.py`) that loops over J values and records GMP-Cor at each J. This directly produces the calibration curve for Comment 2.1. The `main.py` and `competitive_binding.py` simulations are the best candidates as they have a clean J parameter controlling coupling strength.

### 3. Unify the subpopulations vs. dysregulation comparison
`subpopulations.py` currently simulates subpopulations but doesn't compare GMP-Cor to a dysregulated single-population control. A new script (or notebook) should:
- Run `subpopulations.py` with RATIO=1.0 (fully different networks → regulated heterogeneity)
- Run `main.py` with low J (same number of cells, but uniformly dysregulated)
- Run both through the GMP-Cor pipeline
- Plot side-by-side eigenvalue spectra and GMP-Cor values
This produces the core figure addressing Comment 1 and Comment 2.3.

### 4. Enable the entropy comparison in `utils.py`
`utils.calculate_entropy()` imports sklearn but that import is commented out at the top. Once enabled, use it on pseudo-bulk aggregates from both the regulated and dysregulated simulation outputs to produce the entropy vs. GMP-Cor comparison figure (Comment 2.5).

### 5. Fix the global variable pattern in simulation scripts
`main.py`, `subpopulations.py`, etc. reference `degradation_rate`, `noise_amp`, and `J` as module-level globals inside `@njit` functions. This works but makes parameter sweeps awkward. Consider passing them as arguments or defining them as constants before the JIT-compiled function, so that a sweep over J doesn't require re-compiling the Numba kernel.

### 6. Add a subsampling robustness script
Add a loop in any simulation script that varies `n_trajectories` from small (50) to large (1000), computes GMP-Cor at each size, and plots the result. This directly addresses Comment 2.3 robustness question.