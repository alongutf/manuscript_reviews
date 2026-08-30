# Reviewer #1 — subpopulations: state of the analysis and a proposed decisive test

**Date:** 2026-08-30
**Scope:** R1 major comment 1, comment 2.3 (both parts). Related: `reviewer1_response_plan.md`,
`subsampling_analysis.md`, `gmp_cor_provenance_analysis.md`.

---

## 1. What the reviewer actually asks

Three distinct questions, currently answered with one argument:

| | Question | Where |
|---|---|---|
| **Q1** | Could the loss of correlation we call dysregulation instead arise from *distinct cellular subpopulations* rather than uniform dysregulation? | Comment 1 |
| **Q2** | Fig. 1H splits Dis-Arrest into two clusters, yet GMP-Cor is computed on the whole population — is the index sensitive to *which* cells are included? | Comment 2.3a |
| **Q3** | Why not compute correlations *per subpopulation* — that would be the strongest argument for using scRNA-seq. Have you varied the *proportions* of cells? | Comment 2.3b |

Q3 is the one currently answered least directly: the rebuttal answers "we do not observe
subpopulations", then pivots to a cell-number subsampling analysis on synthetic data. The
reviewer asked about *composition*, not *count*.

---

## 2. What already exists

**Simulation.** `src/simulations.py:subpopulation_mixing()`, driven by
`simulations/subpopulation_mixing_run.py` (rho_high = 0.8) and
`simulations/subpopulation_mixing_rho09_run.py` (rho_high = 0.9, seeds chosen to match the two
sub-populations' individual GMP-Cor). Both are **50/50 only, one realization, no error bars**.

| condition (rho = 0.9 run) | sub-pop A | sub-pop B | 50/50 mixture |
|---|---|---|---|
| regulated (rho = 0.9) | 78.03 | 81.74 | **30.71** |
| dysregulated (rho = 0.1) | 0.00 | 0.02 | **4.26** |

**Experimental.** `simulations/dataset_mixing_gmpcor_run.py` — Exp (`EXP_biorep_t0A`) + VapC-2h
(`VAPC_biorep_t2A`), 500 cells each, shared gene space, top-2000 Fano genes -> 1000 x 2000,
**GMP-Cor = 58.88** (lambda_max 32.4, lambda*_scr 6.16, 13 signal modes). This is the "58.88"
quoted in the rebuttal.

**Per-cluster.** `results/cluster_gmp_cor/` — GMP-Cor per sample x leiden cluster, from both
`data_for_paper` and `data/scanpy_shx.h5ad`.

**Robustness.** `documents/subsampling_analysis.md` — four cell/gene subsampling experiments at
rho = 0.9.

---

## 3. Where this is vulnerable

These are the things a determined reviewer (or a reader of the deposited logs) will find.

1. **The headline claim was not what the existing simulation showed.** The rebuttal says
   mixing subpopulations *elevates* GMP-Cor. In the rho = 0.9 run the mixture (30.71) is
   **2.6x lower** than either pure sub-population (78, 82).
   *Resolved (2026-08-30):* two separate causes, both now fixed.
   (a) The pure sub-populations there were evaluated at 500 cells and the mixture at 1000,
   and the scrambled threshold is a pure function of matrix shape — the comparison was
   invalid. (b) The two sub-populations were not actually distinct: they had different
   networks but independently-drawn expression profiles from the same prior, so they barely
   separated in PCA and no between-population mode formed. Mixing then only diluted one
   network with the other. Both are addressed by
   `simulations/inverted_subpopulation_mixing_run.py` (§4.1a below).
2. **A mixture is not distinguishable from a single population by the scalar alone.** The seed
   screen in `subpopulation_mixing_rho09_50_50_20260615_111850.json` gives GMP-Cor from **7.3 to
   176.6 across sigma seeds at fixed rho = 0.9** — a 24x spread driven purely by network topology.
   The mixture value (30.7) sits comfortably inside the range of *single* regulated populations.
   So "GMP-Cor = 30" cannot by itself be read as "not a mixture", and by extension a low
   Dis-Arrest value cannot be defended as "not a mixture" on magnitude alone. This also
   qualifies the calibration curve: mapping an observed GMP-Cor onto a chi value assumes a
   topology we have not constrained.
3. **58.88 had no matched reference.** The number was reported without the two pure
   populations computed in the *same* frame. Since GMP-Cor is extensive in *p* and scales with
   detection depth (`gmp_cor_provenance_analysis.md` §4), an unreferenced absolute value carries
   no argument.
   *Resolved (2026-08-30):* `simulations/dataset_mixing_ratio_run.py`, run
   `dataset_mixing_ratio_20260830_122400`. All points at n = 900 cells on the same fixed
   2000-gene panel (`data_for_umap`, top-Fano on the combined pool, reporter genes
   excluded), 5 independent cell draws each:

   | Exp / VapC-2h | GMP-Cor | centered | dGMP | group-axis AUC |
   |---|---|---|---|---|
   | pure VapC-2h | 29.22 ± 4.20 | — | — | — |
   | 25 / 75 | 41.61 ± 4.17 | 25.48 ± 5.35 | 0.39 | 1.00 |
   | **50 / 50** | **40.00 ± 4.41** | 19.04 ± 3.60 | **0.53** | 1.00 |
   | 75 / 25 | 40.78 ± 1.73 | 23.29 ± 4.19 | 0.43 | 1.00 |
   | pure Exp | 18.66 ± 1.38 | — | — | — |

   **Every mixture exceeds both pure populations** (1.37-1.42x the larger), and
   group-mean centering brings each back into the endpoint range. The rebuttal's premise
   is therefore correct once the endpoints are computed at matched n; the earlier 58.88
   was a 1000-cell mixture compared against nothing. (Note the ordering: at matched n,
   pure VapC-2h is *higher* than pure Exp.)

   Superseded: an earlier single-draw estimate at n ≈ 1000 gave 16.52 / 28.46 / 41.22
   with dGMP 0.55. It agrees with the repeated run to within the scatter, but came from
   an ad-hoc script that was not kept; quote the run-log values above.
4. **The experimental mixture is confounded by batch.** Exp and VapC-2h are different libraries
   with different depth and detection. The leading PCs of that mixture may simply encode
   library, not biology — which would inflate GMP-Cor for a reason that has nothing to do with
   subpopulations. There is no within-library mixture control.
5. **Only one mixing ratio, in both simulation and data.** Q3 explicitly asks about varying
   proportions. 10/90 is the interesting case: a small dysregulated minority inside a regulated
   majority (or vice versa) is exactly the "growing + non-growing" scenario in Fig. 1H.
6. **The per-cluster numbers are not yet usable.** `cluster_gmp_cor` and
   `cluster_gmp_cor_h5ad` disagree by up to 16x for the same cluster, and the clusters differ in
   *n*, *p* and detection — the confounds documented in `gmp_cor_provenance_analysis.md` §4, §8.
   Nothing per-cluster should be quoted until it is run under the `equate_reg_n1000` protocol.
7. **No positive discriminator.** Everything so far argues from magnitude. Nothing yet identifies
   mixture structure *as such*.

---

## 4. The proposed analysis

**The idea.** Mixture and dysregulation leave different signatures in the spectrum, and the
difference is structural, not scalar:

- A **mixture** of two internally-coordinated populations adds a small number of **low-rank,
  between-group** modes. Their eigenvectors align with the group-mean difference, and the cell
  scores on those modes are **bimodal**. Remove the group means and the modes vanish.
- **Dysregulation** removes **within-cell** coordination across the whole spectrum. It is
  distributed (consistent with `supplementary_note_eigenvector.tex`), the cell scores are
  unimodal, and it is unaffected by group-mean centering.

So the decisive observable is not GMP-Cor but **how GMP-Cor responds to group-mean removal**,
together with the modality of the leading cell scores. Define

> **dGMP = 1 − GMP-Cor(group-centered) / GMP-Cor(raw)**

computed with cells centered to their own cluster mean before correlating. A mixture-driven
signal has dGMP ≈ 1; genuine within-cell coordination has dGMP ≈ 0.

This turns Q3 into a strength rather than a defense: it is a *single-cell-only* test — bulk data
cannot compute it — which is exactly the argument the reviewer says would make the case for
scRNA-seq.

### 4.1a Built (2026-08-30): the inverted-subpopulation scenario

`simulations/inverted_subpopulation_mixing_run.py`, on new library code in
`src/simulations.py`. This replaces S1 for the regulated x regulated case and fixes
weakness 1.

**How the two populations are made distinct.** Sub-population B is sub-population A with
its **expression ranking inverted** (`invert_gene_means`): A's most lowly-expressed gene is
B's most highly-expressed, and so on. The *multiset* of gene means is preserved exactly, so
both populations have identical marginal expression distributions, dynamic range and
sparsity — only the assignment of expression level to gene is reversed. Any separation
between them is therefore a real difference in which genes are expressed, and cannot be
dismissed as a depth or library-size artefact. `simulate_scRNA_data()` gained a `gene_mu`
argument so an explicit profile can be supplied (default behaviour unchanged).

**Two network configurations.**
`shared` — both populations use the same hub network, differing only in expression profile.
This is the experimental case (two states of one organism share regulatory architecture).
`distinct` — each also gets its own topology; the two networks dilute one another, so this
is the conservative case.

**Calibration matters more than expected.** At the module's usual `inv_gamma_scale = 0.01`
the simulated cells detect only ~24 genes each and just **6 genes** carry any material
between-population difference — no separating mode can form and the scenario cannot be
tested at all. At `inv_gamma_scale = 0.04` cells detect ~85 genes and **93 genes** differ,
matching the experimental matrices (88–155 detected genes/cell, 91 differing genes in the
Exp/VapC-2h panel). The runner uses 0.04 and states why.

**Every ratio is evaluated at the same total n**, because the scrambled threshold is a pure
function of matrix shape.

**Measured at each ratio:** GMP-Cor; GMP-Cor after group-mean centering; dGMP; group-axis
AUC, bimodality and which spectral mode carries the separation; plus a single dysregulated
population of the same size as reference. Note the separation statistic is computed on the
**group axis**, not PC1 — when a population has strong internal structure the separating
direction is often PC2, and a PC1-only statistic reads as "no separation" for two perfectly
separated groups. (This is why the earlier runs looked like there was no separation.)

**Results** (`inverted_subpopulation_mixing_20260830_105044`, n = 1000 cells x 2000 genes,
rho_high = 0.7, 50/50 mixture, 5 repeats, mean ± SD):

| | shared network | distinct network |
|---|---|---|
| pure sub-pop B | 44.1 ± 3.0 | 42.3 ± 4.2 |
| **50/50 mixture** | **68.8 ± 3.3** | **56.2 ± 2.1** |
| pure sub-pop A | 43.7 ± 3.5 | 43.1 ± 3.2 |
| mixture, group-centered | 45.3 ± 3.1 | 32.3 ± 1.9 |
| dGMP | 0.34 | 0.43 |
| mixture / pure | 1.57x | 1.32x |
| group-axis AUC · bimodality | 1.000 · 0.89 | 1.000 · 0.90 |
| **dysregulated single population** | **3.99 ± 0.47** | |

`rho_high = 0.7` is chosen so a pure sub-population lands at GMP-Cor ~43, inside the
experimental range for regulated samples (~30-50). (At the calibrated expression scale:
rho 0.5 -> 16, 0.6 -> 25-28, 0.7 -> 42-43, 0.8 -> 66-72, 0.9 -> ~130.)

1. **The populations are genuinely distinct**: group-axis AUC = 1.000, bimodality
   coefficient 0.89-0.90 (threshold 5/9), two clearly separated clusters in UMAP.
2. **Mixing inflates GMP-Cor** in both configurations — the mixture exceeds *both* pure
   populations, by 1.57x (shared network) and 1.32x (distinct). Group-mean centering
   returns the shared-network mixture to 45.3, i.e. exactly the pure-population level.
3. **The mixture is nowhere near dysregulation**: the lowest single realization is 54.3,
   **13.6x** the dysregulated reference of 3.99.
4. **The simulation now matches the experiment quantitatively.** Simulated
   mixture/pure = 1.57x vs experimental 41.22/28.46 = 1.45x; simulated dGMP 0.34-0.43 vs
   experimental 0.55. Both the direction and the magnitude of the effect agree.

Note GMP-Cor carries ~2% run-to-run noise from the scramble draw alone (`get_eig_dist`
permutes with the global RNG, unseeded): re-running the identical seeded mixture gave
67.21 vs the recorded 68.57. The ±SD over repeats already covers this.

Files: `results/simulation_results/{logs,raw,figures}/inverted_subpopulation_mixing_20260830_105044.*`
(the `.csv` holds every per-repeat record; the `.json` holds the UMAP coordinates and gene-mean
profiles, so the figure can be re-plotted without re-simulating).

**The dysregulated scenario** (`inverted_subpopulation_mixing_dysregulated_20260830_111402`,
`--scenario dysregulated`: rho_high = 0.1 for *both* sub-populations, shared network, 5 repeats).
Same construction, but neither sub-population has any internal coordination to begin with:

| | GMP-Cor |
|---|---|
| pure sub-pop B | 3.26 ± 0.7 |
| pure sub-pop A | 4.08 ± 0.7 |
| **50/50 mixture** | **29.68 ± 0.6** (8.1x the pure populations) |
| mixture, group-centered | 4.46 ± 0.9 → **dGMP = 0.85** |
| regulated single population (reference) | 43.20 ± 4.67 |

Separation is again complete (group-axis AUC 1.000, bimodality 0.91). Two conclusions:

1. **Mixing never lowers GMP-Cor.** Even two populations with *no* internal coordination,
   sitting at the noise floor, give a mixture 8x higher than either alone — and centering
   returns it exactly to the floor (4.46 vs 3.3/4.1), so essentially all of it is the
   separating mode. This is the reviewer's concern answered from the opposite direction:
   heterogeneity is not a mechanism that can produce a *low* GMP-Cor.
2. **The bound on the claim.** A mixture of two dysregulated populations (29.7) lands in the
   same range as a single regulated population (43.2). The raw index alone therefore cannot
   distinguish them — but dGMP can: 0.85 for the mixture versus 0.01–0.05 for a pure
   population. This is the strongest argument for reporting dGMP alongside GMP-Cor whenever a
   sample shows cluster structure, and it is a single-cell-only measurement.

**Note on the metric's floor.** At this data scale a pure population with rho = 0.1 gives
GMP-Cor ~4, and so does rho = **0.0** (an identity correlation matrix, genes independent by
construction). The ~4 is therefore not residual coupling but a floor of the metric itself —
presumably compositional structure from row normalisation and the lognormal library-size
term, which the scrambled threshold does not fully absorb. Worth stating wherever a
near-zero GMP-Cor is interpreted.

**Still open — the ratio sweep.** An earlier run at rho = 0.9 over seven ratios found the
effect is **non-monotonic**: 10/90 and 90/10 gave the *highest* GMP-Cor and 50/50 the
lowest, so a minority sub-population is the worst case for the reviewer's hypothesis, not
the best — directly relevant to Q3, which asks about varying proportions. That run has been
superseded by the rho = 0.7 configuration and the finding is not currently backed by an
output file. Re-running the sweep at rho = 0.7 would restore it.

### 4.1 Simulations (ground truth known)

All on `simulate_scRNA_data()` / `generate_gram_hub_matrix()`, n = 1000 x 2000, dropout 1.0,
**10 realizations per point, mean ± SD** (fixing weakness 1 and the single-realization issue).

| | Run | Purpose |
|---|---|---|
| **S1** | *(done — see §4.1a)* **Mixing-ratio sweep.** Two regulated sub-pops (rho = 0.9, distinct sigma seeds) mixed at 0/100, 10/90, 25/75, 50/50, 75/25, 90/10, 100/0. Repeat over >= 10 seed pairs. | Gives the **mixture floor**: the *lowest* GMP-Cor two regulated populations can produce at any proportion. The argument becomes "Dis-Arrest lies below the floor", which is falsifiable and survives weakness 2. |
| **S2** | **Divergence sweep.** Interpolate Sigma_A -> Sigma_B (mixing fraction on the covariance, not just independent seeds) at fixed rho; sweep divergence 0 -> 1 at 50/50. | Answers "how different must two subpopulations be before mixing matters". Currently we test only the maximally-divergent case. |
| **S3** | **Asymmetric mixtures.** Regulated (rho = 0.9) + dysregulated (rho = 0.1) at each ratio. | The actual Fig. 1H scenario. Yields the mapping: what fraction of dysregulated cells is needed to pull a population to the observed Dis-Arrest value? A strong, quotable number. |
| **S4** | **Validate dGMP on ground truth.** Apply group-mean centering to (a) S1 mixtures, (b) a single dysregulated population, (c) a single regulated population. Also record dip-test p on PC1 cell scores. | Establishes sensitivity/specificity of the diagnostic and sets the threshold used on real data. Without S4 the dGMP test is unvalidated. |
| **S5** | **Depth/batch confound control.** Two sub-populations with *identical* Sigma and rho but different sequencing depth and dropout. | Quantifies how much GMP-Cor inflation is pure batch — the calibration needed to read E1/E2. Directly addresses weakness 4. |

### 4.2 Experimental analyses

All at **matched n, p and detection**, using the protocol in `gmp_cor_provenance_analysis.md` §7
(`scripts/equate_reg_n1000.py`), rRNA/tRNA removed, 5 draws per point.

| | Run | Purpose |
|---|---|---|
| **E1** | **Mixing-ratio series on real cells.** Extend `dataset_mixing_gmpcor_run.py` with a `--ratio` sweep over the same 7 ratios, on the shared gene panel, total n held fixed at 1000. Pairs: Exp x VapC-2h (the current pair), Exp x Dis-Arrest, Reg-Arrest x Dis-Arrest. Report all pure endpoints in the same frame. | Fixes weaknesses 3 and 5 in one run: gives 58.88 its missing reference points and answers "varying the proportions" literally. |
| **E2** | **Within-library mixture null.** Same sweep but mixing two leiden clusters from the *same* sample (e.g. reg1 cluster 0 x cluster 2), and mixing two biological replicates of the same condition (13b x 15b). | The batch control (weakness 4). Shows how much of E1's elevation is library rather than biology. Calibrated by S5. |
| **E3** | **dGMP + PC1 modality per condition.** For Exp, Reg-Arrest (x2 reps), Dis-Arrest (x2 reps), VapC series: GMP-Cor raw vs. cluster-centered, plus Hartigan dip test on PC1/PC2 cell scores, plus the fraction of GMP-Cor mass carried by modes whose eigenvector correlates with the cluster-mean difference. Include the artificial E1 50/50 mixture as a positive control. | **The decisive panel.** Prediction: artificial mixture dGMP ≈ 1, dip p < 0.01; Dis-Arrest dGMP ≈ 0, dip n.s. That is a positive demonstration that the Dis-Arrest signal is *not* between-group structure — which is what Q1 actually asks. |
| **E4** | **Per-cluster GMP-Cor, done under control.** Re-run `cluster_gmp_cor` for both Fig. 1H Dis-Arrest clusters at matched n/p/detection, bootstrap CI over 5+ draws, both sources reconciled. | Answers Q2 with numbers instead of an assertion. Supersedes the current `results/cluster_gmp_cor/` tables (weakness 6). Pair with the existing marker-gene evidence that the clusters are not growing vs. non-growing. |
| **E5** | **Leverage / eigenvector overlap.** Reuse the `results/eigenvector_analysis/` pipeline: do the top modes of Dis-Arrest load on the same distributed gene set as Reg-Arrest, or on cluster-discriminating genes? | Ties the subpopulation answer to the gene-level note already written for R4 C1, at no extra machinery cost. |

### 4.3 What the figure looks like

One new supplementary figure, five panels:

- **A** — S1/S3: GMP-Cor vs. mixing ratio, regulated x regulated and regulated x dysregulated,
  mean ± SD, with the experimental Dis-Arrest value and the mixture floor as horizontal bands.
- **B** — E1: the same sweep on real cells, all endpoints in a matched frame.
- **C** — E2: within-library mixture null overlaid on B.
- **D** — E3: dGMP per condition, bar chart, artificial mixture as positive control.
- **E** — E3: PC1 cell-score histograms — artificial mixture (bimodal) vs. Dis-Arrest (unimodal).

### 4.4 Text changes that follow

1. **Rewrite the Comment 1 / 2.3 response.** Replace "mixing subpopulations elevates GMP-Cor"
   with the defensible form: *mixing two regulated populations cannot bring GMP-Cor below
   [mixture floor from S1] at any proportion, whereas Dis-Arrest sits at [value]; and the
   group-centering test shows the Dis-Arrest signal is not between-group structure.*
2. **State the topology caveat** where the calibration curve is used: at fixed chi, GMP-Cor varies
   with network topology, so a chi read off the curve is an order-of-magnitude statement.
3. **Answer Q3 affirmatively.** Present dGMP as a per-subpopulation analysis that only
   single-cell data can support, rather than declining on the grounds that no subpopulations
   are observed.
4. Always report *p* (and preferably GMP-Cor/*p*) alongside every GMP-Cor value in this section.

---

## 5. Suggested order of work

1. **S1 + S4** — cheap, self-contained, and they decide whether dGMP is a usable diagnostic
   before any experimental work is spent on it.
2. **E1** — one runner extension; immediately repairs the weakest number in the current rebuttal.
3. **S3, S5, E2** — the controls that make E1 interpretable.
4. **E3** — the decisive panel, once S4 has set the threshold.
5. **E4, E5, S2** — completeness.

## 6. New files this implies

```
simulations/mixing_ratio_sweep_run.py          # S1, S2, S3
simulations/group_centering_validation_run.py  # S4
simulations/depth_confound_mixing_run.py       # S5
simulations/dataset_mixing_gmpcor_run.py       # extend: --ratio sweep, --report-endpoints (E1, E2)
scripts/group_centering_test.py                # E3 (dGMP, dip test, mode overlap)
scripts/cluster_gmp_cor_matched.py             # E4 (equate_reg_n1000 protocol, per cluster)
scripts/figures/figure_sX_subpopulations.py    # the 5-panel figure
```

Outputs follow the project convention: figures -> `results/simulation_results/figures/`,
summaries -> `.../raw/`, parameter logs -> `.../logs/`; experimental outputs ->
`results/cluster_gmp_cor/`.
