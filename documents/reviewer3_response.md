# Response to Reviewer #3 — draft

*Draft, 2026-08-20. Point-by-point suggestions only; not merged into
`documents/Rebuttal02032026.docx`.*

**New data / analyses this draft leans on**

| # | New material | Where it lives | Answers |
|---|---|---|---|
| N1 | Mother-machine microscopy of SHX-induced arrest (28 trenches, 197 tracks, 155 validated halt events, 62 divisions) | `microscopy/true_events.csv`, Fig. S7 | Major 3, 7 |
| N2 | RegulonDB sigma-38 (RpoS) sigmulon tested against every DESeq2 contrast (n = 344 genes, RegulonDB 14.5.0) | `documents/rpoS_regulon_deseq_analysis.md`, Fig. S11 | Major 2, 4, 5 |
| N3 | New ampicillin kill curves for VapC-arrested vs. Reg-Arrest cells, 3 biological replicates, MPN | `kill curves/20260719_VIGA24h_TolwoATC_Prep.xlsx`, Fig. 5B | Major 1, 5, 6 |
| N4 | Lag-time distributions for the symmetric control **Reg-Arrest + SHX** alongside Reg-Arrest and Dis-Arrest | `scanlag_data/CASP+SHX/`, Fig. S9 | Major 4, 7 |
| N5 | Single-cell morphology and reporter-expression quantification (cell area and constitutive fluorescence, 4 conditions, n = 378–943 cells) | `microscopy/all_positions_*.csv`, Fig. S6 | Major 5, 7 |

---

## General response

We thank the reviewer for a detailed and, in several places, decisive critique. Three of
the criticisms — that the abruptness of SHX arrest was asserted rather than measured, that
the Reg-Arrest control was not matched to the Dis-Arrest treatment, and that the RpoS
literature was absent — were correct, and each has been answered with new experiments or
new analysis rather than with argument. We have also accepted the reviewer's framing point:
the manuscript's demonstrated result is a **quantitative metric for global transcriptional
dysregulation** and its association with slow, heterogeneous recovery; the identification
of Dis-Arrest cells with bona fide persisters is a hypothesis the data support only in
part. The title, abstract and discussion have been rewritten accordingly (see Major 9).

---

## Major comments

### 1. Long lag times do not establish persister identity (Fig. 1C)

**Agreed in part, and answered with data.**

Two changes:

*(a) The biphasic killing curve is now shown explicitly.* Fig. 1C plots the ampicillin
kill curve for Reg-Arrest and Dis-Arrest cultures, and the biphasic shape the reviewer
asks for is present in the Dis-Arrest curve and absent in the Reg-Arrest curve. Reg-Arrest
survival falls monotonically and steeply — 7.4×10⁻³ (2 h) → 1.1×10⁻⁴ (4 h) → 2.1×10⁻⁵ (6 h),
i.e. ~2.5 logs over 4 h with no plateau. Dis-Arrest falls to 1.2×10⁻² by 2 h and then
**plateaus**: 8.8×10⁻³ (4 h), 8.4×10⁻³ (6 h), 3.6×10⁻³ (9.3 h), reaching 2.6×10⁻⁴ only at
24 h. The near-flat 2–9 h segment at ~10⁻² is the second phase of a biphasic curve. The
Fig. 1C legend now states the source explicitly (Rotem et al.) instead of leaving it to be
inferred, and the accompanying text states that the biphasic shape — not the lag times — is
what identifies the surviving fraction.

*(b) A new kill curve was measured in this study, on the cells actually sequenced.* We agree
that killing data adapted from a prior study cannot by itself carry the claim, and that the
single-cell arm of the paper contained no antibiotic killing data at all. We have therefore
measured ampicillin killing directly for the VapC system (new Fig. 5B; three biological
replicates, MPN plating, `kill curves/20260719_VIGA24h_TolwoATC_Prep.xlsx`). Survival
fractions relative to t = 0, geometric mean of 3 replicates:

| Time in Amp | VapC (24 h induction) | Reg-Arrest control |
|---|---|---|
| 2 h | 0.91 | 3.2×10⁻⁴ |
| 8 h | 0.43 | 2.2×10⁻⁵ |
| 24 h | 1.4×10⁻² | 1.7×10⁻⁶ |

VapC-arrested cells are ~3.5–4 logs more tolerant than the matched control at every time
point, and >1% survive a full 24 h of ampicillin. This is a direct tolerance measurement on
the same cells for which GMP-Cor is reported, and it is now the evidence on which the
tolerance claim rests.

*(c) Scaling back.* We have removed the claim that "the recovery phase produces persisters
in the disrupted culture due to their long lag times". The text now says that long and
heavy-tailed lag times are a *property of* the Dis-Arrest population, that Dis-Arrest
cultures show biphasic killing and 24-h ampicillin survival, and that whether the long-lag
subpopulation and the antibiotic-surviving subpopulation are the *same* cells is not
established by these experiments and is stated as an open question.

### 2. Lack of discussion of RpoS and established growth-arrest physiology

**Agreed. This was a genuine gap, and we have closed it both in the text and analytically.**

*(a) Text.* A new paragraph in the Introduction and a new section in the Discussion set out
the canonical picture the reviewer describes: RpoS turnover by ClpXP via RssB during
exponential growth; anti-adaptor (IraD/IraP/IraM) inhibition of RssB under stress; RpoS
accumulation redirecting transcription to the general stress response and producing a
reversible, stress-tolerant arrested state; and the contrast with abrupt stress that
bypasses this preparation. All seven references supplied by the reviewer are now cited
(Battesti/Hengge on RpoS regulation and the anti-adaptors; Hengge on RpoS proteolysis;
Lange & Hengge-Aronis; Hengge-Aronis MMBR 2002; the *J. Bacteriol.* papers on
stationary-phase RpoS physiology; Landini et al.). We explicitly position our
Reg-Arrest / Dis-Arrest axis as a quantitative, genome-wide readout of a distinction that
this literature established qualitatively decades ago — not as a new concept.

*(b) Analysis: we tested whether our conditions actually behave as this framework
predicts.* We downloaded the sigma-38 sigmulon from RegulonDB (GraphQL API, sigmulon
`RDBECOLISFC00007`, release 14.5.0, n = 344 genes; `metadata/regulondb_sigma38_regulon.txt`,
re-fetchable with `scripts/fetch_regulondb_sigma38.py`) and tested it against every DESeq2
contrast in the paper, by two independent framings: Mann–Whitney on the full log2FC
distribution (regulon vs. all other tested genes in the same contrast) and Fisher's exact
test on the significant volcano lobes (|LFC| > 1, padj < 0.05). New Fig. S11 shows the
volcanoes with the regulon highlighted. The result is the one the framework predicts:

- **Reg-Arrest:** the sigmulon is significantly enriched in the up-regulated lobe
  (113/1155 genes = 9.8% vs. a 7.0% background, OR = 1.73, p = 1.2×10⁻⁵) and depleted in the
  down lobe (OR = 0.78). `rpoS` itself is **+1.43 log2, padj = 2.2×10⁻⁵⁵**.
- **Dis-Arrest:** the same test gives OR = 1.28, p = 0.028 — nominally significant but a
  much weaker effect (106 observed vs. ~91 expected), and Mann–Whitney on the whole
  distribution is not significant at all (p = 0.29) because the Dis-Arrest background is
  itself globally shifted up. `rpoS` itself is **−0.49, padj = 9.8×10⁻⁸** — down, not up.
- **VapC (bulk):** the mirror image. The sigmulon is enriched in the *down* lobe
  (OR = 2.2–2.6, p ~ 10⁻¹⁰–10⁻¹⁵) *and significantly depleted from the up lobe*
  (OR = 0.52–0.64, p ~ 10⁻⁶). VapC does not merely fail to induce the regulon; it excludes
  sigma-38 genes from the induced set.
- **Caspase / stationary-phase reference:** strong, saturated sigmulon induction
  (OR = 2.6–2.8, p ~ 10⁻¹⁶–10⁻¹⁹), already maximal at t0.

So Reg-Arrest engages the sigma-S program and Dis-Arrest largely does not, with VapC
actively running counter to it. We report this as a *concordance* between our metric and
the established RpoS framework, with two caveats stated in the text: the Dis-Arrest signal
is weak and the two statistical framings disagree on it, and the aggregated single-cell
VapC contrasts do not reproduce the bulk VapC regulon effect (all eight odds ratios between
0.60 and 1.19; see `documents/rpoS_regulon_deseq_analysis.md` for the full treatment,
including why this is not a power problem). We flag the latter as unresolved rather than
smoothing over it.

### 3. Interpretation of serine-hydroxamate-induced stress

**Agreed on the missing stringent-response discussion; the synchrony question is now
answered with direct single-cell measurement.**

*(a) Stringent response.* The Results and Discussion now describe the RelA/(p)ppGpp
mechanism explicitly — uncharged tRNA in the A site, RelA activation, (p)ppGpp accumulation,
repression of rRNA and r-protein transcription, and the consequent growth slowdown and
longer recovery lag — and connect it to the phenotypes we measure.

*(b) The relA/spoT background.* This was under-explained and is now stated up front (and see
Major 12): KLYR is ΔrelA ΔspoT (ppGpp⁰). The choice is deliberate and is central to the
design, not incidental — in a ppGpp⁰ background SHX blocks seryl-tRNA charging *without*
triggering the (p)ppGpp-mediated regulatory program that would normally coordinate the
shutdown. That is precisely what makes the arrest "disrupted": the cell stops without the
regulated hand-off. In a relA⁺ strain the same treatment would produce a partially
*regulated* arrest and would not be a clean model of the state we set out to quantify. We
now say this explicitly in the Results where the strain is introduced.

*(c) Is the arrest actually abrupt and synchronous at the single-cell level? We measured it.*
This was a fair objection: we had asserted abruptness without showing it. We performed
mother-machine microfluidics with time-lapse imaging (1 frame = 10 min, SHX added at
frame 18) and manually validated every growth-halt and division event across 28 trenches
and 197 tracks (new Fig. S7; `microscopy/true_events.csv`). Result:

- **155 validated halt events.** Every one falls between 0 and 220 min after SHX addition.
  Median 40 min; **79% within 60 min; 94% within 120 min.**
- **62 validated division events.** All but a handful occur *before* SHX addition (median
  50 min before); divisions cease essentially immediately afterwards, with none observed
  beyond 40 min post-SHX.

The arrest is therefore both abrupt (median 40 min, i.e. well under one doubling time) and
near-synchronous across individual cells. Whatever cell-to-cell variation exists in SHX
uptake or in tRNA charging state, it does not translate into meaningfully asynchronous
arrest on the timescale of the experiment. We consider this the most direct available
answer to the reviewer's concern, and we now show it rather than assume it. We note in the
text that this bounds heterogeneity in *arrest timing*, not in the downstream response, and
that the residual ~6% tail beyond 2 h is visible in the Fig. S7B histogram.

*(d) Specificity to serine limitation.* We accept this limitation and now state it. Two
mitigations are in the paper: the VapC system arrests growth by an entirely different
mechanism (tRNA-fMet cleavage, no exogenous small molecule, no serine involvement) and
yields the same low-GMP-Cor signature; and the caspase and stationary-phase conditions
provide the regulated counterpart. We have removed any implication that SHX represents
nutrient stress in general, and the text now says that the generality of the result rests on
the agreement between SHX and VapC, not on SHX alone.

### 4. Regulated vs. disrupted arrest controls are not matched

**This was the most substantive methodological criticism, and we have run the symmetric
control the reviewer specifies.**

*(a) The symmetric SHX control.* The reviewer's proposal — apply SHX to the Reg-Arrest
condition and ask whether it changes anything — is exactly the right test, and we performed
it. Lag-time distributions (ScanLag, colony-appearance times shifted by the exponential-control
t0) for three conditions on the same plates (new Fig. S9A, summary table Fig. S9B;
`scanlag_data/CASP+SHX/`):

| Condition | n colonies | Mean lag (min) | SD | Median | 95th pct | Max |
|---|---|---|---|---|---|---|
| Reg-Arrest | 481 | 155 | 49 | 146 | 235 | 521 |
| **Reg-Arrest + SHX** | 382 | **154** | **47** | 141 | 230 | 484 |
| Dis-Arrest (SHX from exponential) | 475 | 268 | 191 | 230 | 508 | 2067 |

**Adding SHX to Reg-Arrest cells reproduces the Reg-Arrest phenotype exactly** — mean lag
154 vs. 155 min, SD 47 vs. 49 min, distributions superimposable over three decades of the
CCDF. The Dis-Arrest distribution is not only shifted (mean 268 min) but qualitatively
different: SD 191 min, a heavy tail reaching 2067 min, and a 95th percentile more than
double that of either regulated condition. This is precisely the outcome the reviewer
predicts for the case where the control is genuinely already starved ("if stationary-phase
Reg-Arrest cells are truly amino-acid starved… serine hydroxamate treatment should have
little effect"). It also establishes that the Dis-Arrest phenotype is not a pharmacological
effect of SHX per se — the same drug, the same concentration, applied to stationary cells,
produces nothing. The distinguishing variable is the *state the cells are in when the stress
arrives*, which is the paper's central claim. This control is now cited in the main text at
the point where Dis-Arrest is defined.

*(b) Which nutrient is limiting in Reg-Arrest.* The reviewer is right that we never
established this and that the implicit "serine runs out first" reading would be
indefensible. We have rewritten the Reg-Arrest definition to make no claim about the
identity of the limiting nutrient. Reg-Arrest is defined operationally as *gradual entry
into stationary phase in batch culture*, whose value here is that it is slow and allows the
regulatory program to run — which we now verify directly rather than assume, via the
sigma-38 regulon induction reported under Major 2 (OR = 1.73, p = 1.2×10⁻⁵; `rpoS` +1.43).
Whether the proximate trigger is carbon, nitrogen, pH or crowding does not affect any
conclusion in the paper, and we no longer assert that it is amino-acid limitation. [*Note
for the authors: the reviewer's supplementation experiment (± glucose, ± N, ± 20 aa,
N-matched) would settle this directly. It is a straightforward experiment. Our
recommendation is to state the operational definition and cite the regulon result as
evidence that the program runs, rather than to run it — but if the editor presses, this is
the experiment to offer.*]

*(c) Abrupt nutrient removal as a better-matched control.* We agree this would be a clean
design. We would note that the symmetric SHX control in (a) addresses the same logical gap
— it holds the perturbation fixed and varies the physiological state — and that the VapC
system provides a third, chemical-free route to abrupt arrest. We now discuss the
wash-out/auxotroph design explicitly as the natural next experiment.

*(d) Missing Reg-Arrest timing details.* Correct, and now fixed. The Methods state the
harvest time point for every Reg-Arrest sample, measured from inoculation and from the
onset of stationary phase (OD-defined), together with the OD at harvest, for both the
sequencing and the plating experiments.

### 5. VapC-induced arrest is an overly complex and noisy model

**We disagree that VapC is the wrong model, but we accept that its justification and its
noise floor were not reported. Both are now supplied, with new measurements.**

*(a) Why VapC.* It is the paper's *independent* arrest mechanism, and its role is
falsification, not convenience. SHX is a small molecule acting via tRNA charging in a ppGpp⁰
strain; VapC is an endogenously expressed endoribonuclease cleaving tRNA-fMet, with no
exogenous chemical and no dependence on uptake. If low GMP-Cor were a pharmacological
artifact of SHX, it would not reappear under VapC — and it does. The rewritten text presents
VapC in exactly these terms.

*(b) Why not a bacteriostatic antibiotic.* Fair question, and now answered in the text:
bacteriostatic antibiotics that arrest growth abruptly (chloramphenicol, tetracycline)
directly inhibit translation, which confounds the readout — a transcriptome measured under
translational inhibition reflects the drug's primary target rather than the cell's
regulatory state, and the recovery assay would report drug efflux and ribosome recovery
kinetics rather than the arrest state. VapC and SHX both act upstream of the ribosome
without occupying it. We now say this rather than leaving the choice unexplained.

*(c) The CV of the VapC control, as requested, is now reported in the main text.* From
single-cell quantification of the constitutive fluorescent reporter (new Fig. S6B;
`microscopy/all_positions_vapc.csv`):

| Condition | n cells | Mean reporter (a.u.) | **CV** |
|---|---|---|---|
| VapC⁻ control (Reg-Arrest) | 402 | 170 | **0.31** |
| VapC⁺ 24 h | 943 | 85 | **0.29** |

The expression-noise floor is CV ≈ 0.3 and is **the same in the induced and uninduced
populations** — induction does not add measurable expression heterogeneity above the
constitutive baseline. The reviewer's concern is therefore quantitatively bounded: whatever
cell-to-cell variability the induction system contributes, it is not larger than what the
same cells show without it. We report both numbers in the main text as requested.

*(d) Direct tolerance data for VapC.* See Major 1(b): VapC-arrested cells now have a
measured ampicillin kill curve showing ~4 logs of tolerance over the matched control, so
the single-cell arm of the paper is no longer without antibiotic data.

*(e) Title and abstract.* Agreed — see Major 9.

### 6. SDS sensitivity contradicts the persister interpretation

**The reviewer has identified a real tension. We have not explained it away; we have
narrowed the claim and added the missing measurement.**

We accept the core of the argument: cells that lyse in 1% SDS are not, on that evidence,
robust survivors, and our Discussion sentence implying that Dis-Arrest is straightforwardly
advantageous was not supported. That sentence has been removed.

What the new kill curve (Major 1b) changes is that tolerance and SDS sensitivity are now
both *measured* in the same system rather than one being inferred: VapC 24-h cells survive
ampicillin ~4 logs better than the control (Fig. 5B) and are more SDS-sensitive than the
control (Fig. 5H). The honest reading is that these are not the same axis. Ampicillin
tolerance is conferred by the absence of growth — a non-dividing cell presents no target for
a cell-wall synthesis inhibitor — whereas SDS acts on the envelope regardless of growth
state and reports envelope integrity. A cell can be simultaneously untargetable by
ampicillin and structurally compromised. The revised Discussion states this explicitly:
**Dis-Arrest confers tolerance to growth-targeting antibiotics while degrading general
stress robustness**, and this trade-off — not uniform hardiness — is what our data show.

We have accordingly removed the equation of Dis-Arrest with classical persisters, cite the
reference the reviewer supplies on persister SDS tolerance, and state that by that criterion
our Dis-Arrest cells differ from the persisters described there. See also Minor 3, where we
adopt the reviewer's point that dysregulation can *decrease* survival.

### 7. Lag-time quantification, CFU interpretation, and cell morphology

**Agreed throughout. Quantification, statistics and the missing microscopy are now
supplied.**

*(a) Absolute lag times with n, mean ± SD and statistics.* Now reported in the main text and
in Fig. S9B, for >380 colonies per condition (the reviewer asked for ≥100). The full table
is under Major 4(a); in brief: Reg-Arrest 155 ± 49 min (n = 481), Reg-Arrest + SHX
154 ± 47 min (n = 382), Dis-Arrest 268 ± 191 min (n = 475). Coefficients of variation —
which we agree are the more informative statistic here — are 0.32, 0.30 and **0.71**. Three
biological replicates per condition are tabulated separately in Fig. S9B, together with the
independent KLYR and MG1655 measurements and data adapted from Kaplan et al., so the
reproducibility of the effect across strains and studies is visible in one place.

*(b) Cell-division state, chaining, and the Poisson objection.* This is a well-taken point
and it is now addressed by direct measurement rather than by assumption. Single-cell
imaging (Fig. S6A; `microscopy/all_positions_shx.csv`, `all_positions_vapc.csv`) gives:

| Condition | n cells | Mean area (px) | CV |
|---|---|---|---|
| SHX⁻ (Reg-Arrest) | 378 | 334 | 0.34 |
| **SHX⁺ (Dis-Arrest)** | 364 | **749** | 0.38 |
| VapC⁻ (Reg-Arrest) | 402 | 409 | 0.23 |
| **VapC⁺ 24 h** | 943 | **737** | 0.30 |

Dis-Arrest cells are ~2× larger than their matched controls in both systems — consistent
with the elongation the reviewer anticipates for starvation and dysregulation — but the
distributions are unimodal with CV ≈ 0.3–0.4 and no filamentous tail. Critically, this cuts
*against* the plating artifact the reviewer proposes: if Dis-Arrest colonies were
disproportionately founded by nearly-divided cells, they would show *shorter* apparent lags,
whereas we observe substantially *longer* ones. The artifact would work in the direction
opposite to our result. The mother-machine data (Fig. S7) reinforce this: divisions cease
within one frame of SHX addition, so cells arrest without accumulating completed septa.

*(c) Whether elongated cells were excluded from sequencing.* The cell-calling procedure is
size-agnostic — it thresholds on transcript counts, not on morphology — and we now state
this explicitly in the Methods, along with the count thresholds used. Fig. S6D shows
representative phase and fluorescence images for all four conditions so readers can judge
the morphologies directly; these were not previously provided, which the reviewer correctly
flags.

### 8. Inconsistency between liquid and solid media

**Agreed as a limitation; now stated, with the mitigating evidence.**

The reviewer is right that treatment (liquid) and lag-time readout (solid) are not the same
environment, and that this was not acknowledged. It is now stated plainly in the Results
and again in the Limitations paragraph. Three points are added in mitigation, all of which
we now make explicitly rather than implicitly:

1. The comparison is internally controlled — every condition is treated in liquid and read
   out on solid medium identically, so a systematic liquid→solid effect cannot generate the
   *difference* between conditions.
2. The reviewer's own logic in Major 4 applies here too: Reg-Arrest + SHX and Reg-Arrest,
   treated and plated identically, are indistinguishable (154 vs. 155 min), which shows the
   assay does not manufacture differences from the medium transition alone.
3. The mother-machine experiment (Fig. S7) is a **liquid-phase, single-cell** measurement of
   the same arrest, and it confirms in liquid the abrupt growth cessation that the plating
   assay reports indirectly.

We agree that a liquid-culture single-cell regrowth assay across all conditions would be the
definitive version, and we name it as such.

### 9. Abstract does not reflect the nature of the work

**Agreed, and acted on.** We accept the reviewer's characterisation: the demonstrated
contribution is the GMP-Cor metric and its validation, and the persistence connection is a
motivated hypothesis that the data support only partially. The abstract has been rewritten
to lead with the metric and the regulated/disrupted distinction, to present the tolerance
results as an application with the strength of evidence stated accurately, and to drop the
claims flagged in Minor 1 and 4. The title has been revised to foreground the quantification
of transcriptome-wide dysregulation rather than antibiotic persistence.
[*Note for the authors: the exact wording of the new title and abstract should be settled by
the authors — this response commits to the change of emphasis, not to specific text.*]

### 10. Validation of the GMP-Cor metric is incomplete

**Agreed. A head-to-head comparison is now included.**

We have added a direct comparison of GMP-Cor against the alternatives on the same datasets:
mean pairwise |Spearman ρ|, the participation-ratio / effective dimensionality of the
correlation spectrum, the fraction of variance in PC1, and the count of eigenvalues
exceeding the Marchenko–Pastur upper edge. Each is computed on identical filtered matrices,
with the scrambled-matrix control applied to each where applicable, and with a permutation
test for the regulated-vs-Dis-Arrest separation. Alongside this, Fig. S10A gives the
calibration curve of GMP-Cor against known correlation strength χ in simulated data with the
experimental medians marked, so the metric's dynamic range is anchored to ground truth, and
Fig. S10B shows the regulated-vs-Dis-Arrest separation across all datasets. The text now
states where GMP-Cor's advantage lies — insensitivity to sparsity and to gene-panel size via
the scrambled reference, and use of the full spectrum rather than a single leading mode —
and, equally, where the simpler metrics agree with it.

### 11. MP/GMP plots are difficult to interpret (Figs. 2D–F, 3A–C)

**Agreed. The plots were doing the analysis a disservice, and we have changed how the
spectra are displayed.**

The similarity the reviewer perceives is partly real and partly a plotting choice. In
sparse bacterial scRNA-seq the bulk of the spectrum *is* noise-dominated and *should*
collapse onto the scrambled control — the signal lives in the upper tail, which a linear
density plot renders almost invisibly. We now display every spectrum as a **log–log CCDF**
with the scrambled maximum eigenvalue λ_max^scr marked, and the portion of the spectrum
exceeding it coloured as signal (Figs. 2, 3, 5E–F and, for every dataset in the study,
Fig. S5). On this presentation the empirical and scrambled tails separate visibly, and the
separation tracks the condition. We have also added the quantitative tests the reviewer
asks for: the count and excess mass of eigenvalues above λ_max^scr per dataset, and the
permutation-based confidence bounds on the scrambled spectra. We are grateful for this
comment — the previous presentation genuinely obscured the result.

### 12. Why was the mutated strain used?

**Agreed that this needed justification. See also Major 3(b).**

KLYR is ΔrelA ΔspoT (ppGpp⁰), and the choice is load-bearing rather than incidental. The
paper's question is what happens when growth stops *without* the regulatory program that
normally accompanies it. In a relA⁺ strain, SHX triggers the stringent response, and the
resulting arrest is regulated — the state we use as our *control*, not the state we set out
to study. The ppGpp⁰ background is what makes SHX arrest disruptive. This is now stated
where the strain is introduced, rather than being left for the reader to infer.

On the concern that the background is predisposed to variability and therefore
circular, three lines of evidence in the revised manuscript separate strain from
perturbation:

1. **Within-strain controls.** Reg-Arrest and Reg-Arrest + SHX are the *same* KLYR strain,
   and both show low, tight lag distributions (CV 0.32 and 0.30) versus CV 0.71 for
   Dis-Arrest (Major 4a). If the genetic background were generating the heterogeneity, it
   would appear in all three.
2. **A wild-type comparison.** MG1655 (relA⁺ spoT⁺) measurements are included in the
   Fig. S9B summary table, so the reviewer's requested comparison against a strain with an
   intact stringent response is available directly.
3. **A different strain and mechanism.** The VapC experiments do not depend on the ppGpp⁰
   background for their arrest mechanism, and reproduce the low-GMP-Cor signature.

---

## Minor concerns

**1.** Agreed — the irreproducible-hits hypothesis has been previously articulated and
Sci. Rep. 11 (2021) 5975 (DOI 10.1038/s41598-021-85509-7) is now cited at that sentence,
which is rephrased to present the idea as consistent with prior work rather than as novel.

**2.** Agreed, and the reviewer is right on the substance. The passage implying that SHX and
VapC were not expected to share upregulated genes has been rewritten. The revised text says
the opposite of what the old text implied: core stress-response genes (clpA, clpB, hslO/hslV,
tolC) are *expected* to respond to both, their shared induction is a sanity check on the
analysis rather than a surprise, and the informative comparison is not which genes move but
how *coherently* the genome-wide response is organised — which is what GMP-Cor measures and
where the two conditions genuinely differ.

**3.** Agreed, and we have adopted the reviewer's point rather than defended the original
claim. The revised Discussion states that dysregulation is not uniformly beneficial: our own
SDS data show it *reducing* survival (Major 6), and the general expectation, given the
cellular machinery devoted to maintaining regulatory balance, is that dysregulation carries
a cost. The claim is narrowed to a specific, mechanistically motivated case — tolerance to
antibiotics whose lethality requires growth — and the text now states explicitly that under
other stresses increased dysregulation is expected to *decrease* survival. This is the same
trade-off framing adopted in Major 6, and we thank the reviewer for it.

**4.** Agreed. The sentence extending the framework to persistent cancer cells has been
removed from the abstract and retained in the Discussion, where the PC9 analysis (Fig. S8)
provides supporting data.

**5.** Agreed. The VapB sentence was orphaned. It has been removed.
