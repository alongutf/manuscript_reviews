# Point-by-point suggested text edits to `main.tex`

Prepared against `documents/Rebuttal20082026.docx` and the current `documents/main.tex`
(line numbers refer to the current file). **No edits have been made to the LaTeX source.**

## Status summary

The rebuttal promises a substantially revised manuscript. Comparing it against the current
`main.tex`, the following have **already landed**: the corrected Zhu et al. summary (l. 69), the
Jensen et al. Discussion sentence (l. 106), the new sum-above-threshold GMP-Cor definition
(Methods, l. 163–166), the permutation test (Methods, l. 170–177), the CCDF presentation of the
spectra (l. 84), the Ma et al. external datasets (l. 84), the PC9/Oren analysis in the Discussion
(l. 115), the eigenvector supplementary note (l. 86), and the microscopy/software methods
(l. 208–213).

The following are promised in the rebuttal but are **not yet in the main text** and are the
substance of this document: the title, the abstract, the killing curves (Fig. 1C, Fig. 5B), the
RpoS framework and its references, the stringent-response/strain justification, the symmetric
SHX control, the calibration curve and expected GMP-Cor ranges, the subsampling/extensivity
caveat, the pseudobulk GO analysis, the scRNA-seq-vs-bulk agreement, the mother-machine
synchrony data, the VapC CV and cell morphology, the liquid-vs-solid caveat, the SDS/tolerance
trade-off rewrite, and the removal of the VapB sentence.

Three global consistency problems are flagged in **Section G** at the end — please read those
first, since two of them affect text you would otherwise write twice.

---

# A. Title and Abstract

### A1. Title — Reviewer #3 (opening remarks, and comment 9)
**Current (l. 31):** `Genome-Wide Dysregulation in Antibiotic Tolerance and Persistence`

**Rebuttal commitment:** "The title has been revised to foreground the quantification of
transcriptome-wide dysregulation."

**Suggested replacement:**
> `Quantifying Transcriptome-Wide Dysregulation in Disrupted Bacterial Growth-Arrest`

Alternatives, if you want tolerance to remain visible:
> `A Correlation-Spectrum Metric Quantifies Transcriptome-Wide Dysregulation in Growth-Arrested Bacteria`
>
> `Quantifying Genome-Wide Dysregulation Distinguishes Regulated from Disrupted Growth-Arrest`

The key change the reviewer asked for is that "persistence" should not be the head noun.

### A2. Abstract — Reviewer #3 comments 9 and minor 4
**Current (l. 51).** Leads with stress-response biology, presents dysregulation as the finding,
and closes with an unsupported extension to cancer persisters.

**Problems to fix:** (i) the metric, not persistence, should lead; (ii) the persistence claim must
be stated at the strength the data now support (kill curves exist, so it can be stated — but as an
application, not the thesis); (iii) Reviewer #3 asked for the cancer sentence to be removed from
the abstract. Since you *now* have the PC9 analysis, you may keep a single, hedged clause, but it
must be clearly marked as a demonstration on external data, not a claim about this system.

**Suggested rewrite (drop-in for l. 51):**
> When cells are exposed to unfavorable conditions, they activate stress-response pathways and
> enter a regulated growth-arrest. When the stress overwhelms these pathways, cells instead enter a
> *disrupted growth-arrest* — a poorly adapted, abrupt cessation of growth \cite{Kaplan2021}. Whether
> these two archetypes of growth-arrest can be distinguished at the molecular level has remained an
> open question, because the distinguishing feature is expected to be global rather than
> pathway-specific. Here we develop GMP-Cor, a metric that quantifies the genome-wide strength of
> gene–gene coordination from a single snapshot of single-cell RNA-seq data. GMP-Cor is derived from
> the eigenvalue spectrum of the gene–gene correlation matrix, measured against a scrambled null that
> absorbs the spurious correlations arising from sparsity and finite sample size, and is calibrated
> against synthetic data with a known ground-truth correlation strength. Applying GMP-Cor to
> *E. coli* scRNA-seq, we find that cells in disrupted growth-arrest — induced either by serine
> hydroxamate or by prolonged expression of the toxin VapC — show a significant loss of gene–gene
> correlation relative to exponentially growing or regulated growth-arrested cells, while a mixture of
> distinct regulated subpopulations does not. Loss of coordination is corroborated at the bulk level by
> a loss of coherence in gene-ontology enrichment. Both disrupted conditions display biphasic killing
> under ampicillin, linking transcriptome-wide dysregulation to antibiotic tolerance and persistence,
> and both are simultaneously more sensitive to membrane-targeting stress — indicating a trade-off
> rather than uniform robustness. GMP-Cor requires no external reference, time series, or prior
> pathway knowledge, and we show it also separates proliferating from arrested cells in human cancer
> single-cell data, suggesting it is applicable wherever single-cell resolution is available.

If the editor pushes back on abstract length, the sentence beginning "GMP-Cor requires no external
reference…" is the one to cut, keeping the final clause.

---

# B. Reviewer #1

### B1. Comment 1 — why single-cell? (mixture vs. dysregulation)
**Where:** end of Introduction (l. 69–70), and Results §2.1 opening (l. 79).

**Add to the Introduction, after the sentence ending "…deep in Dis-Arrest (long after a strong
perturbation)." (l. 69):**
> Measuring the preservation of gene–gene coordination requires single-cell resolution by
> construction: the quantity we define below is a property of the eigenvalue spectrum of the
> gene–gene correlation matrix computed *across cells*, whereas bulk RNA-seq yields a single averaged
> expression vector per condition, from which no correlation spectrum can be defined. Importantly,
> a correlation-based readout also behaves differently from bulk measures of variability when a
> population is heterogeneous: a mixture of two internally regulated subpopulations inflates
> population-level variance and would be read as disorganization by a bulk metric, but *raises*
> the correlation signal, because the mixture itself generates a coordinated mode separating the
> two groups (Supp. Note X, Extended Data Fig. SX). A reduced correlation signal therefore cannot be
> produced by subpopulation structure, and specifically indicates a loss of coordination within cells.

This single paragraph answers Comment 1 and pre-empts Comment 2.3; it is worth placing early
because Reviewer #3 raises a version of the same concern.

**Also add a forward reference in Results §2.1**, after "…we now turn our attention to the network
of gene interactions" (l. 80): a half-sentence noting the mixing control, so a reader who skips the
Introduction still meets it.

### B2. Comments 2.1 / 2.2 — expected ranges and positive control
**Where:** Results §2.1, l. 82 (synthetic data paragraph) and l. 86 (GMP-Cor definition paragraph).

The current l. 82 introduces the synthetic generator but stops at "the low eigenvalue peak
predicted by the GMP distribution is less indicative of true correlations when the data is very
sparse." The calibration sweep is missing entirely from the main text.

**Append to l. 82:**
> Because the correlation strength $\chi$ in this generator is set independently of sparsity and
> dropout, it provides a ground truth against which the metric can be calibrated. Sweeping $\chi$
> from 0 to 1 at fixed sparsity (1000 cells, 2000 genes, matched to the experimental scale; 10
> repeats per value) yields a calibration curve (Fig. XX / Extended Data Fig. S4) in which GMP-Cor
> remains at the noise floor (median $\approx0$) for $\chi\lesssim0.5$ and then rises steeply and
> monotonically. This defines the expected ranges for dysregulated (low), partially regulated
> (intermediate) and regulated (high) states, and shows that the metric tracks genuine coupling
> strength rather than sampling noise or dropout.

**Append to l. 86, after the Mann–Whitney sentence:**
> The calibration curve also anchors the biological scale. Exponentially growing cells, which must
> be tightly coordinated, serve as a biological positive control, and the simulated limits
> $\chi\to1$ and $\chi\to0$ bracket the achievable range. Overlaying the experimental conditions on
> the calibration curve places regulated cells at $\chi\approx0.85$ and Dis-Arrest cells at
> $\chi\approx0.6$ (Fig. XX). The spread of GMP-Cor values within a condition is itself expected:
> in the steep region of the calibration curve, modest differences in coupling strength translate
> into large differences in GMP-Cor.

That last sentence directly answers the reviewer's observation about the "fairly broad spread"
and should not be dropped.

### B3. Comment 2.3 — sensitivity to which cells are included
**Where:** Discussion, l. 109 (the limitations paragraph), plus one sentence in Methods.

**Add to Methods, in "Initial processing of the single cell data" (l. 139), after the sentence on
taking the top 2000 genes by Fano factor:**
> GMP-Cor is an extensive quantity: it is a sum over the eigenvalues of a spectrum whose total mass
> equals the number of genes, so its magnitude scales approximately linearly with the number of
> genes analysed. All comparisons in this study are therefore made between datasets of matched gene
> dimension, as fixed by this filtering step. Where a dimension-free quantity is required, GMP-Cor
> can be divided by the number of genes.

**Add to the Discussion limitations paragraph (l. 109), after "…cannot be determined just from the
sparsity of the data.":**
> The metric is also robust to the number of cells included over the range relevant here. In
> synthetic data with a known correlation structure ($\chi=0.9$), subsampling cells fivefold from
> 1000 to 200 while holding the gene set fixed retains ~63–71\% of the full-size value (Extended
> Data Fig. SX), a decline consistent with the expected statistical cost of fewer observations —
> noisier correlation estimates and a higher noise threshold — rather than with any change in the
> underlying structure.

**On the two Dis-Arrest clusters (Fig. 1H):** the rebuttal states the clustering resolution was
retuned. Confirm that Fig. 1H in the revised figure now shows the retuned clustering, and update
the Methods clustering resolution (l. 139 currently reads `resolution=0.4/0.3`) to the value
actually used. If the clusters persist, add one sentence at l. 80 noting that marker-gene analysis
does not separate them into growing and non-growing cells.

### B4. Comment 2.4 — generalizability and practical use
**Where:** Results §2.1 (l. 84, already contains Ma et al.) — extend; and Discussion (l. 113).

**At l. 84, after "This shift reflects the expected gene-gene correlations present in unperturbed,
exponentially growing bacteria.", add:**
> We note that these external datasets were generated with a fundamentally different capture
> chemistry, library structure and sparsity profile from ProBac-seq. That the expected correlation
> signature is nonetheless recovered establishes that GMP-Cor is not specific to our protocol.
> Because the absolute magnitude of the metric depends on the signal-to-noise ratio and sparsity of
> the technique, conclusions about regulated versus dysregulated states should be drawn by
> comparison within a technique.

**In the Discussion (l. 113 area), add a sentence on practical use** — the reviewer explicitly
asked how the metric would inform therapy, and the manuscript currently gestures at this without
answering:
> Practically, we envisage GMP-Cor as a culture-level biomarker rather than a per-cell diagnostic:
> because dysregulated cultures are enriched for long-lag, tolerant cells, a low GMP-Cor flags
> treatment regimens or strain backgrounds likely to generate persisters, and identifies cultures
> in which the collateral sensitivity of disrupted cells to membrane-targeting agents could be
> exploited.

### B5. Comment 2.5 — Zhu et al., Jensen et al., entropy comparison
Largely **done** (l. 69, l. 106). Two residual gaps:

1. The rebuttal promises an explicit statement of the *advantages and limitations* relative to
   entropy. The current text asserts complementarity but never states the limitation of GMP-Cor.
   **Add to the Discussion limitations paragraph (l. 109):**
   > Conversely, transcriptomic entropy \cite{Zhu2020} requires a defined reference distribution and
   > a time series of perturbation responses, which makes it the natural readout for the acute phase
   > of a stress response but leaves it undefined for a steady-state snapshot such as the one studied
   > here. GMP-Cor requires no external reference but does require single-cell data, which remain
   > less accessible than bulk measurements in bacterial systems, and is extensive in the number of
   > genes analysed. The two metrics therefore view global dysregulation from complementary vantage
   > points.

2. Reviewer #3 comment 10 asks for the same comparison and is answered by the same passage —
   cross-check that you are not writing it twice.

**Reference note:** confirm `Jensen2017` resolves to *Jensen et al., Cell Reports* (2017),
doi:10.1016/j.celrep.2017.07.062, and `Zhu2020` to doi:10.1038/s41467-020-18134-z.

### B6. Comment 3 — GO on the scRNA-seq data, and scRNA-seq vs. bulk agreement
**Where:** Results §2.2 (l. 91), which currently opens by switching to bulk with only the
justification "Due to the higher sensitivity of the bulk measurements".

**Insert after the first sentence of l. 91:**
> Because GO enrichment requires a per-condition differential-expression contrast rather than a
> per-cell quantity, we first confirmed that the single-cell data support the same analysis. We
> aggregated the single-cell counts into per-replicate pseudobulk profiles over the full gene panel
> and applied the identical DESeq2 and GOATOOLS pipeline, using a matched study set of the 500 most
> strongly down-regulated genes per condition. The single-cell data recover the same growth-arrest
> GO terms as the bulk analysis (translation, ribosome biogenesis, flagellar motility, chemotaxis,
> TCA cycle, ATP synthesis), and among the terms significant in both conditions the core
> translation and ribosome-biogenesis program is far more strongly enriched in Reg-Arrest than in
> Dis-Arrest (e.g. "translation", FDR $\approx10^{-37}$ versus $10^{-4}$; Extended Data
> Fig. S9D) — the single-cell counterpart of the loss of coherence reported below.

**And add, either here or in Methods:**
> Across the full gene panel, per-gene mean expression agrees between the two data types
> (Dis-Arrest scRNA-seq versus disrupted bulk: Pearson $r=0.61$, Spearman $\rho=0.56$ over ~3,100
> shared genes), and pseudobulk log$_2$ fold-changes track the bulk DESeq2 log$_2$ fold-changes
> (Dis-Arrest versus control, $r=0.61$ over ~3,500 genes; Reg-Arrest versus control, $r=0.33$;
> Extended Data Fig. S9). Perfect agreement is not expected, since the single-cell protocol relies
> on a pre-synthesized probe set whose composition and capture efficiency differ from bulk RNA-seq.

The honest statement of the probe-set caveat is worth keeping — it also pre-empts the technique-
comparison concern raised in Comment 2.4.

### B7. Comment 4 — Fig. 4A specificity and GO gene-set size
**Where:** Results §2.2, l. 91 (chemotaxis example) and l. 93.

**Add to l. 91, after the chemotaxis example sentence:**
> We show chemotaxis because its breadth makes the effect visually legible; the same qualitative
> contrast holds for any term with a pronounced difference in enrichment score between the
> conditions in Fig. \ref{fig:GO_analysis}B, and the panel is intended to illustrate the concept of
> coherence rather than to claim a chemotaxis-specific program.

**Add to l. 93, after the "all GO terms that have a significantly different enrichment score…"
sentence:**
> This effect is not driven by gene-set size: across the 20 terms significant in both conditions,
> the signed difference in enrichment score between Reg-Arrest and Dis-Arrest shows no significant
> dependence on the number of genes in the term (Spearman $\rho=0.30$, $p=0.2$).

**Figure/caption action (not main text):** Fig. 4B x-axis must show GO term names with gene counts
instead of numeric GO identifiers, per the reviewer's explicit request. Update the caption to match.

### B8. Comment 5 — long lag ≠ low GMP-Cor at 2 h/5 h; persistence assay; comparable lag times
**Where:** Results §2.3, l. 101, and Fig. 1/Fig. 5 text.

The current l. 101 ends with "…suggests that prolonged VapC induction pushes cells into Dis-Arrest,
similarly to the growth-arrest induced by SHX", and never addresses the apparent contradiction the
reviewer identifies.

**Append to l. 101:**
> The 2 h and 5 h time points already show prolonged lag times while retaining a high GMP-Cor
> (Fig. \ref{fig:vapC}F,H). We emphasize that this is not a contradiction but the central
> distinction the metric is designed to draw: GMP-Cor does not predict *whether* cells have long
> lags, but *why* — separating a regulated prolonged arrest from a dysregulated one. Prolonged lag
> and loss of coordination are both consequences of a maladaptive response to acute stress, but
> they need not follow the same time course. In the VapC system in particular, the lag is bounded
> from below by the time required for the toxin to be cleared, independently of the regulatory
> state.

**Add a new paragraph (or extend l. 102–103) reporting the killing curves**, which currently do not
appear in the main text at all despite being new Fig. 1C and Fig. 5B:
> To establish directly that these conditions produce persisters rather than merely slow-growing
> tolerant cells, we measured ampicillin killing curves. Reg-Arrest survival falls monotonically and
> steeply, whereas Dis-Arrest survival falls steeply within 2 h and then plateaus near $10^{-2}$
> through 9 h — the second phase of a biphasic curve (Fig. \ref{fig:preprocessing}C, adapted from
> \cite{Rotem2024}). We performed the same assay on VapC-induced cells and observe biphasic killing
> and significantly enhanced survival after 24 h of induction relative to the regulated control
> (Fig. \ref{fig:vapC}B).

**Lag-time normalization (the reviewer's final sub-question, and an acknowledged error):** confirm
that every lag-time distribution in Figs. 1E, 4D and 5H is now normalized to the fastest-appearing
colony of the exponential control of the same strain, and **add one sentence to the ScanLag
Methods (l. 126)** defining that normalization explicitly, so the panels are stated to be
directly comparable.

### B9. Minor 1 — Fig. 5H legend/opacity
Figure fix only (opacity set to 100%); verify the caption lists exactly three distributions.

### B10. Minor 2 — package versions
Partially done. Present: Cell Ranger 7.0.0, scanpy 1.9.8, PyDESeq2 0.4.4, GOATOOLS 1.6.4, and the
microscopy software block (l. 213). **Missing:** versions for the core numerical analysis (Python,
NumPy, SciPy, pandas, umap-learn, leidenalg), for cutadapt, bowtie2 and htseq-count, and for the
simulation code. **Suggested addition at the end of "Statistical analysis and tests" (l. 219):**
> **Software versions.** Analyses were performed in Python X.Y with NumPy X.Y, SciPy X.Y, pandas
> X.Y, scanpy 1.9.8, umap-learn X.Y, leidenalg X.Y, PyDESeq2 0.4.4 and GOATOOLS 1.6.4. Read
> processing used cutadapt X.Y, bowtie2 X.Y, htseq-count X.Y and Cell Ranger 7.0.0.

---

# C. Reviewer #3

### C1. Comment 1 — biphasic killing curve
Covered by the killing-curve paragraph in **B8**. In addition, **the Fig. 1C legend must state
explicitly** that the killing curve is adapted from Rotem et al. and that the plateau phase is the
persister fraction — the reviewer specifically objected to this being left to inference.

### C2. Comment 2 — RpoS and established growth-arrest physiology (the largest missing block)
**Where:** Introduction, l. 58, and a new short Results passage.

The current l. 58 compresses all of regulated growth-arrest into one clause about the stringent
response and cites only `Potrykus2008`. This is the omission the reviewer calls "a significant
conceptual gap".

**Suggested replacement for the middle of l. 58** (from "For example, when an *E. coli* culture
exhausts amino-acids…" through "…highly tolerant to antibiotic treatments"):
> For example, when an *E. coli* culture gradually exhausts nutrients, it enters stationary phase
> through a well-characterized regulatory program. The central determinant is the alternative sigma
> factor RpoS ($\sigma^S$): during exponential growth RpoS is degraded by ClpXP via the RssB
> adaptor, whereas under stress anti-adaptors sequester RssB and stabilize RpoS, which then
> redirects transcription toward stationary-phase and general stress-response genes
> \cite{Hengge2009,Battesti2011,Battesti2013,HenggeAronis2002,Lange1996,Schellhorn1998,Landini2013}.
> Together with the stringent response \cite{Potrykus2008}, this produces a reversible,
> growth-arrested but stress-tolerant state. During such a regulated growth-arrest (Reg-Arrest),
> bacteria are more protected from various forms of stress \cite{storz-hengge2010,Rittershaus2013}
> and highly tolerant to antibiotic treatments \cite{Lewis2007,Balaban2019,Ledger2023}. Crucially,
> this program requires *time*: an abrupt stress can arrest growth before it can run, and the
> resulting arrest is not the same physiological state.

**New citations required (Reviewer #3's own list — all seven should be added to `references.bib`):**

| Suggested key | Reference | DOI |
|---|---|---|
| `Hengge2009` | Hengge, *Res. Microbiol.* (2009) — proteolysis of $\sigma^S$ and the RssB/anti-adaptor system | 10.1016/j.resmic.2009.08.014 |
| `Battesti2013` | Battesti & Gottesman, *Genes Dev.* (2013) — anti-adaptors and $\sigma^S$ stabilization | 10.1101/gad.229617.113 |
| `Battesti2011` | Battesti, Majdalani & Gottesman, *Annu. Rev. Microbiol.* (2011) — the RpoS-mediated general stress response | 10.1146/annurev-micro-090110-102946 |
| `HenggeAronis2002` | Hengge-Aronis, *Microbiol. Mol. Biol. Rev.* (2002) — signal transduction controlling $\sigma^S$ | 10.1128/mmbr.66.3.373-395.2002 |
| `Lange1996` | Lange & Hengge-Aronis, *J. Bacteriol.* (1996) | 10.1128/jb.178.2.470-476.1996 |
| `Landini2013` | Landini et al. (2013) — $\sigma^S$ regulon, in *Subcell. Biochem.* | 10.1007/978-94-007-5940-4_5 |
| `Schellhorn1998` | Schellhorn et al., *J. Bacteriol.* (1998) | 10.1128/jb.180.5.1154-1158.1998 |

**Then add the novelty framing** — the reviewer's real objection is that the manuscript does not
say what is new relative to this literature. **Insert after the Dis-Arrest definition (l. 60):**
> We stress that the regulated arrest and its RpoS-dependent architecture are not in question here.
> What the disrupted-arrest framework adds is a set of quantitative predictions that do not follow
> from the absence of regulation alone: that regulated and disrupted arrest can be distinguished
> without knowing how the disruption was produced; that the duration of the arrest shapes the
> distribution of lag times \cite{Kaplan2021}; that disruption confers antibiotic tolerance while
> creating collateral sensitivity to membrane-targeting agents \cite{Rotem2024}; and — the subject of
> the present work — that gene–gene correlations are lost genome-wide.

**Add a Results passage reporting the RpoS regulon analysis** (currently absent; Extended Data
Fig. S8 is never cited). Place it at the end of §2.2, after the GO time-course:
> To relate our conditions to this established framework directly, we asked whether the RpoS
> regulon behaves as the framework predicts in each of them. We annotated RpoS-regulated genes
> (RegulonDB) in each DESeq2 contrast against exponential growth (Extended Data Fig. S8). Reg-Arrest
> shows a pronounced up-regulation of the regulon, as expected for a controlled entry into
> stationary phase. SHX-induced Dis-Arrest shows only a weak, marginally significant up-regulation,
> consistent with a partially executed program, and VapC-induced arrest — which has no a priori
> relation to the stringent response — shows a marginal *down*-regulation. Disruption is therefore
> not confined to the impairment of one regulatory axis: it extends to cellular stress-response
> programs more broadly.

*(Note: `documents/rpoS_regulon_deseq_analysis.md` and the RegulonDB provenance should be cited in
Methods; add a short "RpoS regulon annotation" Methods paragraph naming the RegulonDB release/date
and the statistical test used.)*

### C3. Comment 3 — the stringent response, and why a *relA spoT* strain
**Where:** Introduction/Results §2.1 (l. 79) and Methods (l. 121).

The mechanism of SHX action is never spelled out, and the strain choice is never justified in the
main text — only in Methods, as a fact.

**Extend l. 79, replacing "In contrast, Dis-Arrest is induced by abruptly exposing cells to serine
hydroxamate (SHX), a chemical inhibitor that blocks serine processing, causing unnatural serine
starvation":**
> In contrast, Dis-Arrest is induced by abruptly exposing cells to serine hydroxamate (SHX), which
> blocks charging of seryl-tRNA and causes an acute, unnatural serine starvation
> \cite{Tosa1971,Traxler2008,Durfee2008}. In a wild-type background, uncharged tRNA entering the
> ribosomal A site activates RelA and triggers (p)ppGpp synthesis, which broadly represses
> rRNA and ribosomal-protein transcription and produces a *regulated* arrest \cite{Potrykus2008}.
> This is precisely the outcome we wish to avoid here. We therefore use *E. coli* KLYR, which
> carries *relA1* and *spoT1* mutations and mounts only a relaxed stringent response: in this
> background SHX arrests growth without the accompanying regulatory program, which is what makes
> the arrest disruptive. Importantly, the strain is not itself the source of the phenotype — the
> same strain administered SHX gradually, or allowed to enter stationary phase naturally, produces
> a regulated arrest with narrow lag-time distributions (Extended Data Fig. S15) — and Dis-Arrest can
> equally be produced by stresses unrelated to the stringent response, such as chloramphenicol or
> sodium azide \cite{Kaplan2021}, and by VapC induction in a *relA*$^+$ *spoT*$^+$ background
> (below).

This one insertion answers Reviewer #3's comments 3 **and** 12 together.

### C4. Comment 3 (continued) — is the SHX arrest abrupt and synchronous?
**Where:** new sentences in Results §2.1 (l. 79), citing the mother-machine data (Extended Data
Fig. S13 — currently uncited).

**Suggested addition:**
> Because the synchrony of the arrest matters for interpreting the population as a single state, we
> measured it directly by time-lapse imaging in a mother machine, scoring every elongation and
> division event (Extended Data Fig. S13, Methods). Of 234 evaluable tracks, essentially none
> elongated or divided beyond 220 min after SHX addition; the median elongation-arrest time was
> 40 min (95th percentile 120 min), on the order of a single doubling time, and of 62 validated
> division events none occurred more than 40 min after SHX addition. Whatever cell-to-cell variation
> exists in SHX uptake or in tRNA charging \cite{Elowitz2007}, it does not produce meaningfully
> asynchronous arrest on the 24 h timescale relevant here. We cannot exclude that such variation has
> downstream consequences that emerge much later, and this is one candidate mechanism for the global
> dysregulation we observe.

**Reference to add:** `Elowitz2007` — the reviewer's cited example of population-level inference of
single-cell variability, doi:10.1126/science.1141967.

### C5. Comment 4 — the Reg-Arrest control is not defined with respect to the limiting nutrient
**Where:** Results §2.1 (l. 79), Introduction (l. 58), and Methods (l. 121).

This is a claim the rebuttal concedes ("we no longer assert that it is amino-acid limitation"), and
the assertion currently survives in **two** places.

1. **Introduction, l. 58:** change "when an *E. coli* culture exhausts amino-acids, it enters a
   'regulated' stationary phase" → "when an *E. coli* culture gradually exhausts nutrients, it
   enters a 'regulated' stationary phase" (this is already handled by the C2 replacement above).
2. **Results, l. 79:** change "For Reg-Arrest, cells are allowed to gradually deplete amino-acids,
   leading to natural starvation and entry into the stationary phase" →
   > For Reg-Arrest, cells are allowed to grow undisturbed into stationary phase through gradual
   > nutrient depletion. We define this condition operationally, by the gradual and controlled
   > manner of entry rather than by the identity of the limiting nutrient: what matters for the
   > comparison drawn here is that the arrest is slow enough for the regulatory program to run,
   > which we verify directly through the RpoS regulon analysis (Extended Data Fig. S8) rather than
   > assume. Whether the proximal trigger is carbon, nitrogen, pH or crowding does not bear on any
   > conclusion below.
3. **Methods, l. 121:** add the exact time at which Reg-Arrest samples were collected (the reviewer
   asked explicitly; the Methods currently say only "growth-arrest for 24 hours"). State the OD and
   the wall-clock time after the culture ceased growing, e.g. *"Reg-Arrest samples were collected
   24 h after entry into stationary phase, defined as X h after inoculation, at OD$_{600}\approx$Y."*

### C6. Comment 4 (continued) — the symmetric SHX control
**Where:** Results §2.1, immediately after the lag-time sentence in l. 79.

Currently absent from the main text; Extended Data Fig. S15 is uncited.

**Suggested addition:**
> A stronger control for the drug itself is to apply SHX symmetrically. Adding SHX to cells already
> in Reg-Arrest reproduces the Reg-Arrest phenotype exactly: the lag-time distributions with and
> without SHX are indistinguishable (mean X min, CV 0.32 versus mean Y min, CV 0.30), whereas
> Dis-Arrest is not merely shifted (mean 268 min) but qualitatively different (SD 191 min, a heavy
> tail extending to 2067 min, and a 95th percentile more than double that of either regulated
> condition; Extended Data Fig. S15, adapted from \cite{Kaplan2021}). The same drug at the same
> concentration therefore produces nothing when applied to cells that are already arrested. The
> distinguishing variable is the physiological state the cells are in when the stress arrives, and
> whether they have time to adapt.

### C7. Comment 5 — why VapC, and why not a bacteriostatic antibiotic
**Where:** Results §2.3, l. 98–99, which currently justifies VapC only by "a well-characterized
stress condition known to induce high levels of antibiotic tolerance".

**Insert after the first sentence of l. 98:**
> We chose VapC for three reasons. First, it arrests growth through a molecule produced inside the
> cell rather than an exogenous chemical, so a low GMP-Cor under VapC cannot be a pharmacological
> artifact of SHX. Second, unlike SHX-induced arrest, VapC-induced arrest had not previously been
> characterized as regulated or disrupted, making it a genuine test of the metric rather than a
> confirmation. Third, bacteriostatic antibiotics that arrest growth equally abruptly
> (chloramphenicol, tetracycline) act directly on translation, which confounds a transcriptome-wide
> readout; both SHX and VapC act upstream of translation and trigger distinct arrest mechanisms.

**Also delete the VapB sentence at l. 99** ("Interestingly, mutations in its antitoxin, VapB, that
result in a very high level of antibiotic tolerance, have frequently emerged during evolutionary
adaptation under cyclic antibiotic exposure \cite{Levin-Reisman2017, Fridman2014}."). The rebuttal
states this has been removed; it has not. If you wish to retain the point, it needs two sentences of
context connecting VapB to the present results — otherwise remove it and, if `Levin-Reisman2017`
and `Fridman2014` are then orphaned, check whether they are cited elsewhere.

**Report the CV of the VapC control** (explicitly requested for the main text; Extended Data
Fig. S12B is uncited). **Add to l. 100:**
> Because induction itself introduces cell-to-cell variability, we quantified the noise floor of the
> system from single-cell measurements of the constitutive fluorescent reporter, giving a
> coefficient of variation of CV = X in the control (Extended Data Fig. S12B) — well below the
> variability in the response we report.

### C8. Comment 6 — SDS sensitivity apparently contradicts persistence
**Where:** Results §2.3 (l. 102–103) and Discussion (l. 113). This is the single most important
interpretive edit in the revision.

**Rewrite the offending Discussion sentence (l. 113).** Current text:
> "This seemingly unfavorable condition, which results in a damaged membrane and suspension of the
> cell's ability to divide, can actually lead to increased survival, due to the prolonged
> growth-arrested state, when a lethal antibiotic treatment targeting growth is applied."

**Suggested replacement:**
> Dis-Arrest is not a straightforwardly advantageous state, and we do not claim that it is. In the
> same cells we measure both enhanced survival of ampicillin (Fig. \ref{fig:vapC}B) and markedly
> increased sensitivity to 1\% SDS (Fig. \ref{fig:vapC}G). These two properties are not correlated,
> and there is no reason they should be: tolerance to ampicillin arises from the absence of growth,
> whereas SDS acts on the envelope irrespective of growth state and reports envelope integrity. A
> cell can therefore be simultaneously tolerant to a growth-targeting antibiotic and structurally
> compromised. What our data describe is a trade-off — survival under antibiotics whose lethality
> requires growth, purchased at the cost of robustness to stresses that act on the envelope
> \cite{Rotem2024}. Persisters need not excel at both, and the expectation that they should stems
> from studies in which tolerance was assayed only against growth-dependent killing
> \cite{BBRC2024}.

**Reference to add:** `BBRC2024` — the reviewer's cited paper reporting persister tolerance to SDS,
doi:10.1016/j.bbrc.2024.150549. Cite it where the expectation is stated, not as a refutation.

**Correspondingly, extend l. 103** (currently: "the SDS assay provides independent evidence that
extended VapC induction drives cells into a Dis-Arrest state"):
> …a Dis-Arrest state. We note that increased SDS sensitivity and increased ampicillin tolerance
> coexist in these cells, and return to the interpretation of this trade-off in the Discussion.

### C9. Comment 7 — absolute lag times, CFU interpretation, morphology
**Where:** Results §2.1 (l. 79), Methods, and figure captions.

1. **Report absolute numbers in the main text**, as requested:
   > Across $n\geq$ 100 colonies per condition, Dis-Arrest cells show both longer and more variable
   > lag times than either regulated condition (Dis-Arrest: mean 268 min, SD 191 min, CV 0.71;
   > Reg-Arrest: mean X min, CV 0.32; Reg-Arrest + SHX: mean Y min, CV 0.30; two-sided
   > Mann–Whitney $U$, $p=$ Z; Extended Data Fig. S15B). Measurements in the wild-type MG1655
   > background (*relA*$^+$ *spoT*$^+$) are included in the same table.
2. **Address the septum/plating artifact directly** (Extended Data Fig. S12A uncited):
   > Single-cell imaging shows that Dis-Arrest cells are approximately twice the size of their
   > matched controls in both systems (Extended Data Fig. S12A), consistent with the elongation
   > expected under starvation and dysregulation. We note this runs opposite to a plating artifact:
   > if Dis-Arrest colonies were disproportionately founded by cells caught mid-division, they would
   > show *shorter* apparent lags, whereas we observe substantially longer ones.
3. **Add to the Methods cell-calling paragraph (l. 139)**, per the rebuttal:
   > The cell-calling procedure is size- and morphology-agnostic: barcodes are retained on the basis
   > of total transcript counts alone (threshold: X UMI), so elongated or filamentous cells are
   > neither preferentially included nor excluded.
4. **Cite Extended Data Fig. S12C** (representative phase and fluorescence images for all four
   conditions) at the point where morphology is first discussed.

### C10. Comment 8 — liquid versus solid media
**Where:** Results §2.1 (l. 79, after the lag-time sentence) and Methods.

**Suggested addition:**
> We note that treatments are applied in liquid culture whereas lag times are read out on solid
> medium. The comparison is internally controlled — every condition is treated in liquid and plated
> identically — so a systematic liquid-to-solid effect cannot generate the differences between
> conditions. All transcriptomic data, which are the primary readout of this study, derive
> exclusively from cells grown and treated in liquid; lag time and CFU measurements are the only
> assays performed on solid medium, and serve as supporting physiological evidence.

Mirror this in one sentence in the ScanLag Methods paragraph (l. 126).

### C11. Comment 10 — head-to-head comparison with existing metrics
Answered by the passage in **B5(1)**. Add one sentence stating why a numerical head-to-head is not
presented, so the reviewer sees the decision rather than an omission:
> We deliberately present this comparison qualitatively rather than as a numerical benchmark: the
> published metrics were developed for bulk perturbation time-series \cite{Zhu2020} or for mammalian
> single-cell data at very different cell and gene numbers, and a direct numerical comparison across
> such different regimes would be more misleading than informative.

### C12. Comment 11 — the spectra are hard to interpret
Largely **done** — the CCDF presentation and the permutation test are in place. Two additions:

1. **Make the rationale explicit at l. 84**, where the CCDF is introduced:
   > We display the spectra as log-log CCDFs because in sparse bacterial scRNA-seq the bulk of the
   > spectrum is genuinely noise-dominated and is *expected* to collapse onto the scrambled control;
   > the informative signal lives in the upper tail, which a linear density plot renders almost
   > invisible.
2. **State the permutation-test result in the text**, not only in Methods: give the empirical
   $p$-value for the original-versus-scrambled comparison for at least the key conditions, and state
   that the confidence intervals shown on GMP-Cor derive from the spread of $\lambda_{max}^{scr}$
   over 2000 permutations.

### C13. Minor 1 — prior articulation of the "no universal persistence gene" hypothesis
**Where:** Discussion, l. 114.

**Add the citation and soften the claim.** Replace the opening of l. 114:
> This understanding implies that screens designed to pinpoint "persistence genes" may yield
> irreproducible hits whenever persistence arises from broad dysregulation…

with:
> This possibility has been raised before \cite{ScientificReports2021}, and our results bear on it
> from a new direction: screens designed to pinpoint "persistence genes" may yield irreproducible
> hits whenever persistence arises from broad dysregulation…

and add, at the end of the paragraph:
> We do not claim that persistence always arises from global dysregulation, nor that the idea of a
> broad, network-level response is new. Two hypotheses — induction of a defined pathway versus a
> global property of the gene network — remain live, with evidence on both sides. Our results
> indicate that in the conditions studied here the evidence favours the latter, and identify
> dysregulation as the specific form that the broad response takes.

**Reference to add:** `ScientificReports2021` — doi:10.1038/s41598-021-85509-7.

### C14. Minor 2 — shared stress-response genes were *expected*
**Where:** Discussion, l. 111. The current framing ("we would not a priori expect each treatment to
share a common transcriptional program") is what the reviewer objects to.

**Suggested replacement for the first three sentences of l. 111:**
> The two Dis-Arrest conditions studied here disrupt growth by fundamentally different means
> (SHX-induced starvation and VapC-induced tRNA cleavage), and accordingly differ in the pathways
> they engage: SHX-induced starvation enhances expression of the *psp* module, whereas VapC drives a
> metabolic shift toward alternative carbon sources and increased biosynthesis \cite{Ronine2025}.
> They do, as expected, share a core stress-response signature — upregulation of the AAA+ chaperones
> *clpA*/*clpB* and of *hslO*/*hslV* \cite{Kanemori1997,Kim2022}, and of factors maintaining envelope
> integrity including the major lipoproteins Lpp and LolB and the TolC efflux channel
> \cite{Mathelié-Guinlet2020,Zhu2024}. Their shared induction is a sanity check on the analysis
> rather than a surprise; it is the *differences* between the two responses that reflect the distinct
> pathways each stress engages. Dis-Arrest is thus not a single phenotype or even a defined "state",
> but a cellular condition, exhibiting dysregulation across treatments that otherwise have little
> in common.

### C15. Minor 3 — dysregulation is not uniformly beneficial
**Where:** Discussion, l. 114 ("many gene deletions may enhance antibiotic persistence by increasing
dysregulation…").

**Append:**
> This expectation is specific to antibiotics whose lethality requires growth. Dysregulation is not
> beneficial in general — our own SDS data show that it reduces survival of envelope-targeting
> stress (Fig. \ref{fig:vapC}G) — and given the extent of cellular machinery devoted to maintaining
> regulatory balance, the default expectation is that dysregulation carries a cost. Under stresses
> that do not depend on growth for their lethality, increased dysregulation should be expected to
> *decrease* survival.

### C16. Minor 4 — cancer statement in the abstract
Handled in **A2**: removed as a standalone claim, retained as a hedged final clause now that the
PC9 analysis exists. The Discussion treatment (l. 115) is already appropriate and needs no change.

### C17. Minor 5 — VapB sentence
Handled in **C7**: delete l. 99's VapB sentence.

---

# D. Reviewer #4

### D1. Comment 1 — what, specifically, is dysregulated? (eigenvectors)
The one-sentence pointer at l. 86 is in place. Given how substantive the analysis is, it is worth
two more sentences in the main text so the negative result is visible without opening the note:

**Extend l. 86:**
> Analyzing the eigenvectors associated with these highly correlated modes shows that the difference
> in dysregulation between Dis-Arrest and Reg-Arrest is not confined to a few genes but involves a
> large fraction of the transcriptome (Supp. Note \ref{supp_eigenvectors}). Ranking genes by the
> above-noise coordinated variance they carry and testing for GO enrichment recovers almost no
> significant terms in any condition; the only recoverable program is translation and ribosomal
> protein synthesis, and it appears in weakly coordinated Dis-Arrest samples as well as strongly
> coordinated regulated ones, reflecting the dominance of highly expressed ribosomal genes in the
> leading mode rather than the regulatory state. Consistently, the effective number of genes
> participating in a mode — obtained from the Shannon entropy of the squared eigenvector components
> — averages $\approx470$ and does not track GMP-Cor (Spearman $\rho\approx0.35$, n.s.). The loss of
> coordination in Dis-Arrest is thus spread across many weakly defined modes, which is precisely why
> a single scalar integrating most of the transcriptome, rather than a list of dysregulated genes,
> is the appropriate observable here.

**Add the caveat to the Discussion limitations paragraph (l. 109):**
> We cannot exclude that condition-specific coordinated programs would become recoverable in deeper
> or higher-coverage bacterial single-cell data.

### D2. Comment 2 — van Opijnen lab work
Already addressed at l. 106 and l. 69. No further edit needed beyond verifying both bibliography
entries resolve (see B5).

---

# E. Consolidated list of references to add to `references.bib`

| Key | Reference | DOI | Cited for |
|---|---|---|---|
| `Hengge2009` | Hengge, *Res. Microbiol.* 2009 | 10.1016/j.resmic.2009.08.014 | RpoS proteolysis (C2) |
| `Battesti2013` | Battesti & Gottesman, *Genes Dev.* 2013 | 10.1101/gad.229617.113 | anti-adaptors (C2) |
| `Battesti2011` | Battesti et al., *Annu. Rev. Microbiol.* 2011 | 10.1146/annurev-micro-090110-102946 | general stress response (C2) |
| `HenggeAronis2002` | Hengge-Aronis, *MMBR* 2002 | 10.1128/mmbr.66.3.373-395.2002 | $\sigma^S$ signal transduction (C2) |
| `Lange1996` | Lange & Hengge-Aronis, *J. Bacteriol.* 1996 | 10.1128/jb.178.2.470-476.1996 | RpoS regulation (C2) |
| `Landini2013` | Landini et al., *Subcell. Biochem.* 2013 | 10.1007/978-94-007-5940-4_5 | $\sigma^S$ regulon (C2) |
| `Schellhorn1998` | Schellhorn et al., *J. Bacteriol.* 1998 | 10.1128/jb.180.5.1154-1158.1998 | stationary-phase regulon (C2) |
| `Elowitz2007` | *Science* 2007 | 10.1126/science.1141967 | inferred single-cell variability in translation (C4) |
| `BBRC2024` | *BBRC* 2024 | 10.1016/j.bbrc.2024.150549 | persister tolerance to SDS (C8) |
| `ScientificReports2021` | *Sci. Rep.* 2021 | 10.1038/s41598-021-85509-7 | prior articulation of the network-level persistence hypothesis (C13) |
| `Jensen2017` | Jensen et al., *Cell Rep.* 2017 | 10.1016/j.celrep.2017.07.062 | **verify key exists** — cited at l. 106 |
| — | RegulonDB (current release + access date) | — | RpoS regulon annotation (C2 Methods) |

Also verify: `Zhu2020` = doi:10.1038/s41467-020-18134-z; `Rotem2024` — the rebuttal refers to this
work as *Rotem et al., Sci. Adv.* **2026**, while the citation key and in-text usage say 2024.
Reconcile the key, year and journal before submission (see G3).

---

# F. Figure and caption actions implied by the text edits

These are not main-text edits but are required for the text above to be citable.

1. **Fig. 1C** — killing curves (adapted from Rotem et al.); legend must state the source and
   identify the plateau as the persister fraction (C1).
2. **Fig. 4B** — x-axis labels: GO term names + gene counts, replacing numeric GO identifiers (B7).
3. **Fig. 5B** — new VapC ampicillin killing curve (B8, C1).
4. **Fig. 5H** — opacity at 100%; caption must list exactly three distributions (B9).
5. **Figs. 1E, 4D, 5H** — all lag-time distributions normalized to the fastest-appearing colony of
   the same-strain exponential control; state this in the captions (B8).
6. **New/renumbered extended data figures currently uncited in `main.tex`:** S4 (calibration curve),
   S8 (RpoS regulon), S9 (pseudobulk GO + scRNA-vs-bulk agreement), S12A–C (cell size, VapC CV,
   representative images), S13 (mother machine), S14 (PC9), S15 (symmetric SHX control + lag-time
   summary table), plus the subsampling figure. Every one needs a `\ref` in the text at the point
   indicated above.
7. **Supplementary notes** — `supplementary_note_simulation.tex` currently carries no `\label`, so
   it cannot be cross-referenced. Add one (e.g. `\label{supp_simulation}`) and cite it at l. 82.

---

# G. Three global consistency problems — resolve before writing

### G1. "GMP-Cor" currently names two different quantities
In the Methods (l. 158) the GMP model parameter is defined as *"The parameter $\chi\in[0,1)$ is
referred to as the GMP Correlation Index (GMP-Cor)"*. Five lines later (l. 163) GMP-Cor is
defined again, differently, as $\sum_{\lambda_i>\lambda_{max}^{scr}}(\lambda_i-\lambda_{max}^{scr})$.
These are not the same quantity and cannot share a name — particularly now that the rebuttal maps
experimental conditions onto $\chi$ values ("regulated cells fall at $\chi\approx0.85$"), which only
makes sense if the two are distinguishable.

**Recommendation:** reserve **GMP-Cor** for the eigenvalue sum (the new definition), and refer to
$\chi$ throughout as the *correlation strength* or *coupling parameter* of the GMP model and of the
synthetic generator. Edit l. 158 accordingly.

### G2. The correlation-strength parameter is called three different things
`main.tex` uses $\chi$; `documents/figure_captions.txt` uses $\rho$ (Fig. 3b, 3e, 3g);
`supplementary_note_simulation.tex` uses $\alpha$ — which in `main.tex` is already taken for the
aspect ratio $P/N$. Unify on $\chi$ everywhere, and check the simulation note especially, where
the collision with $\alpha=P/N$ is actively confusing.

### G3. Figure panel letters in the text do not match the captions
`main.tex` cites the calibration/summary results as Fig. 3E (`\ref{fig:single_cell_diff_cond_results}E`,
l. 86), while `figure_captions.txt` places the calibration curve at **3g** and the Mann–Whitney
comparison at **3h**. The rebuttal refers variously to "new Fig. 3E", "Fig. 2E-G", and
"supplementary fig. s4" for what appears to be the same calibration curve. Fix the figure files
first, then write the panel letters into the text once — otherwise the edits in **B2** will need
redoing.

Related: the rebuttal's answer to Comment 2.3 reports a mixing GMP-Cor of 58.88 with no stated
gene dimension. Since GMP-Cor is extensive in gene number (B3), that number is only interpretable
alongside the dimension it was computed at — state it wherever the value appears.
