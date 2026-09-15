# Reviewer #4 — Revised Response to Comment 1 (eigenvectors / gene-level information)

**Comment.** The reviewer asks whether gene-level information ("precisely what is
dysregulated / still regulated") can be recovered by examining the eigenvectors
associated with the unusually large eigenvalues, possibly in comparison to a
theoretical compendium of non-dis-arrest conditions.

---

We thank the reviewer for this question. We tested it directly, and the answer is
informative in both directions: the eigenvectors do carry real and interpretable
biology, but they do not separate the two arrest states. We describe both results
below and explain why we keep the scalar GMP-Cor as the main measure.

## What we did

For every sample we computed the gene–gene correlation matrix and kept the five
leading eigenvectors. For each eigenvector we ranked genes by the size of their
loading, took the top-loading genes, and ran GO enrichment against that sample's
own gene panel as background (the same goatools pipeline used for the bulk data).
This is a direct reading of the eigenvectors, without any additional summary
statistic. Results are in `results/eigenvector_analysis/` (`top_genes/` for the
ranked gene lists, `go/` for the enrichment, `summary.csv` for the per-mode
values).

An eigenvector is a *contrast*: genes at its two ends vary in opposite directions
across cells. Ranking by the size of the loading alone pools those two ends into a
single gene list, which can hide a program that occupies only one of them. We
therefore repeated the analysis **separately for each pole** of every mode — the
100 genes with the most positive loading and the 100 with the most negative
loading, enriched independently against the same background
(`scripts/eigenvector_go_signed.py`, results in `go_signed_pos/`, `go_signed_neg/`,
`go_signed_contrast_table.txt`). All numbers quoted below are from this
sign-resolved analysis; the pooled results are retained as a control.

## The eigenvectors do carry real biological signal

**1. The leading mode recovers the translation machinery.** Resolving each mode by
sign, **23 of the 90 leading modes (5 modes × 18 samples) are significantly
enriched, in 13 of the 18 samples** — against 14 modes and 9 samples when genes are
pooled by loading magnitude, so separating the poles roughly doubles what is
recovered. The dominant program is translation: cytoplasmic translation and
translation appear in 19 modes each, ribosomal large- and small-subunit assembly in
16 and 12, and regulation of translation in 12. The genes involved are
ribosomal proteins, elongation factor `fusA`, RNA polymerase subunit `rpoA`, and
tRNA and other stable-RNA probes. This is a coherent, well-defined cellular
program, and it shows that the leading eigenvector is not noise.

**2. In exponentially growing cells, the leading mode separates growth from
stress, by sign.** In the exponential-phase dataset the first eigenvector splits
cleanly into two groups with opposite signs. Among its 30 top-loading genes, 21
are ribosomal proteins and all of them have negative loadings, while 6 are
RpoS-dependent general stress genes (`rpoS`, `dps`, `elaB`, `osmY`, `hdeA`,
`ytjA`) and all of them have positive loadings. Genes with opposite signs in an
eigenvector are anti-correlated across cells. The mode therefore describes a
single axis running from high ribosome content to high stress-response content —
the expected trade-off between growth and stress in a growing population. This is
a real, interpretable biological result obtained directly from an eigenvector.

**2b. The sign structure is meaningful, and in one sample it resolves two
opposing programs.** The sign-resolved analysis above makes this observation
systematic. In 20 of the 23 enriched modes the annotated program occupies a single
pole, the opposite pole carrying no coherent function — a program set against an
unstructured background, which is why pooling the two ends dilutes it. Three modes
carry a distinct program at each end, and are genuine two-programme contrasts:

| sample | mode | one pole | opposite pole |
|---|---|---|---|
| untreated control, dimension-matched subset | 4 | flagellar motility, chemotaxis, flagellum organization | cytoplasmic translation, large-subunit assembly |
| untreated control, dimension-matched subset | 3 | pyrimidine and L-arginine biosynthesis (de novo UMP) | flagellar motility, chemotaxis, swarming |
| Dis-Arrest (SHX biorep 2B) | 1 | lipopolysaccharide biosynthesis | cytoplasmic translation, subunit assembly, regulation of translation |

The motility-versus-translation and motility-versus-nucleotide-biosynthesis axes
are invisible at every cutoff of the pooled analysis and appear only once the poles
are separated. They are consistent with the known trade-off between flagellar
synthesis and biosynthetic/ribosomal investment, and they show that where the data
do contain two anti-correlated programmes, the eigenvector represents both.

**3. Genes transcribed from the same promoter appear together in the same mode,
with the same sign.** Our VapC strains carry a chromosomally integrated cassette
in which `vapC`, `mCherry`, `kanR` and `tetR` are transcribed together. In the
leading mode of these samples all four genes appear among the top-loading genes
and all four carry the same sign (for example +0.13, +0.11, +0.08, +0.07 in
VapC 5h A, and all negative in VapC 24h). Because these genes are known to be
co-transcribed, this acts as an internal positive control: the analysis groups
genuinely co-regulated genes into the same mode, in real data with real dropout.

**4. Simulations confirm the method recovers correlation blocks when they
exist.** We generated synthetic data with a known block structure and matched
sparsity to our experimental data (96.8% zeros, versus a median of 96.6% in our
samples). The leading eigenvector placed 91% of its squared loading in the correct
block, and all 30 of its top-loading genes came from that block. So a failure to
recover a program in the real data is not a failure of sensitivity.

## What the eigenvectors do not do

**They do not distinguish Reg-Arrest from Dis-Arrest.** This is the central
negative result, and the sign-resolved analysis — which recovers roughly twice as
much as the pooled one — does not change it. The translation program appears in
the leading mode of Dis-Arrest samples (Dis-Arrest 1B, 2 and 3) just as it does in
Reg-Arrest samples and in exponentially growing cells: **5 of our 7 Dis-Arrest
samples have at least one enriched leading mode, and in every one the recovered
program is the same translation/ribosome axis** seen in Reg-Arrest (8 of 11
samples). Outside
that axis and the three contrast modes above, only three modes in the entire
dataset return any other program on a pole of their own — glycolysis, pyrimidine
nucleotide biosynthesis and lipopolysaccharide biosynthesis, one term each — and
none of them recurs across replicates of a condition. No program
appears consistently in one arrest state and not the other. Ranking genes by their
loadings, with or without sign resolution, therefore does not tell us which genes
remain regulated and which have lost regulation.

Four further limitations apply:

1. **Only a few programs can be seen at all.** Five modes cannot represent the
   many correlated gene groups present in a cell. In our simulations, where the
   true structure is known, the five leading modes captured only the largest and
   strongest blocks and left the rest invisible, even when the data were perfect.
   Absence of a program from the leading modes is therefore not evidence that the
   program is absent from the cell.

2. **Individual eigenvectors are not stable when eigenvalues are close together.**
   Many of the leading eigenvalues in our data are nearly equal. In simulated
   replicates of the same underlying correlation matrix, only the first two or
   three eigenvectors reproduced (correlation 0.87 and 0.78), while the remaining
   ones did not (0.09–0.32). Gene rankings taken from the lower modes should be
   treated with caution.

3. **Annotation and coverage limit what can be seen.** Between a third and 40% of
   the top-loading genes in each mode are uncharacterised y-genes with no GO
   annotation, so any enrichment test is blind to them. Bacterial single-cell data
   are also sparse and shallow, which limits per-gene resolution.

4. **Sign is informative within a mode, but not across modes.** The two poles of
   one eigenvector are genuinely anti-correlated, which is what the analysis above
   uses. Which pole is called "positive" is arbitrary, however, so signed
   quantities cannot be summed across modes: we verified that a signed
   gene-level score aggregated over the ten leading modes changes the sign of
   about half its genes (median 53% agreement) under an equally valid sign
   convention, and yields no additional enrichment. Gene-level scores that
   aggregate over modes must therefore be built from sign-invariant quantities, as
   ours are.

## Why we keep GMP-Cor as the main measure

The comparison the reviewer suggests — a gene-level read-out compared against a
reference set of well-regulated conditions — would require the eigenvectors to
separate the two states, and they do not. The information that does separate them
is in the eigenvalue spectrum, not in any individual eigenvector: the loss of
coordination in Dis-Arrest is spread across many weak modes rather than
concentrated in a few identifiable gene groups. A scalar that integrates the whole
spectrum is therefore the appropriate measure for our comparison, and a gene list
is not.

We regard this as a meaningful result rather than a limitation of the method.
It supports our broader claim that dysregulation is a property of the whole
correlation network and not the behaviour of a specific set of genes. We have
added a short paragraph to the Discussion stating this, together with the
observation that the leading eigenvector does recover the growth-versus-stress
axis and the translation program, and that a larger and technically matched
reference set would be needed to pursue the compendium comparison the reviewer
describes.

---

## Internal notes (not for the reviewer)

- Source of the numbers in this response: `documents/eigenvector_go_signed.md` and
  `results/eigenvector_analysis/go_signed_comparison.csv` /
  `go_signed_contrast_table.txt` (per-mode sign split, 100 genes per pole, 18
  samples × 5 modes). Category counts use `results/data_metrics/data_metrics.csv`,
  not the stale `test8.csv` the analysis scripts read for their display headers.
- **Dataset naming to confirm before submission.** The two contrast modes in the
  table are from `adam_matrix_filtered2.csv`, which is a strict subset of
  `adam_matrix_filtered.csv` — 998 of 1994 cells and 2007 of 3973 genes, i.e. the
  dimension-matched version. It is called "untreated control, dimension-matched
  subset" above; replace with the label used in the manuscript, and decide whether
  to present a result that appears in the subset but not in the full matrix (modes
  3 and 4 of the full `adam_matrix_filtered` return nothing).
- Limitation 4 summarises `documents/eigenvector_weighted_loading_signed.md`: the
  weighted-loading score without the absolute value (`signed_i = Σ v_{k,i} λ_k / P`)
  gives 4 of 18 samples with terms versus 3 for the published `|v|` score, and its
  gene ranking is not reproducible across sign conventions (median 53% sign
  agreement, Spearman correlation negative in 7 of 18 samples). We do not propose
  using it; it is cited only to justify keeping the score sign-invariant.
- Sample-label check: the SHX -> "Dis-Arrest 1B / 2 / 3" mapping in the body text
  predates this revision and is not recorded anywhere in the repo, so the contrast
  row is labelled "Dis-Arrest (SHX biorep 2B)"; substitute the manuscript's own
  label for that sample.
- `results/eigenvector_analysis/top_genes/` was regenerated for
  `adam_matrix_filtered2.csv`; `deb_Ec_CDS_untreated.csv` and
  `deb_KP_CDS_untreated.csv` are still missing from that folder, and `summary.csv`
  there still covers only the original 15 datasets.
