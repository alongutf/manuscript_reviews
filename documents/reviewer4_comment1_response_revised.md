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

## The eigenvectors do carry real biological signal

**1. The leading mode recovers the translation machinery.** In 8 of the tested
modes the top-loading genes are significantly enriched for cytoplasmic
translation, translation, and ribosomal subunit assembly. The genes involved are
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
negative result. The translation program appears in the leading mode of
Dis-Arrest samples (Dis-Arrest 1B, 2 and 3) just as it does in Reg-Arrest samples
and in exponentially growing cells. No other program appears consistently in one
condition and not the other. Ranking genes by their loadings therefore does not
tell us which genes remain regulated and which have lost regulation.

Three further limitations apply:

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
