# Reviewer #4 — Response to Comment 1 (eigenvectors / gene-level information)

**Comment.** The reviewer asks whether gene-level information ("precisely what is
dysregulated / still regulated") can be recovered by examining the eigenvectors
associated with the unusually large eigenvalues, possibly in comparison to a
theoretical compendium of non-dis-arrest conditions.

---

We thank the reviewer for this thoughtful question, which we investigated directly.

For each condition we computed the eigenvectors (gene loadings) of the gene–gene
correlation matrix, not only the eigenvalues. Individual eigenvectors are, however,
not well defined when eigenvalues are near-degenerate: within a (near-)degenerate
subspace any orthonormal rotation is an equally valid eigenbasis, so the gene
ranking of any single mode is arbitrary and unstable. We therefore summarized the
coordinated structure in a rotation-invariant way. For every gene *i* we computed a
**coordination score**, the diagonal of the de-noised correlation matrix
reconstructed from all modes that exceed the scrambled threshold,

  score_i = Σ_{λ_k > λ_max^scr} (λ_k − λ_max^scr) · v_{k,i}² ,

where v_k is the eigenvector of mode *k*. This score measures how much above-noise
coordinated variance each gene carries; it is stable because it depends on the
signal subspace as a whole rather than on any single (arbitrary) eigenvector. We
then ranked genes by this score and ran GO enrichment (the same goatools pipeline
used for the bulk analysis) on the top-ranked genes against the full gene-panel
background, in every condition.

**Result.** The outcome is informative but, at the level of specific gene programs,
largely negative. In most conditions the top-coordination genes do **not** enrich
for any coherent GO term — and this includes the most strongly coordinated
(highest-GMP-Cor) regulated samples. The single program that is recoverable is the
translation / ribosome-biogenesis module (cytoplasmic translation, ribosomal
subunit assembly), which is driven by the highly expressed, tightly co-varying
ribosomal-protein and tRNA/rRNA probes. Critically, this program appears in **both**
regulated and dis-arrest samples — including weakly coordinated dis-arrest samples —
so it reflects the dominance of the ribosomal probes in the leading mode rather than
the regulatory state. This conclusion is robust to the gene cutoff: re-running the
per-mode enrichment on the top 50, 100, and 200 loading genes recovers the same
translation / ribosome-biogenesis module and no other consistent program (all other
GO terms appear as isolated single hits), and several of the highest-GMP-Cor
regulated samples return no significant terms at any cutoff. Consistent with this, when we correlated the per-gene
coordination scores across conditions (the "compendium" comparison the reviewer
suggests), regulated conditions did not share a common coordinated program any more
than dis-arrest conditions did (mean Spearman ρ between coordination-score vectors:
regulated–regulated 0.19, dis-arrest–dis-arrest 0.39, regulated–dis-arrest 0.28).

To ask the complementary question — *how many genes are dominant in each mode* — we
also measured the localization of the individual eigenvectors directly. Treating the
squared gene loadings of a mode as a probability distribution over genes
(p_i = v_{k,i}²), we computed its Shannon entropy H = −Σ p_i ln p_i and the
corresponding **effective number of participating genes**, exp(H), which ranges from
1 (a single gene carries the entire mode) to *n* (all genes participate equally). The
result is that the leading modes are strongly *delocalized* essentially everywhere:
in all but two of the fifteen conditions the top mode spreads its variance over
hundreds of genes (effective number of dominant genes ≈ 230–980 out of ~1,800–4,000
panel genes), and this holds for the most strongly coordinated (highest-GMP-Cor)
regulated samples just as much as for dis-arrest samples (mean effective number of
genes in the leading mode: regulated ≈ 660, dis-arrest ≈ 500; per-mode entropy is
included in `results/eigenvector_analysis/summary.csv`). GMP-Cor does not track this
localization (Spearman ρ between GMP-Cor and the leading-mode effective gene number
≈ 0.35, not significant). The two exceptions are a single regulated replicate whose
leading modes are localized on ~12–36 genes — top loadings again dominated by the
ribosomal probes rather than a condition-specific program — and the one fully
dysregulated (GMP-Cor ≈ 0) dis-arrest sample, whose leading mode sits at the
scrambled-noise threshold. In short, the coordinated
variance that GMP-Cor integrates is not carried by a small set of dominant genes; it
is distributed across a large fraction of the panel.

**Interpretation.** We read this as a meaningful, and in fact expected, result:
GMP-Cor captures a *global, distributed* property of the gene-correlation network
rather than the activity of a small set of specific gene modules. The loss of
coordination in Dis-Arrest is spread across many weakly defined modes, which is
precisely why a single scalar that integrates the whole eigenvalue spectrum —
rather than a gene-set read-out — is the appropriate measure. This is consistent
with our broader claim that dysregulation is a network-level property and not the
behavior of specific nodes.

**Caveats (now stated explicitly in the text).**
1. Individual eigenvectors are not individually interpretable here because of
   near-degeneracy; any gene-level read-out must be a subspace-level,
   rotation-invariant quantity such as the coordination score above.
2. The probe-based single-cell panel and the modest cell numbers limit per-gene
   resolution. We therefore cannot exclude that recoverable, condition-specific
   programs would emerge in deeper or higher-coverage data.
3. The ribosomal-protein and tRNA/rRNA probes dominate the leading coordinated mode
   and can mask weaker programs.

For these reasons we regard the reviewer's suggested comparison against a
theoretical compendium of non-dis-arrest conditions as a promising direction for
future work — it would require a larger, technically matched reference set than is
currently available — and we have added a sentence to the Discussion to this effect.
