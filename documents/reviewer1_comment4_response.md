# Reviewer #1 — Comment 4: drafted response

> Source analysis: `scripts/go_param_sensitivity.py` →
> `results/GO_results/param_sensitivity/`
> (log: `go_param_sensitivity_log.json`). Figure change: `scripts/figures/figure4.py`
> panel B. Supplementary figure/table numbers in **[brackets]** are placeholders to
> assign on insertion.

---

**Reviewer comment (paraphrased):** In Fig. 4A several genes up-regulated in Dis-Arrest
belong to operons other than *flg*/*fli* and could be assigned to GO terms beyond
chemotaxis — is Fig. 4A really a specific regulatory program? Given that chemotaxis is a
broad GO category, how sensitive is the analysis to the GO terms included (min/max
gene-set size)? Are Dis- vs Reg-Arrest enrichment differences more pronounced for larger
gene sets? Also, please replace the numeric GO IDs on the Fig. 4B x-axis with GO term
names and gene counts.

---

## Response (Comment 4)

**Method clarification.** We first clarify that the Fig. 4 analysis is a gene-set
*over-representation* analysis (goatools; Fisher's exact test with Benjamini–Hochberg
FDR), not GSEA. The quantity plotted in Fig. 4B (labeled "Enrichment score") is
−log10(FDR), and there is no ranked enrichment statistic or explicit min/max gene-set-size
parameter in the original pipeline. To address the reviewer's robustness concern in the
equivalent terms, we re-ran the analysis while restricting the tested GO terms to a range
of gene-set-size windows, where "gene-set size" is the number of background genes
annotated to a term.

**Specificity of Fig. 4A.** We annotated every gene shown in Fig. 4A to all of its GO-term
memberships. All 36 genes in the panel are annotated to a motility/chemotaxis/flagellar
GO term (36/36; 34 to chemotaxis, GO:0006935, and the remainder to flagellar motility and
assembly terms). These genes do collectively map to other GO terms, as the reviewer notes
(16 distinct motility terms vs. 46 other terms across the panel), but those additional
annotations do not constitute a competing program: the most frequent non-motility terms
are "protein secretion by the type III secretion system" (the flagellar export apparatus;
9 genes), "signal transduction" (chemotaxis two-component signaling; 8 genes) and generic
"protein transport"/"transmembrane transport" (≤7 genes each), each shared by only a
handful of genes. Fig. 4A therefore reflects a specific, coherent motility/chemotaxis
program rather than an artifact of a single broad category (new Supplementary Table [S1],
`fig4a_gene_go_annotations.csv` / `fig4a_go_membership_summary.csv`).

**Gene-set-size sensitivity.** We tested robustness to the GO terms included by
restricting the shared growth-arrest program (the terms enriched in both conditions, i.e.
the bars in Fig. 4B; n = 20) to a range of gene-set-size windows (minimum 5/10/15 genes;
maximum 100/200/500/∞ genes) and recomputing the Dis- vs Reg-Arrest comparison. The
stronger enrichment in Reg-Arrest is preserved and significant in **all 8 windows tested**
(paired Wilcoxon signed-rank, Reg > Dis, p = 1.6×10⁻³ to 5.4×10⁻³). In parallel, when we
restrict the full set of tested terms to each window and re-apply FDR correction, the core
program terms — chemotaxis, flagellum-dependent motility, flagellar organization,
flagellar assembly, cytoplasmic translation and translation — remain significant
(FDR < 0.05) in both conditions across every window (the only exception is the extreme
"≤100 genes" cap, which mechanically excludes the 101-gene "translation" term). The
conclusion is thus not an artifact of including very broad or very narrow GO categories
(new Supplementary Fig./Table [S2], `go_size_sensitivity_paired.csv`,
`go_size_keyterm_recovery.csv`).

**Are differences larger for bigger gene sets?** Across the 20 terms significant in both
conditions, the *signed* difference in −log10(FDR) between Reg- and Dis-Arrest shows no
significant dependence on gene-set size (Spearman ρ = 0.30, p = 0.2) — i.e., the direction
of the effect (stronger, more coherent enrichment in Reg-Arrest) is independent of how
many genes a term contains. The *absolute magnitude* of the difference does increase mildly
with gene-set size (Spearman ρ = 0.53, p = 0.017), but this is an expected statistical
property of any −log10(FDR) statistic: larger gene sets can reach smaller FDRs and
therefore span a wider dynamic range, so equal proportional differences appear larger on
the −log10 scale. It does not reflect a size-specific bias in the biological conclusion
(new Supplementary Fig. [S3], `go_size_confound.svg`).

**Fig. 4B labels.** We revised the Fig. 4B x-axis to display GO term names together with
their gene counts (e.g., "chemotaxis (n = 34)") in place of the numeric GO identifiers.

---

## Numbers at a glance (for internal reference)

| Item | Value |
|---|---|
| Fig. 4A genes annotated to a motility/chemotaxis/flagellar term | 36 / 36 (100%) |
| Distinct motility vs. other GO terms across Fig. 4A genes | 16 vs. 46 |
| Top non-motility terms (genes) | type III secretion (9), signal transduction (8), protein transport (7) |
| Shared-program terms compared (Fig. 4B) | 20 |
| Size windows with Reg > Dis significant (paired Wilcoxon) | 8 / 8 (p = 1.6e-3 – 5.4e-3) |
| Key program terms significant in both, all windows | 6/6 (5/6 only at the ≤100-gene cap) |
| Signed (Reg−Dis) diff vs gene-set size | Spearman ρ = 0.30, p = 0.2 (n.s.) |
| \|diff\| vs gene-set size | Spearman ρ = 0.53, p = 0.017 |
