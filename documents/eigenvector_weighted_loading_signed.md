# Sign-resolved weighted-loading score — same formula without the absolute value

**Companion source:** `scripts/eigenvector_weighted_loading_signed.py`
**Outputs:** `results/eigenvector_analysis/weighted_loading_signed/{scores,go_pos,go_neg}/`,
`summary.csv`, `summary.txt`
**Controls:** `results/eigenvector_analysis/weighted_loading/go/` (absolute-value score,
top 50, `scripts/eigenvector_weighted_loading_go.py`)
**Sibling analysis:** `documents/eigenvector_go_signed.md` (the per-mode version, which works)

## What was run

The published per-gene score is

    score_i = sum_{k=1..10} |v_{k,i}| * lambda_k / P

This run drops the absolute value and keeps everything else identical — same 10 leading
modes, same eigenvalue weights, same `1/P` normalizer:

    signed_i = sum_{k=1..10} v_{k,i} * lambda_k / P

Genes are then split by the sign of `signed_i`, and the most positive and most negative
genes of each sample are enriched separately against that sample's own panel, at two
cutoffs (50 per pole, matching the control run, and 100 per pole).

## Result — this does not work, for a reason worth stating

| | samples with ≥1 term (of 18) |
|---|---|
| `\|loading\|` score, top 50 (control) | 3 |
| signed, top 50 per pole | 4 |
| signed, top 100 per pole | 4 |

No sample is enriched on both poles at either cutoff. The signed split gains essentially
nothing over the absolute-value score — in sharp contrast to the per-mode split, which
roughly doubled the number of enriched modes.

**The reason is that a signed sum across modes is not a well-defined quantity.** Each
eigenvector's overall sign is arbitrary: SVD may return `+v_k` or `-v_k`, and flipping any
one mode changes every gene's score. The absolute-value formula is immune to this — that is
why it was written with `|v|` in the first place. The per-mode analysis is also immune,
because there only the *labels* "positive pole" and "negative pole" are arbitrary while the
partition itself is fixed.

To make the score reproducible the script orients every mode the same way as
`eigenvector_go_signed.py` (heavier squared-mass pole = positive), which removes the
dependence on LAPACK's arbitrary choice but does not make the cross-mode sum canonical. The
script quantifies how much this matters by also scoring under the raw, unoriented SVD signs:

- **median sign agreement between the two conventions: 53%** — indistinguishable from a
  coin flip. Per sample it ranges from **13%** (`deb_KP_CDS_untreated`) to 83%
  (`sample_13b_filtered`).
- Spearman correlation between the two rankings has a median near zero and is **negative in
  7 of 18 samples** (as low as −0.64).

So for most samples, "the 50 most positive genes" is a different gene set under a different
but equally valid sign convention. Any enrichment read off this score is a property of the
convention, not of the data. The four samples that do return terms are exactly the ones
where one mode dominates the sum, which makes the convention question moot for them —
`VapC_biorep_t5B` (79% agreement, 1487 positive vs 311 negative genes, λ₁ = 30) simply
reproduces the control's translation terms on its negative pole.

## Per-sample terms

At top-50 per pole, all enrichment is on the negative pole and all of it is the
translation/ribosome axis already reported by the absolute-value score:

| sample | pole | terms |
|---|---|---|
| `VapC_biorep_t5B_filtered` | negative | cytoplasmic translation, translation, LSU/SSU assembly, regulation of translation, response to antibiotic (6) |
| `adam_matrix_filtered` | negative | cytoplasmic translation, translation, LSU assembly, regulation of translation, ribosome biogenesis (7) |
| `adam_matrix_filtered2` | negative | cytoplasmic translation, translation, LSU/SSU assembly, ribosome biogenesis (9) |
| `sample_15b_filtered` | negative | regulation of cellular response to stress (1) |

The only change at top-100 is that `SHX_biorep4A` gains 3 translation terms on the
*positive* pole while `sample_15b` loses its single term — i.e. the pole a term appears on
is not even stable across cutoffs within one sample, which is the same instability from a
second angle.

Note `VapC_biorep_t5A_filtered` had 2 terms under the absolute-value score and has **none**
under either signed cutoff: its translation genes split across the two poles.

## Recommendation

Report the per-mode sign split (`documents/eigenvector_go_signed.md`), not this one. The
weighted-loading score aggregates across modes, and aggregation across modes is only
meaningful for a sign-invariant quantity — `|v|` as published, or `v²` as in
`coordination_score`. If a signed gene-level read-out is wanted, the defensible version is
the per-mode one, where the eigenvector's own contrast structure is preserved and no
cross-mode sign choice is needed.
