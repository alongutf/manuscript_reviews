# Sign-resolved GO enrichment of the leading eigenvectors

**Companion source:** `scripts/eigenvector_go_signed.py`
**Outputs:** `results/eigenvector_analysis/go_signed_pos/`, `go_signed_neg/`,
`go_signed_summary.txt`, `go_signed_comparison.csv`
**Controls:** the existing `|loading|` runs in `results/eigenvector_analysis/go/` (top-50),
`go_top100/`, `go_top200/` (`scripts/eigenvector_analysis.py`, `scripts/eigenvector_go_cutoffs.py`)

## Why

The existing per-mode GO analysis ranks genes by `|loading|`, which pools the two poles of
an eigenvector into one study set. An eigenvector is a *contrast*: genes at opposite ends
anti-correlate with each other. If one pole is a coherent program and the other is not (or
is a different program), pooling dilutes the signal twice over — the study set is twice as
large, and half of it is unrelated to the enriched term.

This run splits each of the top-5 modes into its **top 100 positive** and **top 100
negative** loading genes, enriched separately against the same gene-panel background. The
union is the same 200 genes as the `go_top200/` run, so that run is the exact control.

**Sign convention.** SVD returns `+v` or `-v` arbitrarily. Each mode is oriented so the
pole carrying the larger squared-loading mass is "positive". So *positive* = dominant pole,
*negative* = opposing pole; neither means a direction of expression change, and the label
is comparable across modes only in that sense.

## Result — the sign split matters, and it changes the picture

| | modes with ≥1 significant term (of 90) | datasets (of 18) |
|---|---|---|
| `\|loading\|` top-50 | 9 | 8 |
| `\|loading\|` top-100 | 11 | 9 |
| `\|loading\|` top-200 | 14 | 9 |
| **signed, 100 per pole** | **23** | **13** |

- **22 of 90 modes gain at least one term** that the pooled top-200 run missed; only **2
  modes lose** every term they had (`adam_matrix_filtered` mode 3, `sample_2b_filtered`
  mode 1 — both had terms driven by genes split across the two poles).
- Eleven modes go from *completely blank* under `|loading|` to significantly enriched,
  including whole datasets that had no enriched mode at all: `adam_matrix_filtered2`,
  `sample_13a/13b_filtered`, `sample_15b_filtered`, `SHX_biorep2A`.
- Where both methods find terms, the signed run usually finds **more** of them
  (e.g. `VapC_biorep_t2A` mode 1: 2 → 15 terms; `VapC_biorep_t5B` mode 1: 5 → 11).

**The enrichment is almost always one-sided.** Only **3 of 23** enriched modes have terms
on *both* poles. In 18 of 23, the enriched pole is the *lighter-mass* one — a compact set
of strongly co-varying genes sitting against a broad, unstructured bulk. That is precisely
the configuration `|loading|` ranking destroys: the bulk pole contributes most of the top
`|loading|` genes and carries no program.

**What the programs are.** The recovered terms are dominated by translation/ribosome
biogenesis (cytoplasmic translation 19 modes, translation 19, large/small subunit assembly
16/12, regulation of translation 12), with transcription antitermination / termination /
elongation (6 modes each), response to antibiotic (5), and response to heat (3). The single
clearest new biology is `adam_matrix_filtered2` modes 3–4, invisible at every `|loading|`
cutoff: mode 4 contrasts **flagellar motility and chemotaxis** on one pole against
**translation** on the other, and mode 3 contrasts **pyrimidine and arginine biosynthesis**
against motility/chemotaxis. Those are genuine bipolar contrasts of two coherent programs —
exactly the case the absolute-value ranking cannot represent.

## What this does and does not overturn

It does not overturn the delocalization result. The leading modes still have effective
gene numbers in the hundreds (`summary.csv` is unchanged — eigenvalues and entropies do not
depend on the ranking), and most of the recovered terms are the ribosomal/translation axis
already visible in the pooled run, i.e. growth-rate structure rather than a
condition-specific regulon. Enrichment remains concentrated in the *regulated* category
(14 of 40 modes) over *dis-arrest* (5 of 35), consistent with the paper's claim.

What it does change is the strength of the negative claim. The statement "the leading modes
recover no coherent program" is too strong: with the poles separated, **13 of 18 datasets**
have at least one enriched leading mode, and one dataset yields a clean two-program
contrast. The accurate statement is that the leading modes are delocalized and mostly
report the translation axis, not that they are uninterpretable.

## Caveats

- The split doubles the number of enrichment tests (180 vs 90). FDR is applied *within*
  each study set, as in the original scripts — there is no correction across poles, modes,
  or datasets in either analysis, so the two are comparable to each other but neither is
  corrected globally. The 2 → 15 term jumps should be read as sensitivity, not as 15
  independently established findings.
- Each signed study set is 100 genes vs 200 pooled, so part of the gain is simply a purer
  study set at a smaller size; the top-100 `|loading|` row above (11 modes) is the
  size-matched control and the signed run still beats it roughly two-fold.
- `go_signed_*/` files are written only when a pole has ≥1 significant term, matching the
  convention of the sibling scripts.
- Like its siblings, the script reads condition labels and GMP-Cor from the stale
  `results/data_metrics/test8.csv` for display only; the numbers in the headers do not
  affect any eigenvector or enrichment result.

## Per-sample / per-mode pole contrasts

Full listing: `results/eigenvector_analysis/go_signed_contrast_table.txt` (all 18 samples ×
5 modes, written by `write_contrast_table()` in the companion script).

Of 90 modes, **67 are blank on both poles**, 20 are enriched on one pole only, and **3 are
enriched on both** — those 3 are the only genuine two-program contrasts in the dataset:

| sample | mode | positive pole | negative pole |
|---|---|---|---|
| `adam_matrix_filtered2` | 3 | pyrimidine + L-arginine biosynthesis (de novo UMP) | flagellar motility, chemotaxis, swarming |
| `adam_matrix_filtered2` | 4 | flagellar motility, chemotaxis, flagellum organization | cytoplasmic translation, LSU assembly |
| `SHX_biorep2B` | 1 | lipopolysaccharide biosynthesis | cytoplasmic translation, LSU/SSU assembly, regulation of translation |

Five whole samples are blank at every mode and both poles: `VapC_biorep_tONA_filtered`,
`deb_Ec_CDS_untreated`, `deb_KP_CDS_untreated`, `sample_15a_filtered`, `sample_2b_filtered`
(the last had 3 terms under pooled `|loading|`, all split across the poles). The two `deb_*`
untreated samples and `Expira` modes 1/4 sit at `mass_pos` ≈ 0.99–1.00, i.e. essentially all
loading mass on one pole — those modes are not contrasts at all but near-uniform components,
so a sign split has nothing to separate.

In the remaining 20 one-sided modes the enriched pole is the translation/ribosome axis in
18 of 20; the exceptions are `adam_matrix_filtered` mode 2 (glycolytic process) and
`adam_matrix_filtered2` mode 2 (pyrimidine nucleotide biosynthesis), each a single term.
The recurring picture is therefore *one* coherent program — nearly always translation —
against an unstructured opposing bulk, not two competing programs.
