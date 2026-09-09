# Eigenvector block recovery — validating the gene-level analysis on synthetic data

**Companion source:** `src/eigenvector_validation.py`, `simulations/eigenvector_block_recovery_run.py`
**Outputs:**
- depth-matched (headline): `results/simulation_results/{logs,raw,figures}/eigenvector_block_recovery_20260909_103700.*`
- sparse variant (`inv_gamma_scale` = 0.01): `..._20260909_101852.*`

## Why

`scripts/eigenvector_analysis.py` reads gene-level structure off the top-5 eigenvectors
of the gene–gene correlation matrix, and on the experimental data returns a largely
negative result (Reviewer #4, comment 1; `documents/supplementary_note_eigenvector.tex`):
the leading modes are delocalized (N_eff ≈ 470 of ~2000 genes) and recover no
condition-specific GO program.

That negative result had two incompatible readings, and nothing in the repo
distinguished them:

1. the coordinated signal genuinely is distributed (the paper's claim), or
2. the pipeline cannot recover modular structure at this depth and dropout **even when
   modules exist**.

The synthetic generator settles this, because it builds its correlation matrix out of
known blocks.

## Design

`generate_gram_hub_matrix()` partitions genes into contiguous clusters (Pareto sizes),
gives each cluster a shared loading weight `w_k ~ U(0.1α, α)`, and adds one global hub
factor linking a random subset of per-cluster hub genes. It previously discarded that
structure; it now takes `return_structure=True` and returns
`labels / sizes / strengths / hubs / hub_weights`. The default return value and the RNG
stream are unchanged, so every existing runner is unaffected.

Scoring uses the **raw loadings of the top 5 modes only** — exactly what the
experimental script inspects, no coordination score:

| metric | question |
|---|---|
| `best_block_mass` | fraction of a mode's squared loading mass in its best-matching true block (vs. the mass a random gene set of that size would carry) |
| `effective_n_blocks` | exp(entropy) over blocks: 1 = the mode *is* one block, hundreds = smeared over the panel |
| `top_gene_purity` | fraction of the mode's top-30 loading genes drawn from that block, with a hypergeometric p |
| `projection` | ‖V·u_b‖² of each true block's unit indicator on the 5-mode subspace: 1 = fully captured, 0 = invisible |
| stability | per-mode \|r\| and 5-mode subspace overlap between two independent count replicates of the *same* ground-truth R |

**The noise ladder.** Identical scoring is applied at four layers of one simulation,
which localizes where recovery is lost: `R_oracle` (eigenvectors of the true R — the
5-mode ceiling), `latent` (the correlated Gaussian layer), `true_counts` (NB counts +
library-size scaling), `observed` (after expression-dependent dropout).

Headline condition: α = 0.75, 1000 cells × 2000 genes, dropout rate 1.0, Pareto shape
1.5, hub probability 0.2, σ-seed 31. Ground truth: 853 blocks (554 singletons), 83 of
≥5 genes, largest 77 genes.

**Depth calibration.** `inv_gamma_scale` was raised from 0.01 to **0.03** so the simulated
sparsity matches the experimental data. At 0.01 the observed matrix was 98.7% zeros —
sparser than 16 of the 18 real datasets. At 0.03 it is **96.8% zeros** against a real-data
median of 96.6% (range 87.5–99.0%), with 158 counts and 67 genes per cell, both inside the
real ranges. Zero fraction depends only on the count marginals, so it can be calibrated
without the correlation layer. Note that dropout is a weak sparsity knob here: it moves
zeros only 96.2% → 96.8%, because the inverse-Gamma mean prior already sets the sparsity.

## Result

**The analysis does reconstruct blocks, and dropout is not what limits it.**

- Mode 1 carries **91% of its squared loading mass in the single correct block** in the
  observed (post-dropout) data, and **all 30** of its top-loading genes come from that
  block (p ≈ 7e-46). Modes 2 and 3 map just as cleanly onto the next two blocks (76% and
  67% mass, 30/30 top genes each). At the ceiling mode 1 is at 100%.
- The largest block's projection onto the 5-mode subspace degrades only
  **1.00 → 0.97 → 0.90 → 0.87** across the ladder. Recovery is **not destroyed by dropout
  or count noise**; the whole ladder costs ~13%.
- What limits recovery is **capacity, not noise**: 5 modes cannot represent 83 blocks.
  Only 3 blocks exceed projection 0.5 in the observed data (4 at the ceiling), and the
  median block sits at ~0.002 in *every* stage including the true correlation matrix.
  Blocks are reached in order of size × strength.
- **One of the top-5 modes is the global hub factor, by construction.** In the oracle,
  mode 2 has `effective_n_blocks` = 115 and its best block is a singleton. A delocalized
  leading mode is exactly what a real hub/master-regulator structure produces — it is not
  evidence against modularity.
- **Individual eigenvectors are unstable, as the rebuttal argues.** Across two count
  replicates of the same R: per-mode |r| = 0.87, 0.78, 0.32, 0.09, 0.26, and 5-mode
  subspace overlap 0.38. Only the first two or three modes reproduce; the rest rotate
  freely within the near-degenerate block.

### GMP-Cor on the same matrices

| stage | GMP-Cor | per gene | λ_max | λ_max^scr | modes above | zeros |
|---|---|---|---|---|---|---|
| true R | undefined (no cells to scramble) | — | 54.20 | — | — | — |
| NB counts | 45.0 | 0.0225 | 18.93 | 5.89 | 25 | 96.2% |
| observed | **32.6** | 0.0163 | **15.97** | 5.90 | 24 | 96.8% |

At 32.6 the depth-matched simulation sits squarely inside the Reg-Arrest range of the real
data (13.1–137.3 in `results/data_metrics/data_metrics.csv`), as an α = 0.75 dataset
should. The earlier over-sparse run gave GMP-Cor 8.9, which would have been scored as
Dis-Arrest despite *identical ground-truth coupling* — so **GMP-Cor is strongly confounded
by sequencing depth, and simulated and experimental values are comparable only at matched
sparsity.** (This is the same depth confound already noted for the experimental data in
`documents/gmp_cor_provenance_analysis.md`.)

**GMP-Cor and the eigenvector read-out are limited by different things.** Down the same
ladder GMP-Cor falls 45.0 → 32.6 (−28%) while the largest block's projection falls
0.90 → 0.87 (−3%). Count noise erodes the scalar; the 5-mode capacity limits block
identification.

## What this means for the manuscript

The negative experimental result survives, and the interpretation gets stronger and more
precise. The pipeline demonstrably recovers a 77-gene block through realistic dropout, so
the failure to recover programs in the real data is not a sensitivity failure of the
method — and this now holds at the real data's own sparsity, not only at an idealized
depth — reading #1 above is supported, and caveat #2 of
`documents/reviewer4_comment1_response.md` can be narrowed accordingly.

Two caveats should be stated with it:

1. A top-5 analysis can only ever reach the few largest, strongest modules. "No program
   in the top 5 modes" is not "no modules" — it is "no module large and coherent enough
   to outrank everything else." That is a limit of the read-out, not of the data.
2. The simulation's modules are cleanly separable by construction; real regulons overlap.

Both are limits of the *read-out*, not of dropout, which is a different and more
defensible claim than the current text makes.

## Not yet done

Only the two depth conditions above have been run at a single α. The natural extensions are an α sweep
(0.1–0.9) to find where recovery breaks, a depth sweep across the full real-data sparsity
range (87.5–99.0% zeros) to bound caveat 2 of the reviewer response quantitatively — the
two points so far (98.7% and 96.8% zeros) already move GMP-Cor 8.9 → 32.6 at fixed α — and
a negative control at α → 0 where every metric should collapse to the random-gene-set
null.
