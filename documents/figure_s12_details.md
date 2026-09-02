# Supplementary Figure S12 — details

**Script:** `scripts/supplementary_figures/figure_s12.py`
**Outputs:** `figure_s12.pdf`, `figure_s12_preview.png` (same directory)
**Addresses:** Reviewer #1, comments 1 and 2.3
**Created:** 2026-08-30

A 2 x 2 grid on a 9 x 7.2 inch canvas. Row 1 (A, B) shows that mixing two
transcriptomically distinct sub-populations *raises* GMP-Cor and never lowers it;
row 2 (C, D) shows how the index scales with the dimensions of the matrix.

| Panel | Rect `[left, bottom, width, height]` | Content |
|---|---|---|
| A | `[0.075, 0.575, 0.37, 0.345]` | Simulated regulated sub-populations (χ = 0.7) |
| B | `[0.565, 0.575, 0.37, 0.345]` | Experimental VapC UMAP with mixing annotation |
| C | `[0.095, 0.085, 0.345, 0.355]` | Cell subsampling, experiment vs simulation |
| D | `[0.595, 0.085, 0.345, 0.355]` | Gene subsampling at fixed cell:gene ratio |

The italic strapline sits at figure coordinates `(0.5, 0.535)`, between the rows.

The earlier simulated *dysregulated* mixture (χ = 0.1, run
`inverted_subpopulation_mixing_dysregulated_20260830_111402.json`; pure 4.08 / 3.26,
mixture 29.68 — 7.3x the larger pure population) was dropped from the figure; the log
is still on disk if it is needed again.

---

## Panel A — simulated regulated sub-populations

UMAP of the 50/50 mixture of two inverted sub-populations at χ = 0.7, from
`inverted_subpopulation_mixing_20260830_105044.json`. χ = 0.7 places a pure
sub-population at GMP-Cor ~44, inside the range observed experimentally for
regulated samples (~30-50).

| | GMP-Cor |
|---|---|
| sub-pop. A | 43.70 |
| sub-pop. B | 44.05 |
| 50/50 mixture | **68.75** — 1.56x the larger pure population |

## Panel B — experimental counterpart

The VapC UMAP reproduced from `scripts/figures/figure5.py` panel C
(`scanpy/umap_coordinates_vapc.csv`). Only the two batches that were actually mixed —
**Exponential** (`C_EXP`) and **VapC-2h** (`C_T2`) — are drawn in colour and annotated;
VapC-5h and VapC-24h are drawn as a transparent grey backdrop (`color='0.6'`,
`alpha=0.25`) so they still give the embedding its shape without competing with the
annotated clusters. The panel carries **no legend** — the two coloured clusters are
named by their annotation boxes.

The two pure populations are the **published per-sample GMP-Cor** from
`results/data_metrics/data_metrics.csv` (column `sum_denoised_ev`) — the same numbers
reported per sample elsewhere in the manuscript. The mixture is computed by
`simulations/dataset_mixing_ratio_run.py` on those same two `data_for_paper` matrices
with the **full union gene panel and all cells** — no subsampling, no Fano cut — under
the paper's own exact-case reporter drop list (run
`dataset_mixing_ratio_20260830_131921`):

```bash
python simulations/dataset_mixing_ratio_run.py     --data-dir data_for_paper     --file1 Expira_biorep_t0A_filtered.csv --file2 VapC_biorep_t2A_filtered.csv     --gene-space union --n-genes 99999 --all-cells     --reporter-handling published --norm-sum 1 --repeats 3 --ratios 0.0 0.5 1.0
```

| | n | p | GMP-Cor | published |
|---|---:|---:|---:|---:|
| Exponential (`Expira_biorep_t0A`) | 994 | 1997 | 30.14 ± 1.10 | **32.72** ± 1.59 |
| VapC-2h (`VapC_biorep_t2A`) | 900 | 1793 | 28.96 ± 0.56 | **27.85** ± 0.89 |
| 50/50 mixture | 1894 | 2649 | **49.57 ± 1.15** | — |

The panel annotates the two **published** values and the computed mixture. Because the
full union panel is used without a Fano cut, each pure population is effectively
evaluated on its own genes again — the other dataset's exclusive genes are all-zero for
those cells and are dropped by `get_eig_dist`, leaving p = 1997 and 1793 — so both pure
points reproduce their published values to within the metric's own scatter. This is the
frame the earlier top-2000-Fano version failed to provide (there Exponential collapsed to
6.80 against a published 32.72).

With `--all-cells` the cell set is fixed, so the ± values above are scramble-draw noise,
not cell-sampling variability.

### What still inflates the mixture

The mixture exceeds both pure populations by 1.64x, but two of the reasons are
arithmetic rather than biological and should be stated with the number:

1. **More genes.** GMP-Cor is extensive in p (panel D). The mixture has 2649 genes
   against Exponential's 1997; scaling Exponential by 2649/1997 alone predicts ~40, i.e.
   most of the gap. Corrected for p, the elevation is roughly 1.2x.
2. **More cells.** n = 1894 against 994 and 900. A larger n lowers the Marchenko-Pastur
   edge (lambda*_scr 5.51 in the mixture against 7.11 and 5.97 in the pure points), which
   raises GMP-Cor for free.
3. **1508 of the 2649 union genes (57%) exist in only one dataset** and are structurally
   zero for every cell of the other, so the perfect separation (group-axis AUC = 1.000,
   separating mode PC1) is partly built in by the pairwise panel construction rather than
   by biology. The runner prints this count and repeats it in its summary.

dGMP for this mixture is 0.43 — 43% of its GMP-Cor disappears when each population is
centered on its own gene means.

The version of this comparison with all three points on one panel at matched n and p is
the `data_for_umap` run `dataset_mixing_ratio_20260830_122400` (Exp 18.66 ± 1.38,
VapC-2h 29.22 ± 4.20, mixture 40.00 ± 4.41, 1.37x, dGMP 0.53). Panel A is the
ground-truth version, where n and p are identical by construction.

### Reporter genes in the published values

Establishing the frame above surfaced a provenance issue in the published per-sample
numbers themselves. The paper's gene-removal list is **exact-case**
(`16s_mature`, `16s_unprocessed`, `LELOBEKK`, `kanR`, `mCherry`), but the
`data_for_paper` matrices do not agree on capitalisation:

- `Expira_biorep_t0A_filtered.csv` stores gene names in lower case, so **nothing in the
  drop list matches** and every reporter/plasmid gene survives: `ampr`, `gfp`, `kanr`,
  `laci`, `mcherry`, `tetr`. `laci` alone is **38.4%** of all counts in that matrix
  (ampr 6.2%, gfp 1.1%, kanr 1.0%).
- `VapC_biorep_t2A_filtered.csv` uses mixed case, so its `kanR` and `mCherry` **were**
  dropped.

The same code therefore filtered the two samples differently, and the published
Exponential value (32.72) includes a plasmid transcript that dominates its library.
Removing the reporters case-insensitively (`--reporter-handling clean`) gives
Exponential 24.96 and a mixture of 34.26. Panel B uses `published` so that it remains
comparable to `data_metrics.csv`, but the cleaner numbers are the defensible ones and
this affects the per-sample table, not just this figure.

The UMAP itself is the published embedding of all four batches and was **not**
recomputed on the 1000-cell subsets used for the GMP-Cor values; it is shown to
identify the populations, not as the input to the calculation.

## Panel C — cell subsampling, complete gene panel

`results/subsampling_experimental/logs/subsampling_experimental_2b_20260830_113604.json`
(key `per_size`). Cells drawn uniformly at random without replacement, 5 seeded draws
per size, the full gene panel retained at every size. Each arm is plotted as a
fraction of **its own** value at n = 1000, because the two arms are not comparable in
absolute units.

- **Exponential (sample_2b)** — `data_for_paper/sample_2b_filtered.csv`, 1041 cells x
  2071 genes, rRNA already removed.
- **Simulation (χ = 0.7)** — matched to that matrix's dimensions *and* sparsity
  (`inv_gamma_scale` 0.04, ~87 detected genes per cell against ~100 in the data).

The two curves agree within their draw-to-draw scatter at every size (0.766 vs 0.780 at
n = 500; 0.913 vs 0.944 at n = 800): real data follows the same scaling as the
simulation. The decline is driven by the Marchenko-Pastur edge rising as cells are
removed (λ*_scr 6.21 → 9.62 from n = 1000 → 500) rather than by any loss of structure —
the gene panel is fixed, p stays at 2058-2070 of 2071, and detection per cell is flat.

**Note on naming.** The panel legend reads "Experimental (regulated)", meaning a
transcriptionally *coordinated* sample as opposed to a dysregulated one. Note that
`sample_2b` is the **exponential-growth** sample (`exp` in `scanpy_shx.h5ad`) and not
the Reg-Arrest condition, which is `sample_13b` / `sample_15b` — the label should not be
read as "Reg-Arrest". See `documents/gmp_cor_provenance_analysis.md` §2.

## Panel D — gene subsampling at a fixed cell:gene ratio

`results/simulation_results/logs/subsampling_robustness_rho09_20260615_114715.json`
(Experiment 1 of `documents/subsampling_analysis.md`). Cells and genes are subsampled
together to hold the cell:gene ratio at 1:2 (200c/400g … 1000c/2000g), which holds the
Marchenko-Pastur noise edge fixed and isolates the dependence on gene count.

- Left axis (blue): raw GMP-Cor, which is **extensive** — a proportional fit through the
  origin gives 0.0257 x genes.
- Right axis (grey): GMP-Cor per gene, invariant to ~7% across the 5x range.

This is the reason every comparison in this figure is made at a matched gene dimension,
and it is simulation-only: the experimental matrices have no spare gene dimension to
subsample at a fixed ratio while keeping enough cells.

---

## Regenerating

```bash
python scripts/supplementary_figures/figure_s12.py
```

Reads only the source files listed above, all of which are committed run outputs — no
simulation is re-run and no matrix is re-analysed, so the figure is reproducible in
seconds. Every number in the figure now traces to a run log; nothing is hard-coded. If a
source run is regenerated, update the timestamped filenames at the top of the script
(`REG_LOG`, `GENE_LOG`, `CELL_LOG`, `MIX_LOG`).

To regenerate the experimental mixing values themselves:

```bash
python simulations/dataset_mixing_ratio_run.py          # defaults match panel B
```

## Known cosmetic limits

- In panel B the VapC-2h annotation sits at the right edge of the point cloud and its
  connector crosses the cloud; the four batches overlap heavily in this embedding, so
  there is no fully clear placement.
- Row 1 panels carry ~34% vertical headroom to make space for the mixture box, so the
  point clouds occupy only the lower two-thirds of each panel.
