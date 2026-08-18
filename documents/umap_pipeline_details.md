# SHX/casp UMAP: pipeline parameters and reproduction

**Date:** 2026-08-17
**Scripts:** `scripts/umap_paper_barcodes.py`, `scripts/plot_umap_comparison.py`
**Outputs:** `scanpy/` (coordinates, figures, marker tables) · `scripts/logs/` (JSON logs)

Companion to [`gmp_cor_provenance_analysis.md`](gmp_cor_provenance_analysis.md),
which covers the barcode discrepancy between `data_for_paper/` and the h5ad.

---

## 1. There are two different SHX UMAPs in the repo

They use the same 4995 cells and both have 5 Leiden clusters, so they are easy to
confuse, but their coordinates differ.

| coordinates | h5ad | used by | normalisation |
|---|---|---|---|
| `scanpy/umap_coordinates.csv` | `data/scanpy_shx.h5ad` | **`figures/figure1.py` panels G/H, supplementary fig. s6** | uniform |
| `scanpy/umap_coordinates_shx_scaled.csv` | `data/scanpy_shx_scaled_by_total_rna.h5ad` | nothing in the figure scripts | per-condition 3802 / 3507 / 46 |

**`umap_coordinates.csv` is the published one.** The per-condition scalings
(`exp=3802`, SHX`=3507`, casp`=46`) that appear in `scripts/scanpy_analysis.ipynb`
cell 0 belong to the *other* run. They are hard-coded condition constants,
matching no median or mean of either matrix set, and appear nowhere else in the
repo.

## 2. What `data/scanpy_shx.h5ad` pins down

The h5ad was written by the last cell of `scanpy_analysis.ipynb`, which rebuilds
an AnnData from counts and merges the coordinates back in from CSV. So `uns` and
`obsm` are empty and it records **no** normalisation / PCA / neighbour
parameters. `X` is raw integer counts, identical to `layers['counts']`.

Its matrix (4995 × 4042) *is* reproducible bit-for-bit — same cells, same genes,
same values — from:

- `data_for_umap/<sample>_filtered.csv` barcodes
- gene removal of `['16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR']` —
  **`mCherry` is kept**
- `filter_genes(min_cells=3)`

Adding `mCherry` to the removal list gives 4041 genes, not 4042. The notebook's
current list includes `mCherry`, so it has drifted from what produced the figure.

Note `data_for_umap/*_filtered.csv` is *cell*-filtered only: its gene columns are
identical to `<sample>_unfiltered.csv`.

## 3. Reverse-engineering the rest of the pipeline

66 pipeline variants were scored against `umap_coordinates.csv` by Procrustes
disparity (0 = identical up to translation/rotation/scale).

### The dominant factor is the order of operations

| order | mCherry | HVG | disparity |
|---|---|---|---|
| **remove genes → normalise** | kept | none | **0.054** |
| remove genes → normalise | removed | none | 0.085 |
| remove genes → normalise | removed | 2000 | 0.091 |
| remove genes → normalise | kept | 2000 | 0.113 |
| normalise → remove genes | removed | none | 0.304 |
| normalise → remove genes | kept | none | 0.316 |
| normalise → remove genes | either | 2000 | 0.452 |

Normalising *before* dropping the genes is 4–6× worse in every matched pair.

**Why it matters:** `16s_mature` + `16s_unprocessed` are ~92% of all counts. If
`normalize_total` runs while they are still present, it effectively divides each
cell by its rRNA content, and every gene that survives into the analysis carries
a scale factor that has almost nothing to do with it. The notebook normalises
each sample first and drops the genes afterwards; the published run did the
reverse.

### Secondary factors

- **Normalisation scheme.** Uniform beats per-condition decisively (0.054 vs
  0.298 at best). The published UMAP did **not** use the 3802/3507/46 scalings.
- **HVG subsetting.** Keeping all genes beats `n_top_genes=2000` in every matched
  pair, consistent with the 4042-gene object.
- **mCherry.** Minor (0.054 vs 0.085).

### Closest reproduction found

```python
# on the data_for_umap barcodes
remove ['16s_mature', '16s_unprocessed', 'LELOBEKK', 'kanR']   # mCherry kept
filter_genes(min_cells=3)                                       # -> 4042 genes
normalize_total(target_sum=1e4); log1p()                        # AFTER removal
# no HVG subsetting
scale(max_value=10)
pca(n_comps=50, svd_solver='arpack')
neighbors(n_neighbors=40, n_pcs=40)
umap(min_dist=0.3, random_state=0)
```

Disparity **0.054**; reference cluster sizes are {0: 1522, 1: 1257, 2: 1238,
3: 738, 4: 240}. This is *not* an exact reproduction. Since the saved object
carries no `uns`, a scanpy/umap version difference from when the figure was made
cannot be ruled out.

---

## 4. `scripts/umap_paper_barcodes.py`

Runs the pipeline on either cell set. Genes are removed **before** normalisation.

```bash
python umap_paper_barcodes.py --barcodes paper      # data_for_paper barcodes
python umap_paper_barcodes.py --barcodes umap       # original data_for_umap barcodes
python umap_paper_barcodes.py --barcodes umap --n-top-genes 0 --keep-mcherry
                                                    # published run's configuration
```

### Parameters

| stage | value | notes |
|---|---|---|
| cell set | `--barcodes paper` (default) / `umap` | `paper` = unfiltered matrix subset to `data_for_paper` barcodes |
| matrix source | `rnaseq_correlations/data/<sample>_unfiltered.csv` | for `--barcodes paper`; all paper barcodes are present in it |
| concat | `join='inner'`, `keys=['exp','dis1','dis2','reg1','reg2']`, `index_unique='-'` | order `2b, 13a, 15a, 13b, 15b` |
| genes removed | `16s_mature`, `16s_unprocessed`, `LELOBEKK`, `kanR`, `mCherry` | `--keep-mcherry` retains mCherry |
| gene filter | `filter_genes(min_cells=3)` | |
| `layers['counts']` | raw counts, stored here | |
| normalisation | `normalize_total(target_sum=1e4)` then `log1p()` | `--target-sum`; uniform across samples |
| HVG | `highly_variable_genes(n_top_genes=2000, subset=True)` | `--n-top-genes 0` disables |
| scaling | `scale()` (unclipped) | |
| PCA | `n_comps=50`, `svd_solver='arpack'` | |
| neighbours | `n_neighbors=40`, `n_pcs=40` | |
| UMAP | `min_dist=0.3`, `random_state=0` | |
| clustering | `leiden(resolution=0.4)` | |
| markers | `rank_genes_groups(method='wilcoxon', use_raw=False, layer='counts')` | |
| figure points | `size=4`, `alpha=0.35`, shuffled draw order | |

### Two behaviour changes vs. the notebook

1. **Gene removal precedes normalisation** (section 3).
2. **`layers['counts']` holds genuine raw counts.** The notebook assigns it
   *after* normalisation, so `rank_genes_groups(..., layer='counts')` was running
   Wilcoxon on log-normalised values labelled as counts. Marker tables from this
   script are on real counts.

### Outputs

`<src>` = `paper_barcodes` | `umap_barcodes`, `<ts>` = `YYYYmmdd_HHMMSS`.

```
scanpy/umap_coordinates_shx_<src>_<ts>.csv     UMAP_1, UMAP_2, batch, cluster
scanpy/umap_shx_<src>_<ts>.svg / .png          two panels: clusters, sample
scanpy/marker_genes_shx_<src>_<ts>.xlsx        one sheet per cluster
scripts/logs/umap_<src>_<ts>.json              every parameter + per-sample counts
```

---

## 5. `scripts/plot_umap_comparison.py`

Draws three runs in one format from saved coordinate CSVs, so any visual
difference is a difference in the data rather than in the plotting. Row 1 is the
published embedding; the other two are the most recent output of
`umap_paper_barcodes.py` for each cell set (picked up by glob).

Points are `size=4`, `alpha=0.35`, and all categories go into a **single scatter
in a fixed shuffled order** (`default_rng(0)`). Drawing one category at a time
paints the last one entirely over the others — in the earlier version `reg2`
covered whatever it overlapped, which made the samples look more separated than
they are. Colours are scanpy's `default_20` in pipeline category order, so
`exp/dis1/dis2/reg1/reg2` are blue/orange/green/red/purple in every row. Cluster
colours are **not** comparable across rows: Leiden IDs are arbitrary.

```
scanpy/umap_shx_published_<ts>.svg / .png            row 1 alone
scanpy/umap_shx_three_way_comparison_<ts>.svg / .png 3x2 grid
```

### Similarity to the published embedding

| run | Procrustes | ARI | clusters |
|---|---|---|---|
| original barcodes, ordering fixed | 0.091 | 0.591 | 4 |
| paper barcodes, ordering fixed | 0.178 | 0.488 | 3 |
| original barcodes, normalise-first | 0.452 | 0.517 | 3 |
| paper barcodes, normalise-first | 0.385 | 0.515 | 3 |
| `umap_coordinates_shx_scaled.csv` (per-condition run) | 0.448 | 0.413 | 5 |

The remaining 0.091 for the fixed-ordering run reflects the defaults still using
2000 HVGs and removing mCherry; `--n-top-genes 0 --keep-mcherry` reaches 0.054.

---

## 6. Interpretation

Cluster composition by sample (row fractions):

| | published | orig barcodes | paper barcodes |
|---|---|---|---|
| SHX (dis1/dis2) separate | yes (cl. 1+3: 82–97%) | yes (cl. 2+3: 79–96%) | yes (cl. 2: 30–50%) |
| exp vs casp separable | weakly (exp 54/27 vs reg 47/41) | no (exp 49/49 vs reg ~53/40) | no (exp 46/52 vs reg 57/40) |

- **The SHX-vs-rest separation is robust.** It survives the ordering fix and the
  barcode swap.
- **The exponential-vs-natural-starvation distinction is not.** In the published
  run exp and reg differ only modestly in cluster proportions, and once genes are
  removed before normalisation they become indistinguishable. That part was
  carried by the rRNA-driven normalisation.

## 7. Open items

- No exact reproduction of `umap_coordinates.csv` (best 0.054). Candidate causes:
  scanpy/umap version drift, or an unswept parameter.
- The notebook's gene-removal list (`mCherry` included) and its normalisation
  order both differ from what produced the published figure. If the figure is
  regenerated, the notebook cell should be reconciled with section 3 first.
