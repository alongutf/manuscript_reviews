# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

See @README.md for project overview

Research codebase for scientific paper. The core analysis quantifies gene-gene correlation structure in bacterial scRNA-seq data by comparing empirical eigenvalue spectra against the eigenvalue spectra of a random permutation of the matrix.

## Running Analyses

All analyses run as Jupyter notebooks or Python scripts from the `scripts/` directory:

## Code Architecture

### `src/` — Importable Library

- **`analysis_functions.py`**: Pure numerical routines. Core functions: `get_eig_dist()` computes empirical vs. scrambled eigenvalue distributions; `mp_distribution()` returns the analytic Marchenko-Pastur PDF; `get_entropy()` computes effective dimensionality. Used by `data_functions.py` and notebooks.

- **`data_functions.py`**: Defines the `AnnMat` (annotated matrix) class — the central data structure holding a cell × gene count matrix with row/column names and boolean filter masks. Filters are applied lazily; call `get_filtered_matrix()` to materialize. Also contains `get_annotated_data()` (converts raw probe-count CSV → AnnMat), `plot_eig_dist()`, UMAP embedding, and dataset concatenation.

- **`bulk_functions.py`**: Bulk RNA-seq helpers. `run_deseq()` wraps pydeseq2; `run_go_enrichment()` wraps goatools. **Expects `metadata/` to be a sibling of `os.getcwd()`** — run notebooks from `scripts/` so that `os.path.dirname(os.getcwd())` resolves to the repo root.

- **`reader_functions.py`**: Parses legacy plate reader Excel files (BioTek format) into time-series OD/spectrum DataFrames. Entry points: `get_od_data()`, `get_spectrum_data()`.

- **`tecan_func.py`**: `tecan` class for Tecan-format plate reader files. Supports multi-channel data, background subtraction, spike removal, fluorescence normalization, and derivative calculation.

- **`simulations.py`**: functions for generating the synthetic scRNA-seq data. Additional functions for calculating metrics on the simulated data.
### `scripts/` — Analysis Pipelines

| Notebook                      | Purpose                                                                                           |
|-------------------------------|---------------------------------------------------------------------------------------------------|
| `analysis_notebook.ipynb`     | scRNA-seq processing: probe counts → AnnMat, cell calling, gene filtering, eigenvalue calculation |
| `bulk_analysis.ipynb`         | Bulk RNA-seq: DESeq2 differential expression, GO enrichment                                       |
| `scanpy_analysis.ipynb`       | UMAP, Leiden clustering, marker gene analysis                                                     |
| `supplementary figures.ipynb` | Supplementary figure generation                                                                   |
| `random_matrices.ipynb`       | Synthetic Wishart matrices for GMP model validation                                               |
| `permutation test.ipynb`      | Permutation test for GO term enrichment significance                                              |
| `SDS_experiment.ipynb`        | Plate reader data for SDS sensitivity experiment                                                  |
| `simulated_data.ipynb`        | generating and analyzing the synthetic scRNA-seq data                                             |
| `model_fit.nb`                | Wolfram Mathematica — GMP model fitting (not Python)                                              |

### `scripts/figures/` — Figure Scripts

Each `figureN.py` produces a publication-quality SVG using the `PanelFigure` helper class in `figure_functions.py`. `PanelFigure` uses normalized figure coordinates (`[left, bottom, width, height]`) for precise panel placement. Panel functions (`panel_A`, `panel_B`, …) are defined per figure and called at the bottom of the script.

### Data Flow

```
Raw FASTQ
  → scripts/RNA-seq mapping/single_cell_clean_map_count.bash
  → probe-count CSVs

probe-count CSVs
  → analysis_notebook.ipynb  (AnnMat, filtering, eigenvalues)
  → data_for_paper/*.csv  (filtered cell-gene matrices)
  → data_for_umap/        (matrices for UMAP input)

data_for_paper/*.csv + bulk counts
  → bulk_analysis.ipynb   → results/deseq_results/, results/GO_results/

data_for_umap/ + data_for_paper/
  → scanpy_analysis.ipynb → scanpy/ (UMAP coordinates, marker genes)

data_for_paper/*.csv
  → model fit/pc_data/  (eigenvalue histograms)
  → model_fit.nb        → model fit/*.txt (GMP fit parameters)

scanlag_data/, reader_data/
  → SDS_experiment.ipynb, figure scripts
```

### `simulations/` — Standalone Simulation Runners

Python scripts (not notebooks) that run self-contained simulation experiments and write all outputs automatically.

| Script                             | Purpose                                                    |
|------------------------------------|------------------------------------------------------------|
| `subpopulation_mixing_run.py`      | 50/50 sub-population mixing scenario (Reviewer #1 response)|

### `results/simulation_results/` — Simulation Outputs

All simulation outputs are written here. **Always use these paths when creating new simulation runners.**

```
results/simulation_results/figures/  ← publication figures (.svg, .png)
results/simulation_results/raw/      ← human-readable summaries (.txt)
results/simulation_results/logs/     ← full parameter + results logs (.json)
```

Naming convention: `<experiment_name>_<timestamp>.<ext>`  
Example: `subpopulation_mixing_50_50_20260601_110038.svg`

### Workflow
- When running simulations document every parameter used for reproducibility in a log file
- Simulation figures (SVG, PNG) → `results/simulation_results/figures/`
- Raw/summary text files (.txt) → `results/simulation_results/raw/`
- JSON logs (.json) → `results/simulation_results/logs/`