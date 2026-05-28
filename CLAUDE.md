# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research codebase for the paper **"Genome-wide Dysregulation in Antibiotic Tolerance and Persistence"**. The core analysis quantifies gene-gene correlation structure in bacterial scRNA-seq data by comparing empirical eigenvalue spectra against the Generalized Marchenko-Pastur (GMP) distribution — a theoretical null model for uncorrelated genes.

## Environment Setup

```bash
conda env create -f environment.yml --name rnaseq_correlations
conda activate rnaseq_correlations
pip install -e .   # installs src/ as an importable package
```

The conda environment is named `rnaseq_correlations` (not `rnaseq_env` as the README says).

## Running Analyses

All analyses run as Jupyter notebooks or Python scripts from the `scripts/` directory:

```bash
cd scripts
jupyter lab
```

Figure scripts must be run from `scripts/figures/` (they use `os.getcwd()` as `root_dir` and look for data relative to that):

```bash
cd scripts/figures
python figure1.py   # saves figure1.svg in the current directory
```

There is no test suite or linter configured.

## Code Architecture

### `src/` — Importable Library

- **`analysis_functions.py`**: Pure numerical routines. Core functions: `get_eig_dist()` computes empirical vs. scrambled eigenvalue distributions; `mp_distribution()` returns the analytic Marchenko-Pastur PDF; `get_entropy()` computes effective dimensionality. Used by `data_functions.py` and notebooks.

- **`data_functions.py`**: Defines the `AnnMat` (annotated matrix) class — the central data structure holding a cell × gene count matrix with row/column names and boolean filter masks. Filters are applied lazily; call `get_filtered_matrix()` to materialize. Also contains `get_annotated_data()` (converts raw probe-count CSV → AnnMat), `plot_eig_dist()`, UMAP embedding, and dataset concatenation.

- **`bulk_functions.py`**: Bulk RNA-seq helpers. `run_deseq()` wraps pydeseq2; `run_go_enrichment()` wraps goatools. **Expects `metadata/` to be a sibling of `os.getcwd()`** — run notebooks from `scripts/` so that `os.path.dirname(os.getcwd())` resolves to the repo root.

- **`reader_functions.py`**: Parses legacy plate reader Excel files (BioTek format) into time-series OD/spectrum DataFrames. Entry points: `get_od_data()`, `get_spectrum_data()`.

- **`tecan_func.py`**: `tecan` class for Tecan-format plate reader files. Supports multi-channel data, background subtraction, spike removal, fluorescence normalization, and derivative calculation.

### `scripts/` — Analysis Pipelines

| Notebook | Purpose |
|---|---|
| `analysis_notebook.ipynb` | scRNA-seq processing: probe counts → AnnMat, cell calling, gene filtering, eigenvalue calculation |
| `bulk_analysis.ipynb` | Bulk RNA-seq: DESeq2 differential expression, GO enrichment |
| `scanpy_analysis.ipynb` | UMAP, Leiden clustering, marker gene analysis |
| `supplementary figures.ipynb` | Supplementary figure generation |
| `random_matrices.ipynb` | Synthetic Wishart matrices for GMP model validation |
| `permutation test.ipynb` | Permutation test for GO term enrichment significance |
| `SDS_experiment.ipynb` | Plate reader data for SDS sensitivity experiment |
| `model_fit.nb` | Wolfram Mathematica — GMP model fitting (not Python) |

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

### Key Metadata Files

Located in `metadata/`: `go-basic.obo` (GO ontology), `ecocyc.gaf` (E. coli GO annotations), `genomic.gtf` (gene ID ↔ name mapping used by `get_ID_conversion()`).
