# dysregulated-persistence
This repository contains the code and data used in our paper:

**"Genome-wide Dysregulation in Antibiotic Tolerance and Persistence"**  
# 
## 📌 Overview

We present a computational framework to quantify the global correlation strength and level of dysregulation in bacterial cells under acute stress. We assess gene correlations by calculating the correlation spectrum of the experimental data and comparing it to our theoretical model - the Generalized Marchenko-Pastur (GMP) distribution.<br>
This repository includes:
- Code for reproducing all main figures
- Scripts for RNA-seq preprocessing and gene correlation analysis
- Derivation of the GMP distribution
- processed sequencing data and additional data that appears in the paper 
---

## 📦 Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/alongutf/dysregulated-persistence.git
git pull
conda env create -f environment.yml --name rnaseq_env
```
Installing source packages:
```python
import src.data_functions as df
import src.analysis_functions as af
import src.bulk_functions as bf
import src.reader_functions as rf
```
## 🛠️ Usage
All figures can be generated using the corresponding figureX.py file in scripts/figures.<br>
Genome files and scripts for mapping the scRNA-seq reads are located in scripts/RNA-seq mapping.<br>
Notebooks for different pipelines are available in scripts:
- analysis_notebook.ipynb: initial processing, including transformation of probe count data to cell-gene count matrices, cell calling, gene filtering and calculation of correlation eigenvalues.
- bulk_analysis.ipynb: analysis of the bulk RNA-seq data, including differential gene analysis, GO enrichment analysis and data processing.
- model_fit.nb: Wolfram Mathematica notebook for calculating the theoretical GMP distribution and fitting the empirical correlation eigenvalues to the model.
- permutation test: Perform random sampling test for significance of GO term enrichment difference.
- random_matrices.ipynb: generate random Wishart matrices with non-diagonal covariance - to compare with the GMP model.
- scanpy_analysis.ipynb: perform UMAP dimensionality reduction, marker gene analysis and further analysis of the single-cell data.
- supplementary figures.ipynb: generates additional figures that are not presented in the main text.
- SDS_experiment.ipynb: retrieves and analyzes the plate reader data from the SDS sensitivity experiment
---

