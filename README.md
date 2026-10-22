# IBS Multi-Omics Signature

This repository contains the R analysis pipeline used for multi-omics integration and microbial signature discovery in irritable bowel syndrome (IBS), followed by cross-disease evaluation in inflammatory bowel disease (IBD).

## Overview

The workflow integrates taxonomic, metatranscriptomic, and metabolomic data to identify multi-omics modules associated with IBS and to develop a reduced IBS-associated signature.

The main analytical steps include:

1. Clinical metadata harmonization
2. Taxonomic data preprocessing
3. Metatranscriptomic data preprocessing
4. Metabolomic data preprocessing
5. Multi-omics integration using MintTea
6. Characterization of multi-omics modules
7. Principal component analysis (PCA) and receiver operating characteristic (ROC) analysis
8. Permutation-based evaluation of the main IBS-associated module
9. LASSO feature selection
10. Construction of a reduced multi-omics IBS signature
11. Cross-disease evaluation in an independent IBD cohort

## Repository Structure

- `Dataset/`: Contains all data files used in the study, including processed and unprocessed datasets, input files for MintTea integration and LASSO feature selection, and other files required for the analyses.
- `IBS_multiomics_pipeline.R`: Main R script implementing the analysis pipeline.

## Data Availability

The `Dataset/` directory contains the data files used in the analyses, including unprocessed (publicly available source file)  and processed datasets and input files required to reproduce the main analytical steps.

## Requirements

The analyses were performed in R.

The main R packages used in the pipeline include:

- `dplyr`
- `readr`
- `tidyr`
- `phyloseq`
- `MintTea`
- `ggplot2`
- `igraph`
- `ggraph`
- `tidygraph`
- `pheatmap`
- `pROC`
- `glmnet`
- `caret`
- `car`
- `reticulate`

The required packages should be installed before running the analysis scripts.

## Reproducibility

This repository provides the data files, analysis inputs, and corresponding scripts used for the analyses reported in the associated manuscript, facilitating transparency and reproducibility.

Random seeds are explicitly defined where required, including permutation and cross-validation analyses.

The main analysis can be reproduced using:

`IBS_multiomics_pipeline.R`

The input files are available in the `Dataset/` directory. Scripts should be run with the appropriate working directory and the required R packages installed.

## Citation

If you use this code or the associated datasets, please cite the corresponding manuscript:

Citation to be added upon publication.
