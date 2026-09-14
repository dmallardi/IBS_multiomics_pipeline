# IBS Multi-Omics Signature
This repository contains the R analysis pipeline used for multi-omics integration and microbial signature discovery in irritable bowel syndrome (IBS), followed by cross-disease evaluation in inflammatory bowel disease (IBD).
# Overview
The workflow integrates taxonomic, metatranscriptomic, and metabolomic data to identify multi-omics modules associated with IBS.
The main analytical steps include:
1. Clinical metadata harmonization
2. Taxonomic data preprocessing
3. Metatranscriptomic data preprocessing
4. Metabolomic data preprocessing
5. Multi-omics integration using MintTea
6. Characterization of multi-omics modules
7. PCA and ROC analysis
8. Permutation-based evaluation of the main IBS-associated module
9. LASSO feature selection
10. Construction of a reduced multi-omics IBS signature
11. Cross-disease evaluation in an independent IBD cohort

# Data
The IBS dataset includes taxonomic, metatranscriptomic, metabolomic, and clinical metadata.
The external IBD dataset is used to evaluate the cross-disease behavior of the IBS-derived multi-omics signature.
Raw datasets are not included in this repository. Instructions and references for obtaining the source datasets are provided in the `data` directory.Raw source datasets are not included in this repository. Derived data files prepared for the present analysis are provided in the data directory. Instructions and references for obtaining the original publicly available source datasets are also provided in the data directory.

## Requirements
The analysis was performed in R.
Main R packages:
- dplyr
- readr
- tidyr
- phyloseq
- MintTea
- ggplot2
- igraph
- ggraph
- tidygraph
- pheatmap
- pROC
- glmnet
- caret
- car
- reticulate

# Reproducibility
Random seeds are explicitly defined where required, including permutation and cross-validation analyses.
The main analysis can be reproduced using:
`IBS_multiomics_pipeline.R`

# Citation
If you use this code, please cite the associated manuscript:
[Citation to be added upon publication]
