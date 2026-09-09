# Rheumatoid Arthritis Gut Microbiome Analysis

This repository contains the R scripts used to reproduce the analyses reported in the accompanying manuscript:

**Machine Learning Prediction of Rheumatoid Arthritis Using the Gut Microbiome Shows Limited Cross-Cohort Generalizability**

The study evaluates gut microbiome differences between individuals with rheumatoid arthritis (RA) and healthy controls (HC), benchmarks multiple machine-learning algorithms, assesses model stability across repeated train-test splits, evaluates genus-level feature transformations and permutation importance, performs leakage-free ANCOM-BC2 feature selection, and tests model transportability in an independent external cohort.

## Repository Contents

| Script | Description |
|---|---|
| `00_Main_Analysis_Pipeline.R` | Primary analysis workflow, including microbiome preprocessing, ecological analyses, ASV-level machine-learning benchmarking, repeated internal validation, and generation of primary results. |
| `01_DADA2_processing.R` | Reprocessing of the external paired-end 16S rRNA sequencing data using DADA2, including quality filtering, error learning, denoising, read merging, chimera removal, and taxonomic assignment. |
| `02_Genus_XGBoost_Model_Development.R` | Construction and independent tuning of the genus-level XGBoost workflow used for secondary genus-level analyses and external validation. |
| `03_External_Prediction.R` | Harmonization of discovery and external genus-level features and generation of predictions in the independent external cohort using the finalized genus-level XGBoost workflow. |
| `04_External_Evaluation.R` | Evaluation of external predictive performance, including accuracy, sensitivity, specificity, ROC-AUC, confidence intervals, confusion matrix, and ROC visualization. |
| `05_ANCOMBC2_Discovery.R` | Primary genus-level differential-abundance analysis using ANCOM-BC2. |
| `06_Transformation_Sensitivity.R` | Evaluation of genus-level XGBoost performance across log-transformed counts, relative abundance, and centered log-ratio (CLR) representations across repeated stratified train-test splits. |
| `07_Permutation_Importance_Stability.R` | Repeated held-out permutation-importance analysis used to quantify the stability of genus-level predictors across 30 stratified train-test splits. |
| `08_Nested_ANCOMBC2_Feature_Selection.R` | Leakage-free nested comparison of the full genus-level XGBoost model with models restricted to training-only ANCOM-BC2-selected genera. Reduced-feature models are independently tuned within the corresponding training partitions. |
| `External_Genus_Differential_Abundance.R` | Exploratory genus-level differential-abundance analysis of the external cohort. |

## Recommended Run Order

For the primary and external-validation workflows, the scripts can be run in the following general order:

1. `00_Main_Analysis_Pipeline.R`
2. `05_ANCOMBC2_Discovery.R`
3. `02_Genus_XGBoost_Model_Development.R`
4. `06_Transformation_Sensitivity.R`
5. `07_Permutation_Importance_Stability.R`
6. `08_Nested_ANCOMBC2_Feature_Selection.R`
7. `01_DADA2_processing.R`
8. `03_External_Prediction.R`
9. `04_External_Evaluation.R`
10. `External_Genus_Differential_Abundance.R`

The external-validation scripts require the external FASTQ files, metadata, and SILVA taxonomy reference described below.

## Analysis Overview

The primary machine-learning analysis was performed at the ASV level using 638 prevalence-filtered ASVs. Four classifiers were evaluated:

- LASSO logistic regression
- Random forest
- Support vector machine with a radial basis function kernel (SVM-RBF)
- XGBoost

Model performance was assessed using an independent held-out test set and 30 repeated stratified 80/20 train-test splits.

Secondary machine-learning analyses were performed at the genus level using 447 predictors, comprising 446 classified genera and one unclassified genus-level category. A separate genus-level XGBoost workflow was independently tuned and used for transformation-sensitivity analysis, permutation-importance stability, nested ANCOM-BC2 feature-selection analysis, and external validation.

For nested ANCOM-BC2 analysis, feature selection was restricted to the training partition within each repeated split. XGBoost hyperparameters for the reduced-feature models were independently optimized using cross-validation within the corresponding training data before evaluation on the untouched held-out test set.

## External Validation

The external Shanghai RA cohort was processed from publicly available paired-end FASTQ files using DADA2.

Because the primary cohort was obtained as a preprocessed ASV abundance and taxonomy dataset, the primary and external cohorts could not undergo identical upstream FASTQ-processing pipelines. Cross-cohort harmonization was therefore performed at the genus-feature level.

External predictive evaluation was performed without model retraining or hyperparameter optimization using the external cohort.

## Requirements

Analyses were performed using:

- **R version 4.6.1**
- `vegan`
- `ANCOMBC`
- `DESeq2`
- `dada2`
- `phyloseq`
- `tidymodels`
- `xgboost`
- `vip`
- `pROC`
- `dplyr`
- `tidyr`
- `tibble`
- `stringr`
- `ggplot2`

Users should install the required CRAN and Bioconductor packages before running the scripts. Exact package versions used in the final analysis are reported in Supplementary Table S9 of the manuscript.

## Reproducibility

Random seeds were fixed throughout the analyses where applicable.

Preprocessing operations used for machine-learning evaluation were estimated using training data and applied unchanged to the corresponding held-out test data. Training-only feature selection and nested model tuning were used where indicated to minimize information leakage.

The scripts are intended to be run from the repository root using project-relative file paths.

Required primary-cohort input files should be placed in the repository root:

- `1.ASV.profile.rds`
- `1.taxonomy.info.rds`

External-validation metadata should be placed at:

- `RA_validation/SraRunTable.csv`

External paired-end FASTQ files should be placed in:

- `RA_validation/fastq/`

The SILVA taxonomy reference used for DADA2 processing should be placed at:

- `reference/silva_nr99_v138.1_train_set.fa.gz`

Raw sequencing data and other source datasets are not distributed with this repository and should be obtained from the public repositories described in the manuscript.

## Supplementary Tables

Supplementary tables S1-S11 accompany the manuscript and include sequencing-depth summaries, diversity statistics, differential-abundance results, repeated internal-validation results, external differential-abundance results, nested ANCOM-BC2 feature-selection stability, software versions, transformation-sensitivity comparisons, and permutation-importance stability.

A separate `README_Supplementary_Tables.md` file describes the contents and interpretation of each supplementary table.

## Citation

If you use this code, please cite the associated manuscript once published. Citation details will be updated after publication.

## Contact

For questions regarding the analysis code or reproducibility of the workflow, please open a GitHub Issue.
