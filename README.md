# Rheumatoid Arthritis Gut Microbiome Analysis
This repository contains the R scripts used to reproduce the analyses presented in the accompanying manuscript investigating gut microbiome composition and machine learning–based prediction of rheumatoid arthritis.

## Repository Contents

| Script | Description |
|--------|-------------|
| `00_Main_Analysis_Pipeline.R` | Main workflow that coordinates the complete analysis pipeline. |
| `01_DADA2_processing.R` | Processing of external paired-end 16S rRNA sequencing data using DADA2. |
| `03_External_Prediction.R` | External prediction workflow. |
| `04_External_Evaluation.R` | Evaluation of model performance on the external cohort. |
| `05_ANCOMBC2_Discovery.R` | Differential abundance analysis using ANCOM-BC2. |
| `06_Transformation_Sensitivity.R` | Sensitivity analysis across abundance transformations. |
| `07_Permutation_Importance_Stability.R` | Permutation importance stability analysis. |
| `08_Nested_ANCOMBC2_Feature_Selection.R` | Nested feature selection and model evaluation. |
| `External_Genus_Differential_Abundance.R` | Differential abundance analysis of the external cohort. |

## Requirements

- R version 4.6.1
- Required R packages are loaded within the individual scripts.

## Citation

If you use this code, please cite the associated manuscript once published. Citation details will be updated after publication.

## Contact

For questions regarding this repository, please open a GitHub Issue. 
## Data Availability

The analyses were performed using publicly available 16S rRNA sequencing datasets.

- ** Primary Dataset:** Li et al., Scientific Data (2025). The dataset is available from the public repository referenced in the manuscript.
- **External validation cohort:** Sun et al. The raw sequencing data are available from the public repository referenced in the manuscript.

To comply with the original data providers' terms of use, this repository contains analysis code only. Users should download the datasets from the original repositories before running the scripts.
