# Supplementary Tables README

## Manuscript

**Title:** *Gut Microbiome-Based Prediction of Rheumatoid Arthritis Using Machine Learning: A Benchmarking Study Reveals Limited Cross-Cohort Generalizability*

**Authors:** Sabira Dabeer and Hatim Palitanawala

This folder contains the supplementary tables supporting the analyses reported in the revised manuscript. Files are provided in CSV format so that the complete results can be reviewed, searched, and reused without altering the original values.

## Contents

### Supplementary Table S1 - Primary cohort sequencing depth by sample
Provides sample-level sequencing-depth information for the primary cohort, including rheumatoid arthritis (RA) and healthy control (HC) samples.

### Supplementary Table S2 - Alpha-diversity results
Provides descriptive and inferential statistics for observed ASV richness and Shannon diversity in the RA and HC groups, including group summaries, Wilcoxon rank-sum test results, effect sizes, and confidence intervals where applicable.

### Supplementary Table S3 - Abundant and prevalent genera by community category
Summarizes the most abundant and most prevalent genera classified as shared between RA and HC, RA-specific, or HC-specific.

### Supplementary Table S4 - ANCOM-BC2 genus-level differential-abundance results
Provides the complete genus-level ANCOM-BC2 results for the primary cohort. Positive bias-corrected coefficients indicate higher abundance in RA, whereas negative coefficients indicate higher abundance in HC. The table includes robust differential-abundance and pseudocount-sensitivity results.

### Supplementary Table S5 - DESeq2 ASV-level differential-abundance results
Provides the complementary ASV-level DESeq2 results, including ASV identifiers, taxonomic annotations, log2 fold changes, standard errors, test statistics, raw P values, adjusted P values, and related abundance information where available.

### Supplementary Table S6 - Repeated internal-validation results
Summarizes model performance across 30 independent stratified 80/20 train-test splits for LASSO logistic regression, random forest, SVM-RBF, and XGBoost. Reported metrics include accuracy, sensitivity, specificity, and ROC-AUC, together with measures of variability and confidence intervals.

### Supplementary Table S7 - External-cohort differential-abundance results
Provides the exploratory genus-level differential-abundance results for the external Shanghai cohort after DADA2 processing, genus-level aggregation, prevalence filtering, and multiple-testing correction.

### Supplementary Table S8 - Nested ANCOM-BC2 feature-selection stability
Summarizes genus-level feature-selection stability across 30 repeated stratified train-test splits. ANCOM-BC2 feature selection was performed using the training partition only in each repetition. The table reports how frequently individual genera were selected across repetitions and supports the comparison between the full 447-genus XGBoost model and models restricted to training-only ANCOM-BC2-selected genera. Hyperparameters for the reduced-feature models were independently optimized within the corresponding training partitions.

### Supplementary Table S9 - Software and package versions
Lists the R version and versions of the principal packages used for microbiome analysis, statistical testing, machine learning, visualization, sequence processing, and reproducibility.

### Supplementary Table S10 - Transformation-sensitivity paired comparisons
Provides paired comparisons of genus-level XGBoost ROC-AUC across 30 identical stratified train-test splits using log(x + 1), relative-abundance, and centered log-ratio (CLR) feature representations. Comparisons were performed on matched splits using two-sided Wilcoxon signed-rank tests, with Holm adjustment for multiple comparisons.

### Supplementary Table S11 - Permutation-importance stability across repeated train-test splits
Provides the complete genus-level permutation-importance stability summary across 30 repeated stratified 80/20 train-test splits using the independently tuned, finalized genus-level XGBoost workflow. The table reports mean and median permutation importance, variability across repetitions, mean rank, and the frequency and percentage with which each genus appeared among the top 10 and top 20 predictors. A stable predictor was operationally defined as a genus appearing among the top 20 predictors in at least 50% of repetitions.

## Abbreviations

- **ASV:** Amplicon sequence variant
- **RA:** Rheumatoid arthritis
- **HC:** Healthy control
- **ANCOM-BC2:** Analysis of Compositions of Microbiomes with Bias Correction 2
- **LASSO:** Least absolute shrinkage and selection operator
- **SVM-RBF:** Support vector machine with a radial basis function kernel
- **ROC-AUC:** Area under the receiver operating characteristic curve
- **CLR:** Centered log-ratio
- **FDR:** False discovery rate

## Notes on interpretation

The supplementary tables should be interpreted together with the Methods, Results, figure legends, and limitations reported in the manuscript. Differential-abundance findings from the external cohort are exploratory because of the limited external sample size.

Repeated-validation, transformation-sensitivity, permutation-importance, and nested feature-selection analyses used stratified train-test partitions. Where applicable, preprocessing and feature-selection procedures were estimated using training data and applied unchanged to the corresponding held-out test data. Reduced-feature XGBoost models in the nested ANCOM-BC2 analysis were independently tuned within their corresponding training partitions.

The primary ASV-level machine-learning benchmark and the secondary genus-level analyses represent distinct analytical stages and should not be interpreted as directly interchangeable model-performance estimates.

## File format

All tabular files are supplied as comma-separated value (`.csv`) files. Missing or unavailable values may be represented as blank cells or `NA`, depending on the source analysis. Scientific notation is used where appropriate for small P values.

## Data and code availability

The primary dataset is publicly available through Figshare, and the external paired-end sequencing data are publicly available through the European Nucleotide Archive. Analysis code is available in the GitHub repository `sab-dab/ra-gut-microbiome-analysis`.