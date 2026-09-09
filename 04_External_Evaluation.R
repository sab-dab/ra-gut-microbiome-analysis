# ============================================================
# SCRIPT 4: EXTERNAL VALIDATION METRICS AND FIGURES
# ============================================================



library(dplyr)
library(yardstick)
library(pROC)
library(ggplot2)

# ============================================================
# SCRIPT 4: RETUNED EXTERNAL VALIDATION METRICS AND FIGURES
# ============================================================
project_dir <- getwd() 
output_path <- file.path(
  project_dir,
  "external_validation_results_RETUNED"
)

metadata_file <- file.path(
  project_dir,
  "RA_validation",
  "SraRunTable.csv"
)

if (!file.exists(metadata_file)) {
  stop(
    "Missing RA_validation/SraRunTable.csv"
  )
}

# ============================================================
# 1) Packages
# ============================================================

library(dplyr)
library(yardstick)
library(pROC)
library(ggplot2)

# ============================================================
# 2) Load RETUNED external predictions
# ============================================================

external_results <- read.csv(
  file.path(
    output_path,
    "External_Genus_Predictions_RETUNED.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Prediction rows:",
  nrow(external_results),
  "\n"
)

# ============================================================
# 3) Load full external-study metadata
# ============================================================

metadata_all <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================================
# 4) Derive verified RA / HC labels
# ============================================================

metadata_selected <- metadata_all %>%
  dplyr::transmute(
    Sample = Run,
    Sample_Name = `Sample Name`,
    Group = dplyr::case_when(
      grepl("^RA_", `Sample Name`) ~ "RA",
      grepl("^GUT_", `Sample Name`) ~ "HC",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(
    Sample %in% external_results$Sample
  )

cat(
  "Predicted samples:",
  nrow(external_results),
  "\n"
)

cat(
  "Matched metadata samples:",
  nrow(metadata_selected),
  "\n"
)

cat(
  "RA samples:",
  sum(metadata_selected$Group == "RA"),
  "\n"
)

cat(
  "HC samples:",
  sum(metadata_selected$Group == "HC"),
  "\n"
)

cat(
  "Missing labels:",
  sum(is.na(metadata_selected$Group)),
  "\n"
)

# ============================================================
# 5) Join predictions with true labels
# ============================================================

evaluation_df <- external_results %>%
  dplyr::left_join(
    metadata_selected,
    by = "Sample"
  )

if (any(is.na(evaluation_df$Group))) {
  stop(
    "One or more processed samples do not have a verified RA/HC label."
  )
}

evaluation_df$Group <- factor(
  evaluation_df$Group,
  levels = c("HC", "RA")
)

evaluation_df$.pred_class <- factor(
  evaluation_df$.pred_class,
  levels = c("HC", "RA")
)

write.csv(
  evaluation_df,
  file.path(
    output_path,
    "External_Predictions_With_True_Labels_RETUNED.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6) Point-estimate classification metrics
# ============================================================

external_accuracy <- yardstick::accuracy(
  evaluation_df,
  truth = Group,
  estimate = .pred_class
)

external_sensitivity <- yardstick::sens(
  evaluation_df,
  truth = Group,
  estimate = .pred_class,
  event_level = "second"
)

external_specificity <- yardstick::spec(
  evaluation_df,
  truth = Group,
  estimate = .pred_class,
  event_level = "second"
)

external_auc <- yardstick::roc_auc(
  evaluation_df,
  truth = Group,
  .pred_RA,
  event_level = "second"
)

cat(
  "\nPoint estimates:\n"
)

print(external_accuracy)
print(external_sensitivity)
print(external_specificity)
print(external_auc)

# ============================================================
# 7) Confusion matrix
# ============================================================

external_confusion <- yardstick::conf_mat(
  evaluation_df,
  truth = Group,
  estimate = .pred_class
)

cat(
  "\nConfusion matrix:\n"
)

print(
  external_confusion
)

write.csv(
  as.data.frame(
    external_confusion$table
  ),
  file.path(
    output_path,
    "External_Confusion_Matrix_RETUNED.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8) ROC-AUC and 95% CI
# ============================================================

roc_obj <- pROC::roc(
  response = evaluation_df$Group,
  predictor = evaluation_df$.pred_RA,
  levels = c("HC", "RA"),
  direction = "<",
  quiet = TRUE
)

auc_value <- as.numeric(
  pROC::auc(
    roc_obj
  )
)

auc_ci <- pROC::ci.auc(
  roc_obj
)

cat(
  "\nROC-AUC:",
  round(
    auc_value,
    4
  ),
  "\n"
)

cat(
  "ROC-AUC 95% CI:",
  round(
    as.numeric(
      auc_ci[1]
    ),
    4
  ),
  "-",
  round(
    as.numeric(
      auc_ci[3]
    ),
    4
  ),
  "\n"
)

# ============================================================
# 9) Exact binomial 95% CIs
#    Accuracy, sensitivity, specificity
# ============================================================

n_correct <- sum(
  evaluation_df$Group ==
    evaluation_df$.pred_class
)

true_positive <- sum(
  evaluation_df$Group == "RA" &
    evaluation_df$.pred_class == "RA"
)

n_RA <- sum(
  evaluation_df$Group == "RA"
)

true_negative <- sum(
  evaluation_df$Group == "HC" &
    evaluation_df$.pred_class == "HC"
)

n_HC <- sum(
  evaluation_df$Group == "HC"
)

accuracy_ci <- binom.test(
  n_correct,
  nrow(evaluation_df),
  conf.level = 0.95
)$conf.int

sensitivity_ci <- binom.test(
  true_positive,
  n_RA,
  conf.level = 0.95
)$conf.int

specificity_ci <- binom.test(
  true_negative,
  n_HC,
  conf.level = 0.95
)$conf.int

cat(
  "\nExternal accuracy:",
  round(
    n_correct / nrow(evaluation_df),
    3
  ),
  "95% CI",
  round(
    accuracy_ci[1],
    3
  ),
  "-",
  round(
    accuracy_ci[2],
    3
  ),
  "\n"
)

cat(
  "External sensitivity:",
  round(
    true_positive / n_RA,
    3
  ),
  "95% CI",
  round(
    sensitivity_ci[1],
    3
  ),
  "-",
  round(
    sensitivity_ci[2],
    3
  ),
  "\n"
)

cat(
  "External specificity:",
  round(
    true_negative / n_HC,
    3
  ),
  "95% CI",
  round(
    specificity_ci[1],
    3
  ),
  "-",
  round(
    specificity_ci[2],
    3
  ),
  "\n"
)

# ============================================================
# 10) Final summary table
# ============================================================

external_metrics_summary <- data.frame(
  Samples = nrow(evaluation_df),
  RA = n_RA,
  HC = n_HC,
  
  Accuracy = n_correct / nrow(evaluation_df),
  Accuracy_CI_Lower = accuracy_ci[1],
  Accuracy_CI_Upper = accuracy_ci[2],
  
  Sensitivity = true_positive / n_RA,
  Sensitivity_CI_Lower = sensitivity_ci[1],
  Sensitivity_CI_Upper = sensitivity_ci[2],
  
  Specificity = true_negative / n_HC,
  Specificity_CI_Lower = specificity_ci[1],
  Specificity_CI_Upper = specificity_ci[2],
  
  ROC_AUC = auc_value,
  ROC_AUC_CI_Lower = as.numeric(
    auc_ci[1]
  ),
  ROC_AUC_CI_Upper = as.numeric(
    auc_ci[3]
  )
)

cat(
  "\n============================================\n"
)

cat(
  "RETUNED EXTERNAL VALIDATION METRICS\n"
)

cat(
  "============================================\n"
)

print(
  external_metrics_summary
)

write.csv(
  external_metrics_summary,
  file.path(
    output_path,
    "External_Validation_Metrics_RETUNED_With_95CI.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 11) ROC curve
# ============================================================

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = TPR
  )
) +
  geom_line(
    linewidth = 1.2
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.62,
    y = 0.10,
    label = paste0(
      "AUC = ",
      round(
        auc_value,
        3
      )
    ),
    size = 5
  ) +
  labs(
    title = "External Validation ROC Curve - Retuned Model",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_classic(
    base_size = 14
  )

print(
  p_roc
)

ggsave(
  filename = file.path(
    output_path,
    "Figure9_ROC_RETUNED.pdf"
  ),
  plot = p_roc,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(
    output_path,
    "Figure9_ROC_RETUNED.png"
  ),
  plot = p_roc,
  width = 6,
  height = 5,
  dpi = 600
)

# ============================================================
# 12) Confusion-matrix figure
# ============================================================

cm <- as.data.frame(
  external_confusion$table
)

p_cm <- ggplot(
  cm,
  aes(
    x = Truth,
    y = Prediction,
    fill = Freq
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = Freq
    ),
    size = 7
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    title = "External Validation Confusion Matrix - Retuned Model",
    x = "Observed group",
    y = "Predicted group",
    fill = "Samples"
  )

print(
  p_cm
)

ggsave(
  filename = file.path(
    output_path,
    "Figure9_ConfusionMatrix_RETUNED.pdf"
  ),
  plot = p_cm,
  width = 5,
  height = 5
)

ggsave(
  filename = file.path(
    output_path,
    "Figure9_ConfusionMatrix_RETUNED.png"
  ),
  plot = p_cm,
  width = 5,
  height = 5,
  dpi = 600
)

# ============================================================
# 13) Prediction-probability figure
# ============================================================

p_prob <- ggplot(
  evaluation_df,
  aes(
    x = Group,
    y = .pred_RA,
    fill = Group
  )
) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 2
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    title = "Predicted RA Probability - Retuned Model",
    x = "True Group",
    y = "Predicted Probability (RA)"
  )

print(
  p_prob
)

ggsave(
  filename = file.path(
    output_path,
    "Figure9_Probability_RETUNED.pdf"
  ),
  plot = p_prob,
  width = 5,
  height = 5
)

ggsave(
  filename = file.path(
    output_path,
    "Figure9_Probability_RETUNED.png"
  ),
  plot = p_prob,
  width = 5,
  height = 5,
  dpi = 600
)

# ============================================================
# 14) Sample-flow summary
# ============================================================

external_sample_flow <- data.frame(
  Stage = c(
    "FASTQ pairs downloaded",
    "Samples retained after sequence processing",
    "RA samples retained",
    "Healthy-control samples retained",
    "Discovery genera",
    "External genera",
    "Shared genera",
    "Named shared genera"
  ),
  Value = c(
    40,
    39,
    20,
    19,
    447,
    167,
    114,
    113
  )
)

write.csv(
  external_sample_flow,
  file.path(
    output_path,
    "Table_External_Validation_Sample_Flow_RETUNED.csv"
  ),
  row.names = FALSE
)

print(
  external_sample_flow
)

cat(
  "\n============================================\n"
)

cat(
  "RETUNED EXTERNAL EVALUATION COMPLETE\n"
)

cat(
  "Results saved in:\n",
  output_path,
  "\n"
)

cat(
  "============================================\n"
)
