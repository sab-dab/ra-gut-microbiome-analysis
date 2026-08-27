

library(dplyr)
library(tidymodels)
library(xgboost) 

external_input_dir <- file.path(
  "external_validation_results"
)

retuned_model_dir <- file.path(
  "genus_model_rebuild"
)

retuned_output_dir <- file.path(
  "external_validation_results_RETUNED"
)

if (!dir.exists(retuned_output_dir)) {
  dir.create(
    retuned_output_dir,
    recursive = TRUE
  )
}

# =============================================================================
# 1) Packages
# =============================================================================

library(dplyr)
library(tidymodels)
library(xgboost)

tidymodels::tidymodels_prefer()

# =============================================================================
# 2) Load the SAME external genus counts used previously
# =============================================================================

external_genus_counts <- readRDS(
  file.path(
    external_input_dir,
    "external_genus_counts_combined.rds"
  )
)

cat(
  "External genus matrix dimensions:",
  paste(dim(external_genus_counts), collapse = " x "),
  "\n"
)

# =============================================================================
# 3) Load rebuilt discovery genus dataset
# =============================================================================

discovery_training <- readRDS(
  file.path(
    retuned_model_dir,
    "discovery_genus_training_data_rebuilt.rds"
  )
)

discovery_features <- base::setdiff(
  colnames(discovery_training),
  c("Sample", "Group")
)

cat(
  "Discovery samples:",
  nrow(discovery_training),
  "\n"
)

cat(
  "Discovery predictors:",
  length(discovery_features),
  "\n"
)

# =============================================================================
# 4) Normalize genus names exactly as in original external prediction
# =============================================================================

normalize_genus_name <- function(x) {
  
  x <- tolower(x)
  
  x <- gsub(
    "^g__",
    "",
    x
  )
  
  x <- gsub(
    "\\.",
    "_",
    x
  )
  
  x <- gsub(
    "[^a-z0-9_]",
    "",
    x
  )
  
  trimws(x)
}

disc_clean <- normalize_genus_name(
  discovery_features
)

ext_clean <- normalize_genus_name(
  colnames(external_genus_counts)
)

colnames(external_genus_counts) <- ext_clean

names(disc_clean) <- discovery_features

# =============================================================================
# 5) Construct external prediction matrix with all discovery predictors
# =============================================================================

external_prediction <- matrix(
  0,
  nrow = nrow(external_genus_counts),
  ncol = length(discovery_features)
)

colnames(external_prediction) <- discovery_features

rownames(external_prediction) <- rownames(
  external_genus_counts
)

# =============================================================================
# 6) Populate shared genera
# =============================================================================

shared <- base::intersect(
  disc_clean,
  ext_clean
)

cat(
  "Shared genus predictors:",
  length(shared),
  "\n"
)

for (g in shared) {
  
  discovery_name <- names(
    disc_clean
  )[
    disc_clean == g
  ]
  
  external_prediction[
    ,
    discovery_name
  ] <- external_genus_counts[
    ,
    g
  ]
}

# =============================================================================
# 7) Convert to prediction data frame
# =============================================================================

external_prediction <- as.data.frame(
  external_prediction,
  check.names = FALSE
)

external_prediction$Sample <- rownames(
  external_prediction
)

external_prediction <- external_prediction[
  ,
  c(
    "Sample",
    discovery_features
  )
]

stopifnot(
  identical(
    colnames(external_prediction)[-1],
    discovery_features
  )
)

# =============================================================================
# 8) Load independently RETUNED genus XGBoost workflow
# =============================================================================

workflow_retuned <- readRDS(
  file.path(
    retuned_model_dir,
    "finalized_genus_xgboost_workflow_retuned.rds"
  )
)

cat(
  "\nRetuned model specification:\n"
)

print(
  workflows::extract_spec_parsnip(
    workflow_retuned
  )
)

# =============================================================================
# 9) Fit retuned workflow to COMPLETE discovery cohort
# =============================================================================

set.seed(123)

fitted_workflow_retuned <- fit(
  workflow_retuned,
  data = discovery_training
)

saveRDS(
  fitted_workflow_retuned,
  file.path(
    retuned_output_dir,
    "fitted_genus_xgboost_workflow_RETUNED.rds"
  )
)

# =============================================================================
# 10) Predict external cohort
# =============================================================================

pred_prob <- predict(
  fitted_workflow_retuned,
  new_data = external_prediction,
  type = "prob"
)

pred_class <- predict(
  fitted_workflow_retuned,
  new_data = external_prediction,
  type = "class"
)

external_results_retuned <- cbind(
  external_prediction["Sample"],
  pred_prob,
  pred_class
)

head(
  external_results_retuned
)

write.csv(
  external_results_retuned,
  file.path(
    retuned_output_dir,
    "External_Genus_Predictions_RETUNED.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 11) Save external feature-overlap information
# =============================================================================

feature_overlap_summary <- data.frame(
  Metric = c(
    "Discovery predictors",
    "Shared predictors",
    "Missing discovery predictors"
  ),
  Number = c(
    length(discovery_features),
    length(shared),
    length(discovery_features) - length(shared)
  )
)

write.csv(
  feature_overlap_summary,
  file.path(
    retuned_output_dir,
    "External_Feature_Overlap_RETUNED.csv"
  ),
  row.names = FALSE
)

cat(
  "\n============================================\n"
)

cat(
  "RETUNED EXTERNAL PREDICTION COMPLETE\n"
)

cat(
  "Predictions saved to:\n",
  retuned_output_dir,
  "\n"
)

cat(
  "============================================\n"
)
# =============================================================================
# FIGURE 9: EXTERNAL VALIDATION ROC CURVE
# =============================================================================

library(ggplot2)
library(pROC)

external_dir <- paste0(
  "RA_validation",
  "external_validation_results_RETUNED"
)

evaluation_df <- read.csv(
  file.path(
    external_dir,
    "External_Predictions_With_True_Labels_RETUNED.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

evaluation_df$Group <- factor(
  evaluation_df$Group,
  levels = c(
    "HC",
    "RA"
  )
)

roc_obj <- pROC::roc(
  response = evaluation_df$Group,
  predictor = evaluation_df$.pred_RA,
  levels = c(
    "HC",
    "RA"
  ),
  direction = "<",
  quiet = TRUE
)

roc_auc <- as.numeric(
  pROC::auc(
    roc_obj
  )
)

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p9 <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = TPR
  )
) +
  geom_line(
    linewidth = 1.1
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 0.63,
    y = 0.10,
    label = paste0(
      "AUC = ",
      sprintf(
        "%.3f",
        roc_auc
      )
    ),
    size = 5
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    title = "External Validation ROC Curve",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(
        b = 12
      )
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    plot.margin = margin(
      t = 20,
      r = 25,
      b = 20,
      l = 25
    )
  )

print(p9)

ggsave(
  file.path(
    external_dir,
    "Figure9_External_ROC_Final.png"
  ),
  p9,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    external_dir,
    "Figure9_External_ROC_Final.pdf"
  ),
  p9,
  width = 7.5,
  height = 6.5,
  units = "in",
  bg = "white"
)


file.exists(
  paste0(
    "RA_validation",
    "external_validation_results_RETUNED/",
    "External_Predictions_With_True_Labels_RETUNED.csv"
  )
)

# ============================================================
# FINAL EXTERNAL VALIDATION METRICS + 95% CIs
# ============================================================

library(pROC)

external_dir <- paste0(
  "RA_validation",
  "external_validation_results_RETUNED"
)

evaluation_df <- read.csv(
  file.path(
    external_dir,
    "External_Predictions_With_True_Labels_RETUNED.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

evaluation_df$Group <- factor(
  evaluation_df$Group,
  levels = c("HC", "RA")
)

evaluation_df$.pred_class <- factor(
  evaluation_df$.pred_class,
  levels = c("HC", "RA")
)

# ------------------------------------------------------------
# Confusion matrix counts
# ------------------------------------------------------------

TP <- sum(
  evaluation_df$Group == "RA" &
    evaluation_df$.pred_class == "RA"
)

TN <- sum(
  evaluation_df$Group == "HC" &
    evaluation_df$.pred_class == "HC"
)

FP <- sum(
  evaluation_df$Group == "HC" &
    evaluation_df$.pred_class == "RA"
)

FN <- sum(
  evaluation_df$Group == "RA" &
    evaluation_df$.pred_class == "HC"
)

N <- nrow(evaluation_df)

# ------------------------------------------------------------
# Point estimates
# ------------------------------------------------------------

accuracy <- (TP + TN) / N

sensitivity <- TP / (TP + FN)

specificity <- TN / (TN + FP)

# ------------------------------------------------------------
# Exact binomial 95% confidence intervals
# ------------------------------------------------------------

accuracy_ci <- binom.test(
  TP + TN,
  N
)$conf.int

sensitivity_ci <- binom.test(
  TP,
  TP + FN
)$conf.int

specificity_ci <- binom.test(
  TN,
  TN + FP
)$conf.int

# ------------------------------------------------------------
# ROC-AUC + 95% CI
# ------------------------------------------------------------

roc_obj <- pROC::roc(
  response = evaluation_df$Group,
  predictor = evaluation_df$.pred_RA,
  levels = c("HC", "RA"),
  direction = "<",
  quiet = TRUE
)

auc_value <- as.numeric(
  pROC::auc(roc_obj)
)

auc_ci <- pROC::ci.auc(
  roc_obj,
  conf.level = 0.95,
  method = "delong"
)

# ------------------------------------------------------------
# Final table
# ------------------------------------------------------------

final_metrics <- data.frame(
  Metric = c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "ROC-AUC"
  ),
  
  Estimate = c(
    accuracy,
    sensitivity,
    specificity,
    auc_value
  ),
  
  CI_Lower = c(
    accuracy_ci[1],
    sensitivity_ci[1],
    specificity_ci[1],
    auc_ci[1]
  ),
  
  CI_Upper = c(
    accuracy_ci[2],
    sensitivity_ci[2],
    specificity_ci[2],
    auc_ci[3]
  )
)

final_metrics[, 2:4] <- round(
  final_metrics[, 2:4],
  3
)

print(final_metrics)

cat("\nConfusion matrix counts:\n")

cat(
  "TP =", TP,
  " TN =", TN,
  " FP =", FP,
  " FN =", FN,
  "\n"
)

