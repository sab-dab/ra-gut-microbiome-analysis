# =============================================================================
# 0) Paths
# =============================================================================
project_dir <- getwd()
output_dir <- file.path(
  project_dir,
  "genus_model_rebuild"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# =============================================================================
# 1) Packages
# =============================================================================

library(dplyr)
library(tibble)
library(stringr)
library(tidymodels)
library(xgboost)
library(pROC)

tidymodels::tidymodels_prefer()

set.seed(123)

# =============================================================================
# 2) Load original primary-cohort inputs
# =============================================================================

asv_file <- file.path(
  project_dir,
  "1.ASV.profile.rds"
)

taxonomy_file <- file.path(
  project_dir,
  "1.taxonomy.info.rds"
)

# Check that required input files are available
if (!file.exists(asv_file)) {
  stop(
    "Missing input file: 1.ASV.profile.rds. ",
    "Please place this file in the repository root directory."
  )
}

if (!file.exists(taxonomy_file)) {
  stop(
    "Missing input file: 1.taxonomy.info.rds. ",
    "Please place this file in the repository root directory."
  )
}

# Load input files
asv_profile <- readRDS(asv_file)
tax_info <- readRDS(taxonomy_file)

cat(
  "\nOriginal ASV profile dimensions:\n"
)

print(
  dim(asv_profile)
)

cat(
  "\nTaxonomy dimensions:\n"
)

print(
  dim(tax_info)
)
# =============================================================================
# 3) Convert ASV profile to samples x ASVs
# =============================================================================

# Existing data are ASVs x samples.
# Transpose so rows = samples and columns = ASVs.

asv_counts <- t(
  as.matrix(asv_profile)
)

cat(
  "\nSamples:",
  nrow(asv_counts),
  "\n"
)

cat(
  "ASVs:",
  ncol(asv_counts),
  "\n"
)

# =============================================================================
# 4) Create sample metadata
# =============================================================================

sample_names <- rownames(
  asv_counts
)

metadata <- data.frame(
  Sample = sample_names,
  Group = ifelse(
    grepl(
      "^HC",
      sample_names
    ),
    "HC",
    "RA"
  ),
  stringsAsFactors = FALSE
)

metadata$Group <- factor(
  metadata$Group,
  levels = c(
    "HC",
    "RA"
  )
)

cat(
  "\nGroup counts:\n"
)

print(
  table(metadata$Group)
)

# =============================================================================
# 5) Prepare taxonomy
# =============================================================================

tax_df <- as.data.frame(
  tax_info,
  stringsAsFactors = FALSE
) %>%
  tibble::rownames_to_column(
    "ASV"
  )

# Extract genus from the semicolon-delimited taxonomy field.
tax_df$Genus <- stringr::str_extract(
  tax_df$Taxon,
  "g__[^;]+"
)

# Preserve unclassified ASVs as Unknown.
tax_df$Genus <- ifelse(
  is.na(tax_df$Genus) |
    tax_df$Genus == "g__",
  "Unknown",
  stringr::str_remove(
    tax_df$Genus,
    "^g__"
  )
)

# Match taxonomy exactly to ASV-count columns.
tax_df <- tax_df[
  match(
    colnames(asv_counts),
    tax_df$ASV
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    tax_df$ASV,
    colnames(asv_counts)
  )
)

# =============================================================================
# 6) Normalize genus names for machine learning
# =============================================================================

normalize_genus_name <- function(x) {
  
  x <- as.character(x)
  
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
    "_",
    x
  )
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  x <- gsub(
    "^_|_$",
    "",
    x
  )
  
  x[x == ""] <- "unknown"
  
  x
}

tax_df$Genus_ML <- normalize_genus_name(
  tax_df$Genus
)

# =============================================================================
# 7) Aggregate raw ASV counts to genus level
# =============================================================================

genus_counts <- rowsum(
  t(asv_counts),
  group = tax_df$Genus_ML,
  reorder = FALSE
)

# Return to samples x genera.
genus_counts <- t(
  genus_counts
)

cat(
  "\nGenus-level matrix dimensions before any filtering:\n"
)

print(
  dim(genus_counts)
)

# =============================================================================
# IMPORTANT:
# Do NOT remove Unknown here.
#
# Your predictive workflow contained:
# 446 classified genera + 1 unclassified genus category.
# =============================================================================

unknown_columns <- colnames(
  genus_counts
)[
  grepl(
    "unknown|unclassified",
    colnames(genus_counts),
    ignore.case = TRUE
  )
]

cat(
  "\nUnclassified predictor(s):\n"
)

print(
  unknown_columns
)

# =============================================================================
# 8) Create genus-level ML dataset
# =============================================================================

genus_df <- as.data.frame(
  genus_counts,
  check.names = FALSE
)

genus_df$Sample <- rownames(
  genus_df
)

genus_df <- genus_df %>%
  dplyr::left_join(
    metadata,
    by = "Sample"
  ) %>%
  dplyr::select(
    Sample,
    Group,
    dplyr::everything()
  )

stopifnot(
  !anyNA(genus_df$Group)
)

predictor_names <- base::setdiff(
  colnames(genus_df),
  c(
    "Sample",
    "Group"
  )
)

cat(
  "\nDiscovery samples:",
  nrow(genus_df),
  "\n"
)

cat(
  "Discovery genus predictors:",
  length(predictor_names),
  "\n"
)

cat(
  "Classified predictors:",
  sum(
    !grepl(
      "unknown|unclassified",
      predictor_names,
      ignore.case = TRUE
    )
  ),
  "\n"
)

cat(
  "Unclassified predictors:",
  sum(
    grepl(
      "unknown|unclassified",
      predictor_names,
      ignore.case = TRUE
    )
  ),
  "\n"
)

# =============================================================================
# 9) Save rebuilt discovery files
# =============================================================================

saveRDS(
  genus_df,
  file.path(
    output_dir,
    "discovery_genus_training_data_rebuilt.rds"
  )
)

saveRDS(
  predictor_names,
  file.path(
    output_dir,
    "discovery_genus_training_features_rebuilt.rds"
  )
)

write.csv(
  genus_df,
  file.path(
    output_dir,
    "Discovery_Genus_Training_Data_Rebuilt.csv"
  ),
  row.names = FALSE
)



# =============================================================================
# STOP CONDITION
#
# Ideally the rebuilt matrix should contain exactly 447 predictors.
# If it does NOT, do not proceed blindly.
# =============================================================================

cat(
  "\nFinal rebuilt predictor count:",
  length(predictor_names),
  "\n"
)

if (length(predictor_names) != 447) {
  warning(
    paste0(
      "Rebuilt dataset contains ",
      length(predictor_names),
      " predictors instead of the expected 447. ",
      "Check taxonomy aggregation and input files before proceeding."
    )
  )
}

# =============================================================================
# 11) Stratified 80/20 split
# =============================================================================

set.seed(123)

genus_split <- initial_split(
  genus_df,
  prop = 0.80,
  strata = Group
)

train_genus <- training(
  genus_split
)

test_genus <- testing(
  genus_split
)

cat(
  "\nTraining samples:",
  nrow(train_genus),
  "\n"
)

cat(
  "Test samples:",
  nrow(test_genus),
  "\n"
)

cat(
  "\nTraining groups:\n"
)

print(
  table(train_genus$Group)
)

cat(
  "\nTest groups:\n"
)

print(
  table(test_genus$Group)
)

# =============================================================================
# 12) Genus-level preprocessing recipe
#
# This intentionally matches the saved genus workflow:
# step_zv + step_log(x + 1)
#
# No normalization is added here because the existing saved genus workflow
# did not contain step_normalize().
# =============================================================================

genus_recipe <- recipe(
  Group ~ .,
  data = train_genus
) %>%
  update_role(
    Sample,
    new_role = "id"
  ) %>%
  step_zv(
    all_predictors()
  ) %>%
  step_log(
    all_predictors(),
    offset = 1
  )

# =============================================================================
# 13) Five-fold CV inside training data only
# =============================================================================

set.seed(123)

genus_folds <- vfold_cv(
  train_genus,
  v = 5,
  strata = Group
)

# =============================================================================
# 14) Tunable genus-level XGBoost specification
#
# Keep trees = 1000 to remain comparable with the historical genus workflow.
# Independently tune all other major XGBoost hyperparameters.
# =============================================================================

genus_xgb_spec <- boost_tree(
  trees = 1000,
  learn_rate = tune(),
  mtry = tune(),
  tree_depth = tune(),
  min_n = tune(),
  loss_reduction = tune(),
  sample_size = tune()
) %>%
  set_engine(
    "xgboost",
    eval_metric = "auc",
    nthread = 1
  ) %>%
  set_mode(
    "classification"
  )

genus_xgb_workflow <- workflow() %>%
  add_recipe(
    genus_recipe
  ) %>%
  add_model(
    genus_xgb_spec
  )

# =============================================================================
# 15) Hyperparameter grid
# =============================================================================

set.seed(123)

genus_xgb_grid <- grid_latin_hypercube(
  
  learn_rate(
    range = c(
      -5,
      -1
    )
  ),
  
  mtry(
    range = c(
      5L,
      min(
        300L,
        length(predictor_names)
      )
    )
  ),
  
  tree_depth(
    range = c(
      2L,
      10L
    )
  ),
  
  min_n(
    range = c(
      2L,
      30L
    )
  ),
  
  loss_reduction(
    range = c(
      -5,
      1
    )
  ),
  
  sample_size = sample_prop(
    range = c(
      0.5,
      1.0
    )
  ),
  
  size = 30
)

print(
  genus_xgb_grid
)

# =============================================================================
# 16) Tune XGBoost
# =============================================================================

metric_function <- yardstick::metric_set(
  yardstick::roc_auc,
  yardstick::accuracy,
  yardstick::sens,
  yardstick::spec
)

set.seed(123)

genus_xgb_tuning <- tune::tune_grid(
  genus_xgb_workflow,
  resamples = genus_folds,
  grid = genus_xgb_grid,
  metrics = metric_function,
  control = tune::control_grid(
    save_pred = TRUE,
    verbose = TRUE
  )
)

# =============================================================================
# 17) Save complete tuning object
# =============================================================================

saveRDS(
  genus_xgb_tuning,
  file.path(
    output_dir,
    "genus_xgboost_tuning_results_retuned.rds"
  )
)

# Save tabular tuning results.
tuning_metrics <- tune::collect_metrics(
  genus_xgb_tuning
)

write.csv(
  tuning_metrics,
  file.path(
    output_dir,
    "Genus_XGBoost_Tuning_Metrics_Retuned.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 18) Select best parameters based on mean CV ROC-AUC
# =============================================================================

best_genus_xgb <- tune::select_best(
  genus_xgb_tuning,
  metric = "roc_auc"
)

cat(
  "\n============================================\n"
)

cat(
  "BEST RETUNED GENUS XGBOOST PARAMETERS\n"
)

cat(
  "============================================\n"
)

print(
  best_genus_xgb
)

write.csv(
  best_genus_xgb,
  file.path(
    output_dir,
    "Genus_XGBoost_Best_Parameters_Retuned.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 19) Finalize workflow
# =============================================================================

finalized_genus_xgb_retuned <- finalize_workflow(
  genus_xgb_workflow,
  best_genus_xgb
)

print(
  finalized_genus_xgb_retuned
)

saveRDS(
  finalized_genus_xgb_retuned,
  file.path(
    output_dir,
    "finalized_genus_xgboost_workflow_retuned.rds"
  )
)

# =============================================================================
# 20) Fit finalized workflow to training set
# =============================================================================

set.seed(123)

genus_fit_train <- fit(
  finalized_genus_xgb_retuned,
  data = train_genus
)

# =============================================================================
# 21) Evaluate independent held-out test set
# =============================================================================

test_predictions <- test_genus %>%
  dplyr::select(
    Sample,
    Group
  ) %>%
  bind_cols(
    predict(
      genus_fit_train,
      new_data = test_genus,
      type = "prob"
    ),
    predict(
      genus_fit_train,
      new_data = test_genus,
      type = "class"
    )
  )

test_predictions$Group <- factor(
  test_predictions$Group,
  levels = c(
    "HC",
    "RA"
  )
)

test_metrics <- data.frame(
  
  Accuracy = yardstick::accuracy(
    test_predictions,
    truth = Group,
    estimate = .pred_class
  )$.estimate,
  
  Sensitivity = yardstick::sens(
    test_predictions,
    truth = Group,
    estimate = .pred_class,
    event_level = "second"
  )$.estimate,
  
  Specificity = yardstick::spec(
    test_predictions,
    truth = Group,
    estimate = .pred_class,
    event_level = "second"
  )$.estimate,
  
  ROC_AUC = yardstick::roc_auc(
    test_predictions,
    truth = Group,
    .pred_RA,
    event_level = "second"
  )$.estimate
)

print(test_metrics)
cat(
  "\n============================================\n"
)

cat(
  "HELD-OUT GENUS MODEL PERFORMANCE\n"
)

cat(
  "============================================\n"
)

print(
  test_metrics
)

write.csv(
  test_metrics,
  file.path(
    output_dir,
    "Genus_XGBoost_Heldout_Performance_Retuned.csv"
  ),
  row.names = FALSE
)

write.csv(
  test_predictions,
  file.path(
    output_dir,
    "Genus_XGBoost_Heldout_Predictions_Retuned.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 22) ROC-AUC 95% CI
# =============================================================================

roc_obj <- pROC::roc(
  response = test_predictions$Group,
  predictor = test_predictions$.pred_RA,
  levels = c(
    "HC",
    "RA"
  ),
  direction = "<",
  quiet = TRUE
)

auc_ci <- pROC::ci.auc(
  roc_obj
)

cat(
  "\nHeld-out ROC-AUC:",
  round(
    as.numeric(
      pROC::auc(
        roc_obj
      )
    ),
    4
  ),
  "\n"
)

cat(
  "95% CI:",
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

# =============================================================================
# 23) Fit finalized retuned workflow to all 2,238 discovery samples
#
# This is the model that could later be used for external validation.
# =============================================================================

set.seed(123)

fitted_genus_xgb_full_retuned <- fit(
  finalized_genus_xgb_retuned,
  data = genus_df
)

saveRDS(
  fitted_genus_xgb_full_retuned,
  file.path(
    output_dir,
    "fitted_genus_xgboost_workflow_full_retuned.rds"
  )
)



# =============================================================================
# 24) Session information
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "Genus_Model_Retuning_SessionInfo.txt"
  )
)

cat(
  "\n============================================\n"
)

cat(
  "GENUS MODEL REBUILD AND RETUNING COMPLETE\n"
)

cat(
  "Results directory:\n",
  output_dir,
  "\n"
)

cat(
  "============================================\n"
)
