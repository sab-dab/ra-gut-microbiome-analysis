model_dir <- file.path(
  project_dir,
  "genus_model_rebuild"
)

save_dir <- file.path(
  project_dir,
  "results_final_RETUNED"
)

if (!dir.exists(save_dir)) {
  dir.create(
    save_dir,
    recursive = TRUE
  )
}

# =============================================================================
# 1) Libraries
# =============================================================================

library(ANCOMBC)
library(dplyr)
library(tidyr)
library(tibble)
library(tidymodels)
library(xgboost)
library(ggplot2)

tidymodels::tidymodels_prefer()

# =============================================================================
# 2) Load rebuilt discovery genus-level raw-count data
# =============================================================================

discovery_df <- readRDS(
  file.path(
    model_dir,
    "discovery_genus_training_data_rebuilt.rds"
  )
)

discovery_df$Group <- factor(
  discovery_df$Group,
  levels = c("HC", "RA")
)

predictor_names <- base::setdiff(
  colnames(discovery_df),
  c("Sample", "Group")
)

cat(
  "Samples:",
  nrow(discovery_df),
  "\n"
)

cat(
  "Predictors:",
  length(predictor_names),
  "\n"
)

cat(
  "HC:",
  sum(discovery_df$Group == "HC"),
  "\n"
)

cat(
  "RA:",
  sum(discovery_df$Group == "RA"),
  "\n"
)

stopifnot(
  nrow(discovery_df) == 2238
)

stopifnot(
  length(predictor_names) == 447
)

# =============================================================================
# 3) Load independently retuned reference genus workflow
# =============================================================================

retuned_reference_workflow <- readRDS(
  file.path(
    model_dir,
    "finalized_genus_xgboost_workflow_retuned.rds"
  )
)

retuned_reference_spec <-
  workflows::extract_spec_parsnip(
    retuned_reference_workflow
  )

cat(
  "\nRETUNED reference XGBoost specification:\n"
)

print(
  retuned_reference_spec
)

# =============================================================================
# 4) Repetition settings
# =============================================================================

N_REPEATS <- 30

set.seed(123)

all_repeat_seeds <- sample(
  1000:999999,
  size = 30,
  replace = FALSE
)

repeat_seeds <- all_repeat_seeds[
  seq_len(N_REPEATS)
]

progress_file <- file.path(
  save_dir,
  "Nested_ANCOMBC2_RETUNED_Progress.csv"
)

selected_features_file <- file.path(
  save_dir,
  "Nested_ANCOMBC2_RETUNED_Selected_Features_Progress.csv"
)

tuning_parameters_file <- file.path(
  save_dir,
  "Nested_ANCOMBC2_RETUNED_Selected_Model_Tuning_Progress.csv"
)

# =============================================================================
# 5) Common preprocessing recipe
#
# Matches the independently retuned genus workflow:
# - zero-variance removal
# - log(x + 1)
#
# No normalization.
# =============================================================================

create_model_recipe <- function(
    training_data,
    selected_predictors) {
  
  reduced_training_data <- training_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(
        selected_predictors
      )
    )
  
  recipes::recipe(
    Group ~ .,
    data = reduced_training_data
  ) %>%
    recipes::update_role(
      Sample,
      new_role = "id"
    ) %>%
    recipes::step_zv(
      recipes::all_predictors()
    ) %>%
    recipes::step_log(
      recipes::all_predictors(),
      offset = 1
    )
}

# =============================================================================
# 6) ALL-FEATURE reference model
#
# Uses exactly the independently retuned genus XGBoost specification.
# =============================================================================

fit_and_evaluate_reference_xgb <- function(
    train_data,
    test_data,
    selected_predictors,
    model_name,
    seed_value) {
  
  train_selected <- train_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(
        selected_predictors
      )
    )
  
  test_selected <- test_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(
        selected_predictors
      )
    )
  
  model_recipe <- create_model_recipe(
    training_data = train_selected,
    selected_predictors =
      selected_predictors
  )
  
  model_workflow <-
    workflows::workflow() %>%
    workflows::add_recipe(
      model_recipe
    ) %>%
    workflows::add_model(
      retuned_reference_spec
    )
  
  set.seed(seed_value)
  
  fitted_model <- parsnip::fit(
    model_workflow,
    data = train_selected
  )
  
  predictions <- test_selected %>%
    dplyr::select(
      Sample,
      Group
    ) %>%
    dplyr::bind_cols(
      predict(
        fitted_model,
        new_data = test_selected,
        type = "prob"
      ),
      predict(
        fitted_model,
        new_data = test_selected,
        type = "class"
      )
    ) %>%
    dplyr::mutate(
      Group = factor(
        Group,
        levels = c(
          "HC",
          "RA"
        )
      )
    )
  
  data.frame(
    Seed = seed_value,
    Model = model_name,
    Number_of_Features =
      length(selected_predictors),
    
    Accuracy =
      yardstick::accuracy(
        predictions,
        truth = Group,
        estimate = .pred_class
      )$.estimate,
    
    Sensitivity =
      yardstick::sens(
        predictions,
        truth = Group,
        estimate = .pred_class,
        event_level = "second"
      )$.estimate,
    
    Specificity =
      yardstick::spec(
        predictions,
        truth = Group,
        estimate = .pred_class,
        event_level = "second"
      )$.estimate,
    
    ROC_AUC =
      yardstick::roc_auc(
        predictions,
        truth = Group,
        .pred_RA,
        event_level = "second"
      )$.estimate
  )
}

# =============================================================================
# 7) Independently tuned XGBoost for ANCOM-BC2-selected genera
# =============================================================================

fit_and_evaluate_tuned_xgb <- function(
    train_data,
    test_data,
    selected_predictors,
    model_name,
    seed_value) {
  
  train_selected <- train_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(
        selected_predictors
      )
    )
  
  test_selected <- test_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(
        selected_predictors
      )
    )
  
  model_recipe <- create_model_recipe(
    training_data = train_selected,
    selected_predictors =
      selected_predictors
  )
  
  # All major XGBoost hyperparameters are independently tuned.
  tuned_spec <- parsnip::boost_tree(
    mtry = tune::tune(),
    trees = tune::tune(),
    min_n = tune::tune(),
    tree_depth = tune::tune(),
    learn_rate = tune::tune(),
    loss_reduction = tune::tune(),
    sample_size = tune::tune()
  ) %>%
    parsnip::set_engine(
      "xgboost",
      eval_metric = "auc",
      nthread = 1
    ) %>%
    parsnip::set_mode(
      "classification"
    )
  
  tuned_workflow <-
    workflows::workflow() %>%
    workflows::add_recipe(
      model_recipe
    ) %>%
    workflows::add_model(
      tuned_spec
    )
  
  # Inner CV restricted entirely to the outer training partition.
  set.seed(seed_value)
  
  inner_folds <- rsample::vfold_cv(
    train_selected,
    v = 5,
    strata = Group
  )
  
  set.seed(seed_value)
  
  tuning_results <- tune::tune_grid(
    tuned_workflow,
    resamples = inner_folds,
    grid = 10,
    metrics =
      yardstick::metric_set(
        yardstick::roc_auc
      ),
    control =
      tune::control_grid(
        save_pred = FALSE,
        verbose = FALSE
      )
  )
  
  best_parameters <- tune::select_best(
    tuning_results,
    metric = "roc_auc"
  )
  
  final_workflow <-
    tune::finalize_workflow(
      tuned_workflow,
      best_parameters
    )
  
  set.seed(seed_value)
  
  fitted_model <- parsnip::fit(
    final_workflow,
    data = train_selected
  )
  
  predictions <- test_selected %>%
    dplyr::select(
      Sample,
      Group
    ) %>%
    dplyr::bind_cols(
      predict(
        fitted_model,
        new_data = test_selected,
        type = "prob"
      ),
      predict(
        fitted_model,
        new_data = test_selected,
        type = "class"
      )
    ) %>%
    dplyr::mutate(
      Group = factor(
        Group,
        levels = c(
          "HC",
          "RA"
        )
      )
    )
  
  metrics_result <- data.frame(
    Seed = seed_value,
    Model = model_name,
    Number_of_Features =
      length(selected_predictors),
    
    Accuracy =
      yardstick::accuracy(
        predictions,
        truth = Group,
        estimate = .pred_class
      )$.estimate,
    
    Sensitivity =
      yardstick::sens(
        predictions,
        truth = Group,
        estimate = .pred_class,
        event_level = "second"
      )$.estimate,
    
    Specificity =
      yardstick::spec(
        predictions,
        truth = Group,
        estimate = .pred_class,
        event_level = "second"
      )$.estimate,
    
    ROC_AUC =
      yardstick::roc_auc(
        predictions,
        truth = Group,
        .pred_RA,
        event_level = "second"
      )$.estimate
  )
  
  return(
    list(
      metrics = metrics_result,
      best_parameters =
        best_parameters
    )
  )
}

# =============================================================================
# 8) Training-only ANCOM-BC2 feature selection
# =============================================================================

select_features_ancombc2 <- function(
    train_data,
    seed_value,
    prevalence_threshold = 0.05) {
  
  train_predictors <- as.matrix(
    train_data[
      ,
      predictor_names,
      drop = FALSE
    ]
  )
  
  storage.mode(
    train_predictors
  ) <- "numeric"
  
  # Unknown retained for ML but excluded from DA feature selection.
  valid_predictors <- colnames(
    train_predictors
  )[
    tolower(
      colnames(
        train_predictors
      )
    ) != "unknown"
  ]
  
  train_predictors <-
    train_predictors[
      ,
      valid_predictors,
      drop = FALSE
    ]
  
  # Training-only prevalence filter.
  prevalence_counts <- colSums(
    train_predictors > 0,
    na.rm = TRUE
  )
  
  minimum_samples <- ceiling(
    prevalence_threshold *
      nrow(train_predictors)
  )
  
  retained_predictors <- names(
    prevalence_counts[
      prevalence_counts >=
        minimum_samples
    ]
  )
  
  train_predictors_filt <-
    train_predictors[
      ,
      retained_predictors,
      drop = FALSE
    ]
  
  cat(
    "  ANCOM-BC2 training genera:",
    ncol(train_predictors_filt),
    "\n"
  )
  
  # taxa x samples
  ancom_matrix <- t(
    train_predictors_filt
  )
  
  ancom_metadata <- data.frame(
    Sample = train_data$Sample,
    Group = factor(
      train_data$Group,
      levels = c(
        "HC",
        "RA"
      )
    ),
    row.names =
      train_data$Sample
  )
  
  colnames(ancom_matrix) <-
    train_data$Sample
  
  stopifnot(
    identical(
      colnames(ancom_matrix),
      rownames(ancom_metadata)
    )
  )
  
  set.seed(seed_value)
  
  ancom_fit <- ANCOMBC::ancombc2(
    data = ancom_matrix,
    meta_data = ancom_metadata,
    fix_formula = "Group",
    p_adj_method = "BH",
    prv_cut = 0,
    lib_cut = 0,
    group = "Group",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE
  )
  
  ancom_results <- ancom_fit$res
  
  selected_features <-
    ancom_results %>%
    dplyr::filter(
      diff_robust_GroupRA ==
        TRUE
    ) %>%
    dplyr::pull(
      taxon
    )
  
  selected_features <-
    base::intersect(
      selected_features,
      predictor_names
    )
  
  list(
    selected_features =
      selected_features,
    ancom_results =
      ancom_results,
    genera_tested =
      nrow(ancom_results)
  )
}

# =============================================================================
# 9) Start fresh RETUNED progress
# =============================================================================

if (file.exists(progress_file)) {
  
  comparison_results <- read.csv(
    progress_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  completed_seeds <-
    comparison_results %>%
    dplyr::count(
      Seed
    ) %>%
    dplyr::filter(
      n >= 2
    ) %>%
    dplyr::pull(
      Seed
    )
  
  cat(
    "Completed RETUNED repetitions found:",
    length(completed_seeds),
    "\n"
  )
  
} else {
  
  comparison_results <- data.frame()
  
  completed_seeds <-
    numeric(0)
  
  cat(
    "No completed RETUNED nested repetitions found.\n"
  )
}

if (
  file.exists(
    selected_features_file
  )
) {
  
  selected_feature_results <-
    read.csv(
      selected_features_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  
} else {
  
  selected_feature_results <-
    data.frame()
}

if (
  file.exists(
    tuning_parameters_file
  )
) {
  
  tuning_parameter_results <-
    read.csv(
      tuning_parameters_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  
} else {
  
  tuning_parameter_results <-
    data.frame()
}

remaining_seeds <- repeat_seeds[
  !repeat_seeds %in%
    completed_seeds
]

cat(
  "Repetitions remaining:",
  length(remaining_seeds),
  "\n"
)

# =============================================================================
# 10) Fully nested repeated comparison
# =============================================================================

for (
  current_seed in
  remaining_seeds
) {
  
  repetition_number <- match(
    current_seed,
    all_repeat_seeds
  )
  
  cat(
    "\n====================================================\n"
  )
  
  cat(
    "Starting RETUNED nested repetition",
    repetition_number,
    "of",
    N_REPEATS,
    "- seed:",
    current_seed,
    "\n"
  )
  
  # ---------------------------------------------------------------------------
  # Outer train-test split
  # ---------------------------------------------------------------------------
  
  set.seed(
    current_seed
  )
  
  split_rep <- rsample::initial_split(
    discovery_df,
    prop = 0.8,
    strata = Group
  )
  
  train_rep <- rsample::training(
    split_rep
  )
  
  test_rep <- rsample::testing(
    split_rep
  )
  
  # ---------------------------------------------------------------------------
  # Model 1: all 447 features
  # ---------------------------------------------------------------------------
  
  cat(
    "  Fitting RETUNED all-features model...\n"
  )
  
  all_features_result <-
    fit_and_evaluate_reference_xgb(
      train_data = train_rep,
      test_data = test_rep,
      selected_predictors =
        predictor_names,
      model_name =
        "All features",
      seed_value =
        current_seed
    )
  
  cat(
    "  All-features ROC-AUC:",
    round(
      all_features_result$ROC_AUC,
      4
    ),
    "\n"
  )
  
  # ---------------------------------------------------------------------------
  # Training-only ANCOM-BC2
  # ---------------------------------------------------------------------------
  
  cat(
    "  Running training-only ANCOM-BC2...\n"
  )
  
  selection_result <-
    select_features_ancombc2(
      train_data =
        train_rep,
      seed_value =
        current_seed,
      prevalence_threshold =
        0.05
    )
  
  selected_features <-
    selection_result$selected_features
  
  cat(
    "  Robust selected genera:",
    length(
      selected_features
    ),
    "\n"
  )
  
  # Save selected genera.
  if (
    length(selected_features) >
    0
  ) {
    
    selected_features_current <-
      data.frame(
        Seed =
          current_seed,
        Repetition =
          repetition_number,
        Feature =
          selected_features,
        stringsAsFactors =
          FALSE
      )
    
    selected_feature_results <-
      dplyr::bind_rows(
        selected_feature_results,
        selected_features_current
      ) %>%
      dplyr::distinct(
        Seed,
        Feature,
        .keep_all = TRUE
      )
    
    write.csv(
      selected_feature_results,
      selected_features_file,
      row.names = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Model 2: ANCOM-selected model with INDEPENDENT tuning
  # ---------------------------------------------------------------------------
  
  if (
    length(selected_features) >=
    2
  ) {
    
    cat(
      "  Independently tuning ANCOM-selected model...\n"
    )
    
    tuned_selected <-
      fit_and_evaluate_tuned_xgb(
        train_data =
          train_rep,
        test_data =
          test_rep,
        selected_predictors =
          selected_features,
        model_name =
          "Nested ANCOM-BC2 selected",
        seed_value =
          current_seed
      )
    
    selected_result <-
      tuned_selected$metrics
    
    best_parameters_current <-
      tuned_selected$best_parameters %>%
      dplyr::mutate(
        Seed =
          current_seed,
        Repetition =
          repetition_number,
        Number_of_Features =
          length(
            selected_features
          )
      )
    
    tuning_parameter_results <-
      dplyr::bind_rows(
        tuning_parameter_results,
        best_parameters_current
      )
    
    write.csv(
      tuning_parameter_results,
      tuning_parameters_file,
      row.names = FALSE
    )
    
    cat(
      "  Selected-feature ROC-AUC:",
      round(
        selected_result$ROC_AUC,
        4
      ),
      "\n"
    )
    
  } else {
    
    warning(
      paste0(
        "Fewer than two robust ANCOM-BC2 genera selected for seed ",
        current_seed
      )
    )
    
    selected_result <-
      data.frame(
        Seed =
          current_seed,
        Model =
          "Nested ANCOM-BC2 selected",
        Number_of_Features =
          length(
            selected_features
          ),
        Accuracy =
          NA_real_,
        Sensitivity =
          NA_real_,
        Specificity =
          NA_real_,
        ROC_AUC =
          NA_real_
      )
  }
  
  # ---------------------------------------------------------------------------
  # Save completed outer repetition
  # ---------------------------------------------------------------------------
  
  current_results <-
    dplyr::bind_rows(
      all_features_result,
      selected_result
    )
  
  comparison_results <-
    dplyr::bind_rows(
      comparison_results,
      current_results
    ) %>%
    dplyr::distinct(
      Seed,
      Model,
      .keep_all = TRUE
    )
  
  write.csv(
    comparison_results,
    progress_file,
    row.names = FALSE
  )
  
  cat(
    "  Repetition",
    repetition_number,
    "completed and saved.\n"
  )
  
  rm(
    split_rep,
    train_rep,
    test_rep,
    all_features_result,
    selection_result,
    selected_features,
    selected_result,
    current_results
  )
  
  invisible(
    gc()
  )
}

# =============================================================================
# 11) Final results
# =============================================================================

nested_results <-
  comparison_results %>%
  dplyr::filter(
    Seed %in%
      repeat_seeds
  ) %>%
  dplyr::arrange(
    match(
      Seed,
      repeat_seeds
    ),
    Model
  )

write.csv(
  nested_results,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_RETUNED_All_Features_vs_Selected_All_Runs.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCompleted seeds:",
  dplyr::n_distinct(
    nested_results$Seed
  ),
  "\n"
)

cat(
  "Total result rows:",
  nrow(nested_results),
  "\n"
)

stopifnot(
  dplyr::n_distinct(
    nested_results$Seed
  ) == 30
)

stopifnot(
  nrow(nested_results) == 60
)

# =============================================================================
# 12) Performance summary
# =============================================================================

nested_summary <-
  nested_results %>%
  tidyr::pivot_longer(
    cols = c(
      Accuracy,
      Sensitivity,
      Specificity,
      ROC_AUC
    ),
    names_to =
      "Metric",
    values_to =
      "Value"
  ) %>%
  dplyr::group_by(
    Model,
    Metric
  ) %>%
  dplyr::summarise(
    Repetitions =
      sum(
        !is.na(Value)
      ),
    
    Mean =
      mean(
        Value,
        na.rm = TRUE
      ),
    
    SD =
      stats::sd(
        Value,
        na.rm = TRUE
      ),
    
    Median =
      stats::median(
        Value,
        na.rm = TRUE
      ),
    
    Q1 =
      stats::quantile(
        Value,
        0.25,
        na.rm = TRUE
      ),
    
    Q3 =
      stats::quantile(
        Value,
        0.75,
        na.rm = TRUE
      ),
    
    Minimum =
      min(
        Value,
        na.rm = TRUE
      ),
    
    Maximum =
      max(
        Value,
        na.rm = TRUE
      ),
    
    CI_Lower =
      Mean -
      stats::qt(
        0.975,
        df =
          Repetitions - 1
      ) *
      SD /
      sqrt(
        Repetitions
      ),
    
    CI_Upper =
      Mean +
      stats::qt(
        0.975,
        df =
          Repetitions - 1
      ) *
      SD /
      sqrt(
        Repetitions
      ),
    
    .groups = "drop"
  )

print(
  nested_summary
)

write.csv(
  nested_summary,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_RETUNED_Summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\nROC-AUC summary:\n"
)

print(
  nested_summary %>%
    dplyr::filter(
      Metric ==
        "ROC_AUC"
    )
)

# =============================================================================
# 13) Paired ROC-AUC comparison
# =============================================================================

auc_wide <-
  nested_results %>%
  dplyr::select(
    Seed,
    Model,
    ROC_AUC
  ) %>%
  tidyr::pivot_wider(
    names_from =
      Model,
    values_from =
      ROC_AUC
  )

colnames(
  auc_wide
) <- make.names(
  colnames(
    auc_wide
  )
)

auc_wide <-
  auc_wide %>%
  dplyr::mutate(
    AUC_Difference_Selected_minus_All =
      Nested.ANCOM.BC2.selected -
      All.features
  )

write.csv(
  auc_wide,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_RETUNED_Paired_AUC_Differences.csv"
  ),
  row.names = FALSE
)

paired_auc_test <-
  stats::wilcox.test(
    auc_wide$Nested.ANCOM.BC2.selected,
    auc_wide$All.features,
    paired = TRUE,
    exact = FALSE,
    conf.int = TRUE
  )

cat(
  "\nPaired ROC-AUC test:\n"
)

print(
  paired_auc_test
)

paired_test_results <-
  data.frame(
    Test =
      "Paired Wilcoxon signed-rank test",
    
    Repetitions =
      sum(
        complete.cases(
          auc_wide$Nested.ANCOM.BC2.selected,
          auc_wide$All.features
        )
      ),
    
    Mean_AUC_Difference_Selected_minus_All =
      mean(
        auc_wide$AUC_Difference_Selected_minus_All,
        na.rm = TRUE
      ),
    
    Median_AUC_Difference_Selected_minus_All =
      stats::median(
        auc_wide$AUC_Difference_Selected_minus_All,
        na.rm = TRUE
      ),
    
    Wilcoxon_V =
      unname(
        paired_auc_test$statistic
      ),
    
    P_Value =
      paired_auc_test$p.value
  )

print(
  paired_test_results
)

write.csv(
  paired_test_results,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_RETUNED_Paired_AUC_Test.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 14) Feature-selection stability
# =============================================================================

feature_selection_stability <-
  selected_feature_results %>%
  dplyr::filter(
    Seed %in%
      repeat_seeds
  ) %>%
  dplyr::count(
    Feature,
    name =
      "Selection_Count"
  ) %>%
  dplyr::mutate(
    Repetitions =
      N_REPEATS,
    
    Selection_Percentage =
      100 *
      Selection_Count /
      N_REPEATS
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Selection_Count
    ),
    Feature
  )

write.csv(
  feature_selection_stability,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_RETUNED_Feature_Selection_Stability.csv"
  ),
  row.names = FALSE
)

cat(
  "\nMost consistently selected genera:\n"
)

print(
  head(
    feature_selection_stability,
    30
  )
)

# =============================================================================
# 15) Number of selected genera
# =============================================================================

selected_feature_count_summary <-
  nested_results %>%
  dplyr::filter(
    Model ==
      "Nested ANCOM-BC2 selected"
  ) %>%
  dplyr::summarise(
    Repetitions =
      dplyr::n(),
    
    Mean_Selected =
      mean(
        Number_of_Features,
        na.rm = TRUE
      ),
    
    SD_Selected =
      stats::sd(
        Number_of_Features,
        na.rm = TRUE
      ),
    
    Median_Selected =
      stats::median(
        Number_of_Features,
        na.rm = TRUE
      ),
    
    Minimum_Selected =
      min(
        Number_of_Features,
        na.rm = TRUE
      ),
    
    Maximum_Selected =
      max(
        Number_of_Features,
        na.rm = TRUE
      )
  )

cat(
  "\nSelected feature count:\n"
)

print(
  selected_feature_count_summary
)

write.csv(
  selected_feature_count_summary,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_RETUNED_Selected_Feature_Count_Summary.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 16) Figure 8
# =============================================================================

figure8_data <-
  nested_results %>%
  dplyr::mutate(
    Model = factor(
      Model,
      levels = c(
        "All features",
        "Nested ANCOM-BC2 selected"
      ),
      labels = c(
        "All 447 genera",
        "Nested ANCOM-BC2-selected genera"
      )
    )
  )

p_figure8 <-
  ggplot2::ggplot(
    figure8_data,
    ggplot2::aes(
      x = Model,
      y = ROC_AUC,
      fill = Model
    )
  ) +
  ggplot2::geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.65
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y =
      max(
        figure8_data$ROC_AUC,
        na.rm = TRUE
      ) +
      0.012,
    
    label =
      paste0(
        "Paired Wilcoxon P = ",
        format.pval(
          paired_auc_test$p.value,
          digits = 3,
          eps = 0.001
        )
      ),
    
    size = 4.2
  ) +
  ggplot2::labs(
    title =
      "Predictive Performance Using All Versus ANCOM-BC2-Selected Genera",
    subtitle =
      "ANCOM-BC2 feature selection and selected-model tuning were restricted to each training split",
    x = NULL,
    y = "Test ROC-AUC",
    fill = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 13
  ) +
  ggplot2::theme(
    legend.position =
      "none"
  )

print(
  p_figure8
)

ggplot2::ggsave(
  filename = file.path(
    save_dir,
    "Figure8_All_Features_vs_Nested_ANCOMBC2_RETUNED.png"
  ),
  plot = p_figure8,
  width = 8,
  height = 6,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(
    save_dir,
    "Figure8_All_Features_vs_Nested_ANCOMBC2_RETUNED.pdf"
  ),
  plot = p_figure8,
  width = 8,
  height = 6
)

# =============================================================================
# 17) Session info
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    save_dir,
    "Nested_ANCOMBC2_RETUNED_SessionInfo.txt"
  )
)

cat(
  "\n============================================\n"
)

cat(
  "RETUNED NESTED ANCOM-BC2 ANALYSIS COMPLETE\n"
)

cat(
  "Results saved in:\n",
  save_dir,
  "\n"
)

cat(
  "============================================\n"
)
# =============================================================================
# FIGURE 8: ALL GENERA VS NESTED ANCOM-BC2 SELECTED GENERA
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

results_dir <- "C:/Users/sabir/Downloads/gut microbiota/results_final_RETUNED"

nested_results <- read.csv(
  file.path(
    results_dir,
    "Nested_ANCOMBC2_RETUNED_All_Features_vs_Selected_All_Runs.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Paired test using identical repeated splits
auc_wide <- nested_results %>%
  dplyr::select(
    Seed,
    Model,
    ROC_AUC
  ) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = ROC_AUC
  )

paired_auc_test <- wilcox.test(
  auc_wide$`Nested ANCOM-BC2 selected`,
  auc_wide$`All features`,
  paired = TRUE,
  exact = FALSE
)

figure8_data <- nested_results %>%
  dplyr::mutate(
    Model = factor(
      Model,
      levels = c(
        "All features",
        "Nested ANCOM-BC2 selected"
      ),
      labels = c(
        "All 447 genera",
        "Nested ANCOM-BC2-selected genera"
      )
    )
  )

p8 <- ggplot(
  figure8_data,
  aes(
    x = Model,
    y = ROC_AUC,
    fill = Model
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.65
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(
      figure8_data$ROC_AUC,
      na.rm = TRUE
    ) + 0.012,
    label = paste0(
      "Paired Wilcoxon P ",
      ifelse(
        paired_auc_test$p.value < 0.001,
        "< 0.001",
        paste0(
          "= ",
          format.pval(
            paired_auc_test$p.value,
            digits = 3
          )
        )
      )
    ),
    size = 4.3
  ) +
  labs(
    title = "Predictive Performance Using All Versus ANCOM-BC2-Selected Genera",
    subtitle = stringr::str_wrap(
      paste0(
        "ANCOM-BC2 feature selection and selected-model tuning ",
        "were restricted to each training split"
      ),
      width = 90
    ),
    x = NULL,
    y = "Test ROC-AUC",
    fill = NULL
  ) +
  coord_cartesian(
    ylim = c(
      min(
        figure8_data$ROC_AUC,
        na.rm = TRUE
      ) - 0.015,
      max(
        figure8_data$ROC_AUC,
        na.rm = TRUE
      ) + 0.030
    ),
    clip = "off"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15,
      margin = margin(
        b = 8
      )
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10.5,
      lineheight = 1.1,
      margin = margin(
        b = 12
      )
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 11
    ),
    legend.position = "none",
    plot.margin = margin(
      t = 20,
      r = 25,
      b = 20,
      l = 25
    )
  )

print(p8)

ggsave(
  file.path(
    results_dir,
    "Figure8_All_vs_ANCOMBC2_Final.png"
  ),
  p8,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    results_dir,
    "Figure8_All_vs_ANCOMBC2_Final.pdf"
  ),
  p8,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)
