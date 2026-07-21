# =============================================================================
# FULLY NESTED ANCOM-BC2 FEATURE-SELECTION COMPARISON
#
# Comparison:
# 1. XGBoost using all genus predictors
# 2. XGBoost using genera selected by ANCOM-BC2 on training data only
#
# Feature selection is repeated independently within every train-test split.
# =============================================================================



save_dir <- file.path(
  getwd(),
  "results_final"
)

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
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

tidymodels::tidymodels_prefer()

# =============================================================================
# 2) Load discovery genus-level raw-count data
# =============================================================================

discovery_df <- readRDS(
  file.path(
    model_dir,
    "discovery_genus_training_data.rds"
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

cat("Samples:", nrow(discovery_df), "\n")
cat("Predictors:", length(predictor_names), "\n")
cat("HC:", sum(discovery_df$Group == "HC"), "\n")
cat("RA:", sum(discovery_df$Group == "RA"), "\n")

# =============================================================================
# 3) Repetition settings
# =============================================================================

# Test with 5.
# After successful completion, change only this line to 30.
N_REPEATS <- 30

# Always create the same complete set of 30 seeds.
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
  "Nested_ANCOMBC2_XGBoost_Progress.csv"
)

selected_features_file <- file.path(
  save_dir,
  "Nested_ANCOMBC2_Selected_Features_Progress.csv"
)

# For a fresh 5-repeat test, remove old test files.
# Comment these lines out before changing N_REPEATS to 30,
# so the 30-repeat analysis can resume from the completed five.
#if (N_REPEATS == 5) {
  
#  if (file.exists(progress_file)) {
#    file.remove(progress_file)
#  }
  
#  if (file.exists(selected_features_file)) {
#    file.remove(selected_features_file)
#  }
#}

# =============================================================================
# 4) XGBoost specification
#
# These values match your finalized genus-level XGBoost model.
# mtry is adjusted downward when fewer features are selected.
# =============================================================================

create_xgb_spec <- function(number_of_predictors) {
  
  adjusted_mtry <- min(
    269L,
    as.integer(number_of_predictors)
  )
  
  parsnip::boost_tree(
    mtry = adjusted_mtry,
    trees = 1000,
    min_n = 7,
    tree_depth = 9,
    learn_rate = 0.00788046281566991,
    loss_reduction = 0.0126896100316792,
    sample_size = 0.913793103448276
  ) %>%
    parsnip::set_engine(
      "xgboost",
      eval_metric = "auc",
      nthread = 1
    ) %>%
    parsnip::set_mode("classification")
}

# =============================================================================
# 5) Recipe function
# =============================================================================
create_model_recipe <- function(
    training_data,
    selected_predictors) {
  
  # Build a training table containing only:
  # Sample, Group, and the selected predictors
  reduced_training_data <- training_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(selected_predictors)
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
    ) %>%
    recipes::step_normalize(
      recipes::all_predictors()
    )
}

# =============================================================================
# 6) Model evaluation function
# =============================================================================
fit_and_evaluate_xgb <- function(
    train_data,
    test_data,
    selected_predictors,
    model_name,
    seed_value) {
  
  train_selected <- train_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(selected_predictors)
    )
  
  test_selected <- test_data %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(selected_predictors)
    )
  
  model_recipe <- create_model_recipe(
    training_data = train_selected,
    selected_predictors = selected_predictors
  )
  
  model_spec <- create_xgb_spec(
    number_of_predictors = length(selected_predictors)
  )
  
  model_workflow <- workflows::workflow() %>%
    workflows::add_recipe(model_recipe) %>%
    workflows::add_model(model_spec)
  
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
        levels = c("HC", "RA")
      )
    )
  
  data.frame(
    Seed = seed_value,
    Model = model_name,
    Number_of_Features = length(selected_predictors),
    
    Accuracy = yardstick::accuracy(
      predictions,
      truth = Group,
      estimate = .pred_class
    )$.estimate,
    
    Sensitivity = yardstick::sens(
      predictions,
      truth = Group,
      estimate = .pred_class,
      event_level = "second"
    )$.estimate,
    
    Specificity = yardstick::spec(
      predictions,
      truth = Group,
      estimate = .pred_class,
      event_level = "second"
    )$.estimate,
    
    ROC_AUC = yardstick::roc_auc(
      predictions,
      truth = Group,
      .pred_RA,
      event_level = "second"
    )$.estimate
  )
}


# =============================================================================
# 7) Training-only ANCOM-BC2 feature selection
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
  
  storage.mode(train_predictors) <- "numeric"
  
  # Remove the unclassified genus from DA selection, if present.
  valid_predictors <- colnames(train_predictors)[
    tolower(colnames(train_predictors)) != "unknown"
  ]
  
  train_predictors <- train_predictors[
    ,
    valid_predictors,
    drop = FALSE
  ]
  
  # Training-only prevalence filtering.
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
      prevalence_counts >= minimum_samples
    ]
  )
  
  train_predictors_filt <- train_predictors[
    ,
    retained_predictors,
    drop = FALSE
  ]
  
  cat(
    "  ANCOM-BC2 training genera:",
    ncol(train_predictors_filt),
    "\n"
  )
  
  # ANCOM-BC2 expects taxa in rows and samples in columns.
  ancom_matrix <- t(
    train_predictors_filt
  )
  
  ancom_metadata <- data.frame(
    Sample = train_data$Sample,
    Group = factor(
      train_data$Group,
      levels = c("HC", "RA")
    ),
    row.names = train_data$Sample
  )
  
  # Ensure matrix sample columns match metadata rows exactly.
  colnames(ancom_matrix) <- train_data$Sample
  
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
  
  selected_features <- ancom_results %>%
    dplyr::filter(
      diff_robust_GroupRA == TRUE
    ) %>%
    dplyr::pull(taxon)
  
  # Keep only names that exist in the original ML matrix.
  selected_features <- base::intersect(
    selected_features,
    predictor_names
  )
  
  list(
    selected_features = selected_features,
    ancom_results = ancom_results,
    genera_tested = nrow(ancom_results)
  )
}

# =============================================================================
# 8) Load completed progress
# =============================================================================

if (file.exists(progress_file)) {
  
  comparison_results <- read.csv(
    progress_file,
    stringsAsFactors = FALSE
  )
  
  completed_seeds <- comparison_results %>%
    dplyr::count(Seed) %>%
    dplyr::filter(n >= 2) %>%
    dplyr::pull(Seed)
  
  cat(
    "Completed repetitions found:",
    length(completed_seeds),
    "\n"
  )
  
} else {
  
  comparison_results <- data.frame()
  completed_seeds <- numeric(0)
  
  cat("No completed nested repetitions found.\n")
}

if (file.exists(selected_features_file)) {
  
  selected_feature_results <- read.csv(
    selected_features_file,
    stringsAsFactors = FALSE
  )
  
} else {
  
  selected_feature_results <- data.frame()
}

remaining_seeds <- repeat_seeds[
  !repeat_seeds %in% completed_seeds
]

cat(
  "Repetitions remaining:",
  length(remaining_seeds),
  "\n"
)

# =============================================================================
# 9) Fully nested repeated comparison
# =============================================================================

for (current_seed in remaining_seeds) {
  
  repetition_number <- match(
    current_seed,
    all_repeat_seeds
  )
  
  cat(
    "\n====================================================\n"
  )
  
  cat(
    "Starting nested repetition",
    repetition_number,
    "of",
    N_REPEATS,
    "- seed:",
    current_seed,
    "\n"
  )
  
  # ---------------------------------------------------------------------------
  # Train-test split
  # ---------------------------------------------------------------------------
  
  set.seed(current_seed)
  
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
  # Model 1: All features
  # ---------------------------------------------------------------------------
  
  cat("  Fitting all-features model...\n")
  
  all_features_result <- fit_and_evaluate_xgb(
    train_data = train_rep,
    test_data = test_rep,
    selected_predictors = predictor_names,
    model_name = "All features",
    seed_value = current_seed
  )
  
  cat(
    "  All-features ROC-AUC:",
    round(all_features_result$ROC_AUC, 4),
    "\n"
  )
  
  # ---------------------------------------------------------------------------
  # Training-only ANCOM-BC2
  # ---------------------------------------------------------------------------
  
  cat("  Running training-only ANCOM-BC2...\n")
  
  selection_result <- select_features_ancombc2(
    train_data = train_rep,
    seed_value = current_seed,
    prevalence_threshold = 0.05
  )
  
  selected_features <- selection_result$selected_features
  
  cat(
    "  Robust selected genera:",
    length(selected_features),
    "\n"
  )
  
  # Save the features chosen in this repetition.
  if (length(selected_features) > 0) {
    
    selected_features_current <- data.frame(
      Seed = current_seed,
      Repetition = repetition_number,
      Feature = selected_features,
      stringsAsFactors = FALSE
    )
    
    selected_feature_results <- dplyr::bind_rows(
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
  # Model 2: ANCOM-BC2-selected features
  # ---------------------------------------------------------------------------
  
  if (length(selected_features) >= 2) {
    
    cat("  Fitting ANCOM-selected model...\n")
    
    selected_result <- fit_and_evaluate_xgb(
      train_data = train_rep,
      test_data = test_rep,
      selected_predictors = selected_features,
      model_name = "Nested ANCOM-BC2 selected",
      seed_value = current_seed
    )
    
    cat(
      "  ANCOM-selected ROC-AUC:",
      round(selected_result$ROC_AUC, 4),
      "\n"
    )
    
  } else {
    
    warning(
      paste0(
        "Fewer than two robust ANCOM-BC2 genera were selected ",
        "for seed ",
        current_seed,
        ". Selected-feature model was not fitted."
      )
    )
    
    selected_result <- data.frame(
      Seed = current_seed,
      Model = "Nested ANCOM-BC2 selected",
      Number_of_Features = length(selected_features),
      Accuracy = NA_real_,
      Sensitivity = NA_real_,
      Specificity = NA_real_,
      ROC_AUC = NA_real_
    )
  }
  
  # ---------------------------------------------------------------------------
  # Save after both models
  # ---------------------------------------------------------------------------
  
  current_results <- dplyr::bind_rows(
    all_features_result,
    selected_result
  )
  
  comparison_results <- dplyr::bind_rows(
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
  
  invisible(gc())
}

# =============================================================================
# 10) Final results
# =============================================================================

nested_results <- comparison_results %>%
  dplyr::filter(
    Seed %in% repeat_seeds
  ) %>%
  dplyr::arrange(
    match(Seed, repeat_seeds),
    Model
  )

print(nested_results)

write.csv(
  nested_results,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_All_Features_vs_Selected_All_Runs.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCompleted seeds:",
  dplyr::n_distinct(nested_results$Seed),
  "\n"
)

cat(
  "Total result rows:",
  nrow(nested_results),
  "\n"
)

# =============================================================================
# 11) Performance summary
# =============================================================================

nested_summary <- nested_results %>%
  tidyr::pivot_longer(
    cols = c(
      Accuracy,
      Sensitivity,
      Specificity,
      ROC_AUC
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::group_by(
    Model,
    Metric
  ) %>%
  dplyr::summarise(
    Repetitions = sum(!is.na(Value)),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(
      Value,
      0.25,
      na.rm = TRUE
    ),
    Q3 = quantile(
      Value,
      0.75,
      na.rm = TRUE
    ),
    Minimum = min(Value, na.rm = TRUE),
    Maximum = max(Value, na.rm = TRUE),
    
    CI_Lower = Mean -
      qt(
        0.975,
        df = Repetitions - 1
      ) *
      SD / sqrt(Repetitions),
    
    CI_Upper = Mean +
      qt(
        0.975,
        df = Repetitions - 1
      ) *
      SD / sqrt(Repetitions),
    
    .groups = "drop"
  )

print(nested_summary)

write.csv(
  nested_summary,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_All_Features_vs_Selected_Summary.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 12) Paired ROC-AUC difference across identical splits
# =============================================================================

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

colnames(auc_wide) <- make.names(
  colnames(auc_wide)
)

auc_wide <- auc_wide %>%
  dplyr::mutate(
    AUC_Difference_Selected_minus_All =
      Nested.ANCOM.BC2.selected -
      All.features
  )

print(auc_wide)

write.csv(
  auc_wide,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_Paired_AUC_Differences.csv"
  ),
  row.names = FALSE
)

paired_auc_test <- wilcox.test(
  auc_wide$Nested.ANCOM.BC2.selected,
  auc_wide$All.features,
  paired = TRUE,
  exact = FALSE,
  conf.int = TRUE
)

print(paired_auc_test)

paired_test_results <- data.frame(
  Test = "Paired Wilcoxon signed-rank test",
  Repetitions = sum(
    complete.cases(
      auc_wide$Nested.ANCOM.BC2.selected,
      auc_wide$All.features
    )
  ),
  Median_AUC_Difference_Selected_minus_All =
    median(
      auc_wide$AUC_Difference_Selected_minus_All,
      na.rm = TRUE
    ),
  Wilcoxon_V = unname(
    paired_auc_test$statistic
  ),
  P_Value = paired_auc_test$p.value
)

write.csv(
  paired_test_results,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_Paired_AUC_Test.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 13) Feature-selection stability
# =============================================================================

feature_selection_stability <- selected_feature_results %>%
  dplyr::filter(
    Seed %in% repeat_seeds
  ) %>%
  dplyr::count(
    Feature,
    name = "Selection_Count"
  ) %>%
  dplyr::mutate(
    Repetitions = N_REPEATS,
    Selection_Percentage =
      100 * Selection_Count / N_REPEATS
  ) %>%
  dplyr::arrange(
    dplyr::desc(Selection_Count),
    Feature
  )

print(
  head(
    feature_selection_stability,
    30
  )
)

write.csv(
  feature_selection_stability,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_Feature_Selection_Stability.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 14) Number of selected features per repetition
# =============================================================================

selected_feature_count_summary <- nested_results %>%
  dplyr::filter(
    Model == "Nested ANCOM-BC2 selected"
  ) %>%
  dplyr::summarise(
    Repetitions = dplyr::n(),
    Mean_Selected = mean(
      Number_of_Features,
      na.rm = TRUE
    ),
    SD_Selected = sd(
      Number_of_Features,
      na.rm = TRUE
    ),
    Median_Selected = median(
      Number_of_Features,
      na.rm = TRUE
    ),
    Minimum_Selected = min(
      Number_of_Features,
      na.rm = TRUE
    ),
    Maximum_Selected = max(
      Number_of_Features,
      na.rm = TRUE
    )
  )

print(selected_feature_count_summary)

write.csv(
  selected_feature_count_summary,
  file.path(
    save_dir,
    "Nested_ANCOMBC2_Selected_Feature_Count_Summary.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# Figure 8: All features versus nested ANCOM-BC2-selected features
# =============================================================================

library(ggplot2)
library(dplyr)

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

p_figure8 <- ggplot(
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
  scale_fill_manual(
    values = c(
      "All 447 genera" = "#5B8FF9",
      "Nested ANCOM-BC2-selected genera" = "#E07A5F"
    )
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(figure8_data$ROC_AUC, na.rm = TRUE) + 0.012,
    label = "Paired Wilcoxon P = 1.83 × 10⁻⁶",
    size = 4.2
  ) +
  labs(
    title = "Predictive Performance Using All Versus ANCOM-BC2-Selected Genera",
    subtitle = "ANCOM-BC2 feature selection was performed independently within each training split",
    x = NULL,
    y = "Test ROC-AUC",
    fill = NULL
  ) +
  coord_cartesian(
    ylim = c(
      min(figure8_data$ROC_AUC, na.rm = TRUE) - 0.02,
      max(figure8_data$ROC_AUC, na.rm = TRUE) + 0.025
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text.x = element_text(
      size = 11
    ),
    legend.position = "none"
  )

print(p_figure8)

ggsave(
  filename = file.path(
    save_dir,
    "Figure8_All_Features_vs_Nested_ANCOMBC2.png"
  ),
  plot = p_figure8,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  filename = file.path(
    save_dir,
    "Figure8_All_Features_vs_Nested_ANCOMBC2.pdf"
  ),
  plot = p_figure8,
  width = 8,
  height = 6
)
