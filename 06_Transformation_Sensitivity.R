# =============================================================================
# Transformation sensitivity analysis
# Raw-count log transformation vs relative abundance vs CLR
# =============================================================================


save_dir <- file.path(getwd(), "results_final")


library(dplyr)
library(tidyr)
library(tibble)
library(tidymodels)
library(xgboost)

tidymodels::tidymodels_prefer()

# =============================================================================
# 1) Load saved discovery genus data and finalized XGBoost workflow
# =============================================================================

discovery_genus_df <- readRDS(
  file.path(
    model_dir,
    "discovery_genus_training_data.rds"
  )
)

finalized_xgb_workflow <- readRDS(
  file.path(
    model_dir,
    "finalized_genus_xgboost_workflow.rds"
  )
)

discovery_genus_df$Group <- factor(
  discovery_genus_df$Group,
  levels = c("HC", "RA")
)

predictor_names <- base::setdiff(
  colnames(discovery_genus_df),
  c("Sample", "Group")
)

cat("Samples:", nrow(discovery_genus_df), "\n")
cat("Predictors:", length(predictor_names), "\n")
cat("HC:", sum(discovery_genus_df$Group == "HC"), "\n")
cat("RA:", sum(discovery_genus_df$Group == "RA"), "\n")

# =============================================================================
# 2) Transformation functions
# =============================================================================

make_raw_counts <- function(data) {
  
  data
}

make_relative_abundance <- function(data) {
  
  feature_matrix <- as.matrix(
    data[, predictor_names, drop = FALSE]
  )
  
  sample_totals <- rowSums(feature_matrix)
  
  relative_matrix <- sweep(
    feature_matrix,
    1,
    sample_totals,
    "/"
  )
  
  relative_matrix[!is.finite(relative_matrix)] <- 0
  
  result <- as.data.frame(
    relative_matrix,
    check.names = FALSE
  )
  
  result$Sample <- data$Sample
  result$Group <- data$Group
  
  result %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(predictor_names)
    )
}

make_clr <- function(
    data,
    pseudocount = 0.5) {
  
  feature_matrix <- as.matrix(
    data[, predictor_names, drop = FALSE]
  )
  
  log_matrix <- log(
    feature_matrix + pseudocount
  )
  
  clr_matrix <- log_matrix -
    rowMeans(log_matrix)
  
  result <- as.data.frame(
    clr_matrix,
    check.names = FALSE
  )
  
  result$Sample <- data$Sample
  result$Group <- data$Group
  
  result %>%
    dplyr::select(
      Sample,
      Group,
      dplyr::all_of(predictor_names)
    )
}

raw_count_df <- make_raw_counts(
  discovery_genus_df
)

relative_abundance_df <- make_relative_abundance(
  discovery_genus_df
)

clr_df <- make_clr(
  discovery_genus_df,
  pseudocount = 0.5
)

cat(
  "Raw-count range:",
  range(
    as.matrix(raw_count_df[, predictor_names]),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Relative-abundance range:",
  range(
    as.matrix(relative_abundance_df[, predictor_names]),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "CLR range:",
  range(
    as.matrix(clr_df[, predictor_names]),
    na.rm = TRUE
  ),
  "\n"
)

# =============================================================================
# 3) Reuse the same finalized XGBoost model specification
# =============================================================================

final_xgb_spec <- workflows::extract_spec_parsnip(
  finalized_xgb_workflow
)

print(final_xgb_spec)

# =============================================================================
# 4) Create recipes
# =============================================================================

create_recipe <- function(
    training_data,
    transformation) {
  
  base_recipe <- recipes::recipe(
    Group ~ .,
    data = training_data
  ) %>%
    recipes::update_role(
      Sample,
      new_role = "id"
    ) %>%
    recipes::step_zv(
      recipes::all_predictors()
    )
  
  if (transformation == "Raw log(x+1)") {
    
    base_recipe <- base_recipe %>%
      recipes::step_log(
        recipes::all_predictors(),
        offset = 1
      ) %>%
      recipes::step_normalize(
        recipes::all_predictors()
      )
    
  } else {
    
    # Relative abundance and CLR are already transformed.
    # Standardization is learned from training data only.
    base_recipe <- base_recipe %>%
      recipes::step_normalize(
        recipes::all_predictors()
      )
  }
  
  base_recipe
}

# =============================================================================
# 5) Repeated comparison of transformations
# =============================================================================
N_REPEATS <- 30
set.seed(123)

repeat_seeds <- sample(
  1000:999999,
  size = N_REPEATS,
  replace = FALSE
)

transformation_data <- list(
  "Raw log(x+1)" = raw_count_df,
  "Relative abundance" = relative_abundance_df,
  "CLR" = clr_df
)

evaluate_transformation <- function(
    full_data,
    transformation_name,
    seed_value) {
  
  set.seed(seed_value)
  
  split_obj <- rsample::initial_split(
    full_data,
    prop = 0.8,
    strata = Group
  )
  
  train_df <- rsample::training(split_obj)
  test_df  <- rsample::testing(split_obj)
  
  rec <- create_recipe(
    training_data = train_df,
    transformation = transformation_name
  )
  
  wf <- workflows::workflow() %>%
    workflows::add_recipe(rec) %>%
    workflows::add_model(final_xgb_spec)
  
  fitted_wf <- parsnip::fit(
    wf,
    data = train_df
  )
  
  predictions <- test_df %>%
    dplyr::select(Sample, Group) %>%
    dplyr::bind_cols(
      predict(
        fitted_wf,
        new_data = test_df,
        type = "prob"
      ),
      predict(
        fitted_wf,
        new_data = test_df,
        type = "class"
      )
    ) %>%
    dplyr::mutate(
      Group = factor(Group, levels = c("HC", "RA"))
    )
  
  data.frame(
    Seed = seed_value,
    Transformation = transformation_name,
    
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

all_transformation_results <- list()

result_index <- 1

for (current_seed in repeat_seeds) {
  
  cat(
    "\nStarting seed:",
    current_seed,
    "\n"
  )
  
  for (transformation_name in names(transformation_data)) {
    
    cat(
      "  Running:",
      transformation_name,
      "\n"
    )
    
    current_result <- evaluate_transformation(
      full_data = transformation_data[[transformation_name]],
      transformation_name = transformation_name,
      seed_value = current_seed
    )
    
    all_transformation_results[[result_index]] <- current_result
    result_index <- result_index + 1
    
    write.csv(
      dplyr::bind_rows(all_transformation_results),
      file.path(
        save_dir,
        "Transformation_Sensitivity_Progress.csv"
      ),
      row.names = FALSE
    )
    
    invisible(gc())
  }
}

transformation_results <- dplyr::bind_rows(
  all_transformation_results
)

print(transformation_results)

write.csv(
  transformation_results,
  file.path(
    save_dir,
    "Transformation_Sensitivity_All_Runs.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 6) Summarize transformation sensitivity
# =============================================================================

transformation_summary <- transformation_results %>%
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
    Transformation,
    Metric
  ) %>%
  dplyr::summarise(
    Repetitions = dplyr::n(),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Minimum = min(Value, na.rm = TRUE),
    Maximum = max(Value, na.rm = TRUE),
    CI_Lower = Mean -
      qt(0.975, df = Repetitions - 1) *
      SD / sqrt(Repetitions),
    CI_Upper = Mean +
      qt(0.975, df = Repetitions - 1) *
      SD / sqrt(Repetitions),
    .groups = "drop"
  )

print(transformation_summary)

write.csv(
  transformation_summary,
  file.path(
    save_dir,
    "Transformation_Sensitivity_Summary.csv"
  ),
  row.names = FALSE
)
