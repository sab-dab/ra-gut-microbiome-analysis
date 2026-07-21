# =============================================================================
# Permutation importance stability analysis
# =============================================================================



save_dir <- file.path(getwd(), "results_final")



library(tidymodels)
library(vip)
library(dplyr)

tidymodels_prefer()

install.packages("vip") 

packageVersion("vip")

discovery_df <- readRDS(
  file.path(
    model_dir,
    "discovery_genus_training_data.rds"
  )
)

workflow <- readRDS(
  file.path(
    model_dir,
    "finalized_genus_xgboost_workflow.rds"
  )
)

discovery_df$Group <- factor(
  discovery_df$Group,
  levels=c("HC","RA")
)

predictor_names <- setdiff(
  colnames(discovery_df),
  c("Sample","Group")
)

cat(
  "Predictors:",
  length(predictor_names),
  "\n"
)
###############################################################################
# =============================================================================
# 3) Repeated permutation-importance stability analysis
# =============================================================================

# Use 5 for testing.
# After it completes successfully, change this to 30.
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
  "Permutation_Importance_Stability_Progress.csv"
)

# Remove an incomplete test file before restarting the 5-repeat run.
# Comment this block out when you later resume from 5 to 30 repeats.
if (N_REPEATS == 5 && file.exists(progress_file)) {
  file.remove(progress_file)
}

# -----------------------------------------------------------------------------
# Prediction wrapper required by vip
# -----------------------------------------------------------------------------

predict_ra_probability <- function(object, newdata) {
  
  newdata <- as.data.frame(
    newdata,
    check.names = FALSE
  )
  
  # vip supplies predictors but may omit the Sample ID column.
  # The fitted workflow still expects Sample because it has an ID role.
  if (!"Sample" %in% colnames(newdata)) {
    newdata$Sample <- paste0(
      "Permutation_",
      seq_len(nrow(newdata))
    )
  }
  
  # Retain exactly the columns expected by the fitted workflow.
  newdata <- newdata %>%
    dplyr::select(
      Sample,
      dplyr::all_of(predictor_names)
    )
  
  predict(
    object,
    new_data = newdata,
    type = "prob"
  )$.pred_RA
}

# -----------------------------------------------------------------------------
# Load completed repetitions, if present
# -----------------------------------------------------------------------------

if (file.exists(progress_file)) {
  
  permutation_results <- read.csv(
    progress_file,
    stringsAsFactors = FALSE
  )
  
  completed_seeds <- unique(
    permutation_results$Seed
  )
  
  cat(
    "Completed repetitions found:",
    length(completed_seeds),
    "\n"
  )
  
} else {
  
  permutation_results <- data.frame()
  completed_seeds <- numeric(0)
  
  cat("No previous permutation results found.\n")
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
# 4) Repeated train-test splits and permutation importance
# =============================================================================

for (current_seed in remaining_seeds) {
  
  repetition_number <- match(
    current_seed,
    all_repeat_seeds
  )
  
  cat(
    "\nStarting permutation repetition",
    repetition_number,
    "of",
    N_REPEATS,
    "- seed:",
    current_seed,
    "\n"
  )
  
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
  
  set.seed(current_seed)
  
  fitted_workflow_rep <- parsnip::fit(
    workflow,
    data = train_rep
  )
  
  baseline_predictions <- predict(
    fitted_workflow_rep,
    new_data = test_rep,
    type = "prob"
  )
  
  baseline_auc <- yardstick::roc_auc_vec(
    truth = test_rep$Group,
    estimate = baseline_predictions$.pred_RA,
    event_level = "second"
  )
  
  cat(
    "  Baseline ROC-AUC:",
    round(baseline_auc, 4),
    "\n"
  )
  
  set.seed(current_seed)
  
  importance_rep <- vip::vi_permute(
    object = fitted_workflow_rep,
    
    feature_names = predictor_names,
    
    # vip receives the outcome and predictors.
    # Sample is recreated inside the prediction wrapper.
    train = test_rep %>%
      dplyr::select(
        Group,
        dplyr::all_of(predictor_names)
      ),
    
    target = "Group",
    
    metric = "roc_auc",
    
    event_level = "second",
    
    pred_wrapper = predict_ra_probability,
    
    type = "difference",
    
    nsim = 1,
    
    keep = FALSE,
    
    verbose = FALSE,
    
    parallel = FALSE
  )
  
  importance_rep <- as.data.frame(
    importance_rep
  ) %>%
    dplyr::mutate(
      Seed = current_seed,
      Repetition = repetition_number,
      Baseline_ROC_AUC = baseline_auc
    ) %>%
    dplyr::select(
      Seed,
      Repetition,
      Baseline_ROC_AUC,
      Variable,
      Importance
    ) %>%
    dplyr::arrange(
      dplyr::desc(Importance)
    ) %>%
    dplyr::mutate(
      Rank = dplyr::row_number()
    )
  
  permutation_results <- dplyr::bind_rows(
    permutation_results,
    importance_rep
  ) %>%
    dplyr::distinct(
      Seed,
      Variable,
      .keep_all = TRUE
    )
  
  write.csv(
    permutation_results,
    progress_file,
    row.names = FALSE
  )
  
  cat(
    "  Repetition",
    repetition_number,
    "completed and saved.\n"
  )
  
  rm(
    fitted_workflow_rep,
    split_rep,
    train_rep,
    test_rep,
    baseline_predictions,
    importance_rep
  )
  
  invisible(gc())
}

# =============================================================================
# 5) Save completed results
# =============================================================================

permutation_results_final <- permutation_results %>%
  dplyr::filter(
    Seed %in% repeat_seeds
  ) %>%
  dplyr::arrange(
    Repetition,
    Rank
  )

write.csv(
  permutation_results_final,
  file.path(
    save_dir,
    "Permutation_Importance_All_Runs.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCompleted repetitions:",
  dplyr::n_distinct(
    permutation_results_final$Seed
  ),
  "\n"
)

cat(
  "Total importance rows:",
  nrow(permutation_results_final),
  "\n"
)
#######################################################################
# Keep only repetitions requested in the current run
permutation_results_final <- permutation_results %>%
  dplyr::filter(
    Seed %in% repeat_seeds
  ) %>%
  dplyr::arrange(
    Repetition,
    Rank
  )

write.csv(
  permutation_results_final,
  file.path(
    save_dir,
    "Permutation_Importance_All_Runs.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCompleted repetitions:",
  dplyr::n_distinct(
    permutation_results_final$Seed
  ),
  "\n"
)

cat(
  "Total importance rows:",
  nrow(permutation_results_final),
  "\n"
) 

# =============================================================================
# 5) Permutation-importance stability summary
# =============================================================================

permutation_stability_summary <- permutation_results_final %>%
  dplyr::group_by(
    Variable
  ) %>%
  dplyr::summarise(
    Repetitions = dplyr::n(),
    
    Mean_Importance = mean(
      Importance,
      na.rm = TRUE
    ),
    
    SD_Importance = sd(
      Importance,
      na.rm = TRUE
    ),
    
    Median_Importance = median(
      Importance,
      na.rm = TRUE
    ),
    
    Minimum_Importance = min(
      Importance,
      na.rm = TRUE
    ),
    
    Maximum_Importance = max(
      Importance,
      na.rm = TRUE
    ),
    
    Mean_Rank = mean(
      Rank,
      na.rm = TRUE
    ),
    
    Median_Rank = median(
      Rank,
      na.rm = TRUE
    ),
    
    Top_10_Count = sum(
      Rank <= 10,
      na.rm = TRUE
    ),
    
    Top_20_Count = sum(
      Rank <= 20,
      na.rm = TRUE
    ),
    
    Positive_Importance_Count = sum(
      Importance > 0,
      na.rm = TRUE
    ),
    
    Positive_Importance_Percentage =
      100 * Positive_Importance_Count /
      Repetitions,
    
    Top_10_Percentage =
      100 * Top_10_Count /
      Repetitions,
    
    Top_20_Percentage =
      100 * Top_20_Count /
      Repetitions,
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_Importance)
  )

print(
  head(
    permutation_stability_summary,
    30
  )
)

write.csv(
  permutation_stability_summary,
  file.path(
    save_dir,
    "Permutation_Importance_Stability_Summary.csv"
  ),
  row.names = FALSE
) 

# =============================================================================
# 6) Stable top predictors
# =============================================================================

stable_top_predictors <- permutation_stability_summary %>%
  dplyr::filter(
    Mean_Importance > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(Top_20_Percentage),
    dplyr::desc(Mean_Importance)
  ) %>%
  dplyr::slice_head(
    n = 30
  )

print(stable_top_predictors)

write.csv(
  stable_top_predictors,
  file.path(
    save_dir,
    "Permutation_Importance_Top30_Stable_Genera.csv"
  ),
  row.names = FALSE
) 

# =============================================================================
# 7) Figure: top 20 stable permutation predictors
# =============================================================================
p_permutation <- ggplot2::ggplot(
  top20_permutation,
  ggplot2::aes(
    x = Variable,
    y = Mean_Importance
  )
) +
  ggplot2::geom_col(
    width = 0.75,
    fill = "#5A5A5A"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = pmax(
        Mean_Importance - SD_Importance,
        0
      ),
      ymax = Mean_Importance + SD_Importance
    ),
    width = 0.18,
    linewidth = 0.5
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Stable Permutation Importance",
    x = "Genus",
    y = "Mean decrease in test ROC-AUC"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    axis.text.x = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 55
    )
  ) 

ggplot2::ggsave(
  filename = file.path(
    save_dir,
    "Figure_Permutation_Importance_Stability.png"
  ),
  plot = p_permutation,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
) 

