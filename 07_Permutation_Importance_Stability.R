project_dir <- getwd()
# =============================================================================
# Permutation importance stability analysis
# =============================================================================
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
# 1) Packages
# =============================================================================

library(tidymodels)
library(vip)
library(dplyr)
library(ggplot2)

tidymodels::tidymodels_prefer()

cat(
  "vip version:",
  as.character(packageVersion("vip")),
  "\n"
)

# =============================================================================
# 2) Load rebuilt genus dataset + RETUNED workflow
# =============================================================================

discovery_df <- readRDS(
  file.path(
    model_dir,
    "discovery_genus_training_data_rebuilt.rds"
  )
)

workflow_retuned <- readRDS(
  file.path(
    model_dir,
    "finalized_genus_xgboost_workflow_retuned.rds"
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

cat(
  "\nRETUNED XGBoost specification:\n"
)

print(
  workflows::extract_spec_parsnip(
    workflow_retuned
  )
)

# =============================================================================
# 3) Repeated permutation-importance stability settings
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
  "Permutation_Importance_Stability_RETUNED_Progress.csv"
)

# =============================================================================
# 4) Prediction wrapper required by vip
# =============================================================================

predict_ra_probability <- function(
    object,
    newdata) {
  
  newdata <- as.data.frame(
    newdata,
    check.names = FALSE
  )
  
  # vip may omit Sample because it is not a predictor.
  # The fitted workflow still expects this ID-role column.
  if (
    !"Sample" %in%
    colnames(newdata)
  ) {
    
    newdata$Sample <- paste0(
      "Permutation_",
      seq_len(
        nrow(newdata)
      )
    )
  }
  
  # Retain exactly the columns expected by the workflow.
  newdata <- newdata %>%
    dplyr::select(
      Sample,
      dplyr::all_of(
        predictor_names
      )
    )
  
  predict(
    object,
    new_data = newdata,
    type = "prob"
  )$.pred_RA
}

# =============================================================================
# 5) Load completed repetitions if resuming
# =============================================================================

if (
  file.exists(
    progress_file
  )
) {
  
  permutation_results <- read.csv(
    progress_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
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
  
  cat(
    "No previous RETUNED permutation results found.\n"
  )
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
# 6) Repeated train-test splits + permutation importance
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
    "\nStarting permutation repetition",
    repetition_number,
    "of",
    N_REPEATS,
    "- seed:",
    current_seed,
    "\n"
  )
  
  # ---------------------------------------------------------------------------
  # Outer stratified split
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
  # Fit RETUNED specification on this repetition's training partition
  # ---------------------------------------------------------------------------
  
  set.seed(
    current_seed
  )
  
  fitted_workflow_rep <- parsnip::fit(
    workflow_retuned,
    data = train_rep
  )
  
  # ---------------------------------------------------------------------------
  # Baseline held-out ROC-AUC
  # ---------------------------------------------------------------------------
  
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
    round(
      baseline_auc,
      4
    ),
    "\n"
  )
  
  # ---------------------------------------------------------------------------
  # Permutation importance on held-out test partition
  # ---------------------------------------------------------------------------
  
  set.seed(
    current_seed
  )
  
  importance_rep <- vip::vi_permute(
    object = fitted_workflow_rep,
    
    feature_names = predictor_names,
    
    train = test_rep %>%
      dplyr::select(
        Group,
        dplyr::all_of(
          predictor_names
        )
      ),
    
    target = "Group",
    
    metric = "roc_auc",
    
    event_level = "second",
    
    pred_wrapper =
      predict_ra_probability,
    
    type = "difference",
    
    nsim = 1,
    
    keep = FALSE,
    
    verbose = FALSE,
    
    parallel = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Add repetition metadata + ranks
  # ---------------------------------------------------------------------------
  
  importance_rep <- as.data.frame(
    importance_rep
  ) %>%
    dplyr::mutate(
      Seed = current_seed,
      Repetition =
        repetition_number,
      Baseline_ROC_AUC =
        baseline_auc
    ) %>%
    dplyr::select(
      Seed,
      Repetition,
      Baseline_ROC_AUC,
      Variable,
      Importance
    ) %>%
    dplyr::arrange(
      dplyr::desc(
        Importance
      )
    ) %>%
    dplyr::mutate(
      Rank =
        dplyr::row_number()
    )
  
  permutation_results <-
    dplyr::bind_rows(
      permutation_results,
      importance_rep
    ) %>%
    dplyr::distinct(
      Seed,
      Variable,
      .keep_all = TRUE
    )
  
  # Save after every repetition.
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
  
  invisible(
    gc()
  )
}

# =============================================================================
# 7) Save completed results
# =============================================================================

permutation_results_final <-
  permutation_results %>%
  dplyr::filter(
    Seed %in%
      repeat_seeds
  ) %>%
  dplyr::arrange(
    Repetition,
    Rank
  )

write.csv(
  permutation_results_final,
  file.path(
    save_dir,
    "Permutation_Importance_RETUNED_All_Runs.csv"
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
  nrow(
    permutation_results_final
  ),
  "\n"
)

# Must be 30 repetitions x 447 predictors.
stopifnot(
  dplyr::n_distinct(
    permutation_results_final$Seed
  ) == 30
)

stopifnot(
  nrow(
    permutation_results_final
  ) == 30 * 447
)

# =============================================================================
# 8) Baseline ROC-AUC summary across repetitions
# =============================================================================

baseline_summary <-
  permutation_results_final %>%
  dplyr::distinct(
    Seed,
    Repetition,
    Baseline_ROC_AUC
  ) %>%
  dplyr::summarise(
    Repetitions =
      dplyr::n(),
    
    Mean_ROC_AUC =
      mean(
        Baseline_ROC_AUC
      ),
    
    SD_ROC_AUC =
      stats::sd(
        Baseline_ROC_AUC
      ),
    
    Median_ROC_AUC =
      stats::median(
        Baseline_ROC_AUC
      ),
    
    Minimum_ROC_AUC =
      min(
        Baseline_ROC_AUC
      ),
    
    Maximum_ROC_AUC =
      max(
        Baseline_ROC_AUC
      )
  )

cat(
  "\nRETUNED baseline ROC-AUC summary:\n"
)

print(
  baseline_summary
)

write.csv(
  baseline_summary,
  file.path(
    save_dir,
    "Permutation_Importance_RETUNED_Baseline_AUC_Summary.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 9) Permutation-importance stability summary
# =============================================================================

permutation_stability_summary <-
  permutation_results_final %>%
  dplyr::group_by(
    Variable
  ) %>%
  dplyr::summarise(
    
    Repetitions =
      dplyr::n(),
    
    Mean_Importance =
      mean(
        Importance,
        na.rm = TRUE
      ),
    
    SD_Importance =
      stats::sd(
        Importance,
        na.rm = TRUE
      ),
    
    Median_Importance =
      stats::median(
        Importance,
        na.rm = TRUE
      ),
    
    Minimum_Importance =
      min(
        Importance,
        na.rm = TRUE
      ),
    
    Maximum_Importance =
      max(
        Importance,
        na.rm = TRUE
      ),
    
    Mean_Rank =
      mean(
        Rank,
        na.rm = TRUE
      ),
    
    Median_Rank =
      stats::median(
        Rank,
        na.rm = TRUE
      ),
    
    Top_10_Count =
      sum(
        Rank <= 10,
        na.rm = TRUE
      ),
    
    Top_20_Count =
      sum(
        Rank <= 20,
        na.rm = TRUE
      ),
    
    Positive_Importance_Count =
      sum(
        Importance > 0,
        na.rm = TRUE
      ),
    
    Positive_Importance_Percentage =
      100 *
      Positive_Importance_Count /
      Repetitions,
    
    Top_10_Percentage =
      100 *
      Top_10_Count /
      Repetitions,
    
    Top_20_Percentage =
      100 *
      Top_20_Count /
      Repetitions,
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Mean_Importance
    )
  )

cat(
  "\nTop 30 RETUNED permutation predictors:\n"
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
    "Permutation_Importance_Stability_RETUNED_Summary.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 10) Stable top predictors
# =============================================================================

stable_top_predictors <-
  permutation_stability_summary %>%
  dplyr::filter(
    Mean_Importance > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Top_20_Percentage
    ),
    dplyr::desc(
      Mean_Importance
    )
  ) %>%
  dplyr::slice_head(
    n = 30
  )

cat(
  "\nStable top predictors:\n"
)

print(
  stable_top_predictors
)

write.csv(
  stable_top_predictors,
  file.path(
    save_dir,
    "Permutation_Importance_RETUNED_Top30_Stable_Genera.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 11) Operational stable-predictor definition
# Top-20 in >=50% of repeated splits
# =============================================================================

stable_internal_predictors <-
  permutation_stability_summary %>%
  dplyr::filter(
    Top_20_Percentage >= 50
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Top_20_Percentage
    ),
    dplyr::desc(
      Mean_Importance
    )
  )

cat(
  "\n============================================\n"
)

cat(
  "STABLE INTERNAL PREDICTORS\n"
)

cat(
  "Top-20 in >=50% of repetitions\n"
)

cat(
  "============================================\n"
)

print(
  stable_internal_predictors
)

cat(
  "\nNumber of stable predictors:",
  nrow(
    stable_internal_predictors
  ),
  "\n"
)

write.csv(
  stable_internal_predictors,
  file.path(
    save_dir,
    "Permutation_Importance_RETUNED_Stable_Internal_Predictors.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 12) Figure: Top 20 stable permutation predictors
# =============================================================================

top20_permutation <-
  stable_top_predictors %>%
  dplyr::slice_head(
    n = 20
  )

# Order factors so highest importance appears at the top after coord_flip.
top20_permutation$Variable <-
  factor(
    top20_permutation$Variable,
    levels =
      rev(
        top20_permutation$Variable
      )
  )

p_permutation <-
  ggplot2::ggplot(
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
      ymin =
        pmax(
          Mean_Importance -
            SD_Importance,
          0
        ),
      
      ymax =
        Mean_Importance +
        SD_Importance
    ),
    width = 0.18,
    linewidth = 0.5
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title =
      "Stable Permutation Importance - Retuned Model",
    x = "Genus",
    y =
      "Mean decrease in test ROC-AUC"
  ) +
  ggplot2::theme_classic(
    base_size = 13
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      ),
    
    axis.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 11
      ),
    
    axis.text.x =
      ggplot2::element_text(
        size = 11
      ),
    
    plot.margin =
      ggplot2::margin(
        t = 10,
        r = 15,
        b = 10,
        l = 55
      )
  )

print(
  p_permutation
)

ggplot2::ggsave(
  filename = file.path(
    save_dir,
    "Figure7_Permutation_Importance_Stability_RETUNED.png"
  ),
  plot = p_permutation,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    save_dir,
    "Figure7_Permutation_Importance_Stability_RETUNED.pdf"
  ),
  plot = p_permutation,
  width = 7.5,
  height = 6.5,
  units = "in",
  bg = "white"
)

# =============================================================================
# 13) Save session information
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    save_dir,
    "Permutation_Importance_RETUNED_SessionInfo.txt"
  )
)

cat(
  "\n============================================\n"
)

cat(
  "RETUNED PERMUTATION IMPORTANCE COMPLETE\n"
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
# FIGURE 7: STABLE PERMUTATION IMPORTANCE
# =============================================================================

library(dplyr)
library(ggplot2)

results_dir <- "C:/Users/sabir/Downloads/gut microbiota/results_final_RETUNED"

perm_summary <- read.csv(
  file.path(
    results_dir,
    "Permutation_Importance_Stability_RETUNED_Summary.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Same stability ranking used in the analysis:
# positive mean importance, then Top-20 stability, then mean importance
top20_permutation <- perm_summary %>%
  dplyr::filter(
    Mean_Importance > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(Top_20_Percentage),
    dplyr::desc(Mean_Importance)
  ) %>%
  dplyr::slice_head(
    n = 20
  )

# Make names publication-friendly
top20_permutation$Display_Genus <- gsub(
  "_",
  " ",
  top20_permutation$Variable
)

# Highest-ranking genus at top
top20_permutation$Display_Genus <- factor(
  top20_permutation$Display_Genus,
  levels = rev(
    top20_permutation$Display_Genus
  )
)

p7 <- ggplot(
  top20_permutation,
  aes(
    x = Display_Genus,
    y = Mean_Importance
  )
) +
  geom_col(
    width = 0.72,
    fill = "grey35"
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(
        Mean_Importance - SD_Importance,
        0
      ),
      ymax = Mean_Importance + SD_Importance
    ),
    width = 0.16,
    linewidth = 0.45
  ) +
  coord_flip(
    clip = "off"
  ) +
  labs(
    title = "Stable Permutation Importance",
    x = "Genus",
    y = "Mean decrease in test ROC-AUC"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(
        b = 12
      )
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.text.x = element_text(
      size = 10
    ),
    plot.margin = margin(
      t = 20,
      r = 30,
      b = 20,
      l = 30
    )
  )

print(p7)

ggsave(
  file.path(
    results_dir,
    "Figure7_Permutation_Importance_Final.png"
  ),
  p7,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    results_dir,
    "Figure7_Permutation_Importance_Final.pdf"
  ),
  p7,
  width = 9,
  height = 7,
  units = "in",
  bg = "white"
)

# ============================================================
# Stable predictor overlap with external cohort
# ============================================================

stable_file <- paste0(
  "C:/Users/sabir/Downloads/gut microbiota/",
  "results_final_RETUNED/",
  "Permutation_Importance_RETUNED_Stable_Internal_Predictors.csv"
)

external_file <- paste0(
  "C:/Users/sabir/RA_validation/",
  "external_validation_results/",
  "external_genus_counts_combined.rds"
)

stable <- read.csv(
  stable_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

external_genus_counts <- readRDS(
  external_file
)

normalize_genus_name <- function(x) {
  x <- tolower(x)
  x <- gsub("^g__", "", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("[^a-z0-9_]", "", x)
  trimws(x)
}

stable_clean <- normalize_genus_name(
  stable$Variable
)

external_clean <- normalize_genus_name(
  colnames(external_genus_counts)
)

shared_stable <- base::intersect(
  stable_clean,
  external_clean
)

missing_stable <- base::setdiff(
  stable_clean,
  external_clean
)

cat("\nTotal stable predictors:", length(stable_clean), "\n")

cat(
  "Stable predictors represented externally:",
  length(shared_stable),
  "\n"
)

cat(
  "Percentage represented externally:",
  round(
    100 * length(shared_stable) /
      length(stable_clean),
    1
  ),
  "%\n"
)

cat("\nSHARED STABLE PREDICTORS:\n")
print(shared_stable)

cat("\nNOT REPRESENTED EXTERNALLY:\n")
print(missing_stable)



