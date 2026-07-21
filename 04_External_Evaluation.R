# ============================================================
# SCRIPT 4: EXTERNAL VALIDATION METRICS AND FIGURES
# ============================================================



library(dplyr)
library(yardstick)
library(pROC)
library(ggplot2)

# Load predictions created in the external prediction script
external_results <- read.csv(
  file.path(
    output_path,
    "External_Genus_Predictions.csv"
  ),
  stringsAsFactors = FALSE
)

# Load the full study metadata
metadata_all <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Derive verified disease labels
metadata_selected <- metadata_all %>%
  transmute(
    Sample = Run,
    Sample_Name = `Sample Name`,
    Group = case_when(
      grepl("^RA_", `Sample Name`) ~ "RA",
      grepl("^GUT_", `Sample Name`) ~ "HC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    Sample %in% external_results$Sample
  )

# Check that each prediction has exactly one true label
cat("Predicted samples:", nrow(external_results), "\n")
cat("Matched metadata samples:", nrow(metadata_selected), "\n")
cat("RA samples:", sum(metadata_selected$Group == "RA"), "\n")
cat("HC samples:", sum(metadata_selected$Group == "HC"), "\n")
cat("Missing labels:", sum(is.na(metadata_selected$Group)), "\n")

evaluation_df <- external_results %>%
  left_join(
    metadata_selected,
    by = "Sample"
  )

if (any(is.na(evaluation_df$Group))) {
  stop("One or more processed samples do not have a verified RA/HC label.")
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
    "External_Predictions_With_True_Labels.csv"
  ),
  row.names = FALSE
)


# Accuracy
external_accuracy <- yardstick::accuracy(
  evaluation_df,
  truth = Group,
  estimate = .pred_class
)

# Sensitivity for RA
external_sensitivity <- yardstick::sens(
  evaluation_df,
  truth = Group,
  estimate = .pred_class,
  event_level = "second"
)

# Specificity for HC
external_specificity <- yardstick::spec(
  evaluation_df,
  truth = Group,
  estimate = .pred_class,
  event_level = "second"
)

# ROC-AUC
external_auc <- yardstick::roc_auc(
  evaluation_df,
  truth = Group,
  .pred_RA,
  event_level = "second"
)

print(external_accuracy)
print(external_sensitivity)
print(external_specificity)
print(external_auc)


external_confusion <- yardstick::conf_mat(
  evaluation_df,
  truth = Group,
  estimate = .pred_class
)

print(external_confusion)

roc_obj <- pROC::roc(
  response = evaluation_df$Group,
  predictor = evaluation_df$.pred_RA,
  levels = c("HC", "RA"),
  direction = "<",
  quiet = TRUE
)

auc_ci <- pROC::ci.auc(roc_obj)

cat("ROC-AUC:", as.numeric(pROC::auc(roc_obj)), "\n")
cat(
  "95% CI:",
  as.numeric(auc_ci[1]),
  "-",
  as.numeric(auc_ci[3]),
  "\n"
)


external_metrics_summary <- data.frame(
  Samples = nrow(evaluation_df),
  RA = sum(evaluation_df$Group == "RA"),
  HC = sum(evaluation_df$Group == "HC"),
  Accuracy = external_accuracy$.estimate,
  Sensitivity = external_sensitivity$.estimate,
  Specificity = external_specificity$.estimate,
  ROC_AUC = external_auc$.estimate,
  ROC_AUC_CI_Lower = as.numeric(auc_ci[1]),
  ROC_AUC_CI_Upper = as.numeric(auc_ci[3])
)

print(external_metrics_summary)

write.csv(
  external_metrics_summary,
  file.path(
    output_path,
    "External_Validation_Metrics_Final.csv"
  ),
  row.names = FALSE
)

##############################################################################
#checks
roc_obj_reverse <- pROC::roc(
  response = evaluation_df$Group,
  predictor = evaluation_df$.pred_RA,
  levels = c("HC", "RA"),
  direction = ">",
  quiet = TRUE
)

pROC::auc(roc_obj_reverse)

head(pred_prob)
levels(evaluation_df$Group)
table(evaluation_df$.pred_class)
all(colnames(external_prediction)[-1] ==
      setdiff(colnames(discovery_training), c("Sample","Group")))
summary(discovery_training[, 2:6])
length(shared)
sum(colSums(external_prediction[, discovery_features]) > 0)
#############################################################################
#ROC curve 
library(ggplot2)
library(pROC)

roc_obj <- pROC::roc(
  response = evaluation_df$Group,
  predictor = evaluation_df$.pred_RA,
  levels = c("HC", "RA"),
  direction = "<"
)

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

roc_auc <- as.numeric(pROC::auc(roc_obj))

ggplot(roc_df,
       aes(FPR, TPR)) +
  geom_line(size = 1.2,
            color = "#0072B2") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey60"
  ) +
  annotate(
    "text",
    x = 0.62,
    y = 0.10,
    label = paste0(
      "AUC = ",
      round(roc_auc,3)
    ),
    size = 5
  ) +
  labs(
    title = "External Validation ROC Curve",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_classic(base_size = 14)

ggsave(
  file.path(
    output_path,
    "Figure8A_ROC.pdf"
  ),
  width = 6,
  height = 5
)

ggsave(
  file.path(
    output_path,
    "Figure8A_ROC.png"
  ),
  width = 6,
  height = 5,
  dpi = 600
)

###########################################################################
#confusion matrix 
cm <- as.data.frame(external_confusion$table)

print(cm)
colnames(cm) 
cm <- as.data.frame(external_confusion$table)

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
    aes(label = Freq),
    size = 7
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#0072B2"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "External Validation Confusion Matrix",
    x = "Observed group",
    y = "Predicted group",
    fill = "Samples"
  )

print(p_cm)

ggsave(
  file.path(
    output_path,
    "Figure8B_ConfusionMatrix.pdf"
  ),
  plot = p_cm,
  width = 5,
  height = 5
)

ggsave(
  file.path(
    output_path,
    "Figure8B_ConfusionMatrix.png"
  ),
  plot = p_cm,
  width = 5,
  height = 5,
  dpi = 600
)

##########################################################################
#prediction probabilities
ggplot(
  evaluation_df,
  aes(
    x=Group,
    y=.pred_RA,
    fill=Group
  )
)+
  geom_boxplot(
    alpha=.6,
    outlier.shape=NA
  )+
  geom_jitter(
    width=.15,
    size=2
  )+
  theme_classic(base_size=14)+
  labs(
    title="Predicted RA Probability",
    x="True Group",
    y="Predicted Probability (RA)"
  )

ggsave(
  file.path(
    output_path,
    "Figure8C_Probability.pdf"
  ),
  width=5,
  height=5
)

ggsave(
  file.path(
    output_path,
    "Figure8C_Probability.png"
  ),
  width=5,
  height=5,
  dpi=600
)
#################################################################################
#save confusion matrix table 
write.csv(
  as.data.frame(external_confusion$table),
  file.path(
    output_path,
    "External_Confusion_Matrix_Final.csv"
  ),
  row.names = FALSE
)

#save labeled prediction table 
file.exists(
  file.path(
    output_path,
    "External_Predictions_With_True_Labels.csv"
  )
)

#create one sample-flow table 
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
    "Table_External_Validation_Sample_Flow.csv"
  ),
  row.names = FALSE
)

print(external_sample_flow)
