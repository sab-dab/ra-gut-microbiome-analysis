#############################################################
# Rheumatoid Arthritis Gut Microbiome Analysis
#
# Main Analysis Pipeline
#
# This script coordinates the analyses performed in the
# accompanying manuscript. Individual analysis steps are
# implemented in the supporting R scripts contained in this
# repository.
#
# Execution order:
# 1. 01_DADA2_processing.R
# 2. 05_ANCOMBC2_Discovery.R
# 3. 06_Transformation_Sensitivity.R
# 4. 07_Permutation_Importance_Stability.R
# 5. 08_Nested_ANCOMBC2_Feature_Selection.R
# 6. 03_External_Prediction.R
# 7. 04_External_Evaluation.R
# 8. External_Genus_Differential_Abundance.R 
#
# Author: Sabira Dabeer
# Year: 2026
# =============================================================================
# 0) Setup
# =============================================================================
setwd("C:/Users/sabir/Downloads/gut microbiota")

save_dir <- file.path(getwd(), "results_final")
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

install.packages(c(
  "dplyr",
  "tibble",
  "stringr",
  "tidyr",
  "ggplot2",
  "vegan",
  "yardstick",
  "pROC",
  "pheatmap",
  "xgboost"
))

install.packages("BiocManager")

BiocManager::install("DESeq2") 
install.packages("tidymodels", dependencies = TRUE)
# =============================================================================
# 1) Libraries
# =============================================================================
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(ggplot2)
library(vegan)
library(yardstick)
library(DESeq2)
library(pROC)

library(pheatmap)
library(xgboost)

tidymodels::tidymodels_prefer()
conflicted::conflicts_prefer(dplyr::desc)
# =============================================================================
# 2) Load data
# =============================================================================
asv_profile <- readRDS("1.ASV.profile.rds")
tax_info    <- readRDS("1.taxonomy.info.rds")

asv_profile_df <- asv_profile %>%
  as.data.frame() %>%
  rownames_to_column("ASV")

tax_info_df <- tax_info %>%
  as.data.frame() %>%
  rownames_to_column("ASV")

sample_cols <- base::setdiff(colnames(asv_profile_df), "ASV")
hc_cols <- grep("^HC", sample_cols, value = TRUE)
ra_cols <- grep("^RA", sample_cols, value = TRUE)

cat("HC samples:", length(hc_cols), "\n")
cat("RA samples:", length(ra_cols), "\n")
cat("Total ASVs:", nrow(asv_profile_df), "\n\n")

# =============================================================================
# 3) Dataset suitability and quality checks
# =============================================================================
cat("Taxonomy rows:", nrow(tax_info_df), "\n")
cat("ASV overlap between profile and taxonomy:",
    sum(asv_profile_df$ASV %in% tax_info_df$ASV), "\n\n")

cat("Missing values in ASV profile:", sum(is.na(asv_profile_df)), "\n")
cat("Missing values in taxonomy:", sum(is.na(tax_info_df)), "\n\n")

# =============================================================================
# 4) Sequencing depth
# =============================================================================
asv_mat_counts <- asv_profile_df %>% select(-ASV) %>% as.matrix()

seq_depth <- colSums(asv_mat_counts, na.rm = TRUE)

depth_df <- data.frame(
  Sample = names(seq_depth),
  Depth = as.numeric(seq_depth),
  Group = ifelse(grepl("^HC", names(seq_depth)), "HC", "RA")
)

p_depth <- ggplot(depth_df, aes(x = Group, y = Depth, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Sequencing Depth per Sample",
    x = NULL,
    y = "Total Reads (Sum of ASVs)"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(save_dir, "Figure_Sequencing_Depth.png"),
  plot = p_depth, width = 7, height = 5, dpi = 300
)

write.csv(
  depth_df,
  file.path(save_dir, "Supplementary_Table_S4_Discovery_Sequencing_Depth_By_Sample.csv"),
  row.names = FALSE
)

low_depth_cutoff <- quantile(depth_df$Depth, 0.05)

cat("Low-depth cutoff (5th percentile):", low_depth_cutoff, "\n")
cat("Number of low-depth samples:", sum(depth_df$Depth < low_depth_cutoff), "\n\n")

depth_stats <- depth_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_depth = mean(Depth, na.rm = TRUE),
    median_depth = median(Depth, na.rm = TRUE),
    sd_depth = sd(Depth, na.rm = TRUE),
    min_depth = min(Depth, na.rm = TRUE),
    max_depth = max(Depth, na.rm = TRUE),
    .groups = "drop"
  )

print(depth_stats)
write.csv(
  depth_stats,
  file.path(save_dir, "Table_Sequencing_Depth_Summary.csv"),
  row.names = FALSE
)

# =============================================================================
# 5) Alpha diversity (Richness + Shannon)
# =============================================================================
richness <- colSums(asv_mat_counts > 0, na.rm = TRUE)
shannon  <- vegan::diversity(t(asv_mat_counts), index = "shannon")

alpha_df <- data.frame(
  Sample = names(richness),
  Richness = as.numeric(richness),
  Shannon = as.numeric(shannon),
  Group = ifelse(grepl("^HC", names(richness)), "HC", "RA")
)

p_richness <- ggplot(alpha_df, aes(x = Group, y = Richness, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "ASV Richness by Group",
    x = NULL,
    y = "Number of ASVs"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(save_dir, "Figure_Richness_Boxplot.png"),
  plot = p_richness, width = 7, height = 5, dpi = 300
)

p_shannon <- ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Shannon Diversity by Group",
    x = NULL,
    y = "Shannon Index"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(save_dir, "Figure_Shannon_Boxplot.png"),
  plot = p_shannon, width = 7, height = 5, dpi = 300
)

# =============================================================================
# 5A) Complete alpha-diversity statistical reporting
# =============================================================================

# Ensure disease group is consistently ordered
alpha_df$Group <- factor(
  alpha_df$Group,
  levels = c("HC", "RA")
)

# -----------------------------------------------------------------------------
# Descriptive statistics by group
# -----------------------------------------------------------------------------

alpha_summary <- alpha_df %>%
  tidyr::pivot_longer(
    cols = c(Shannon, Richness),
    names_to = "Diversity_Metric",
    values_to = "Value"
  ) %>%
  dplyr::group_by(Diversity_Metric, Group) %>%
  dplyr::summarise(
    n = sum(!is.na(Value)),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q1 = unname(quantile(Value, 0.25, na.rm = TRUE)),
    Q3 = unname(quantile(Value, 0.75, na.rm = TRUE)),
    IQR = IQR(Value, na.rm = TRUE),
    Minimum = min(Value, na.rm = TRUE),
    Maximum = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(alpha_summary)

write.csv(
  alpha_summary,
  file.path(
    save_dir,
    "Table_Alpha_Diversity_Summary.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Cliff's delta
#
# Negative value:
# metric tends to be lower in RA than HC
#
# Positive value:
# metric tends to be higher in RA than HC
# -----------------------------------------------------------------------------

calculate_cliffs_delta <- function(ra_values, hc_values) {
  
  ra_values <- as.numeric(ra_values)
  hc_values <- as.numeric(hc_values)
  
  ra_values <- ra_values[is.finite(ra_values)]
  hc_values <- hc_values[is.finite(hc_values)]
  
  pairwise_comparison <- outer(
    ra_values,
    hc_values,
    FUN = "-"
  )
  
  delta <- (
    sum(pairwise_comparison > 0) -
      sum(pairwise_comparison < 0)
  ) / length(pairwise_comparison)
  
  return(as.numeric(delta))
}

# -----------------------------------------------------------------------------
# Bootstrap CI for the difference in group medians
#
# Difference is always:
# median RA minus median HC
# -----------------------------------------------------------------------------

bootstrap_median_difference <- function(
    ra_values,
    hc_values,
    n_boot = 2000,
    seed = 123) {
  
  ra_values <- as.numeric(ra_values)
  hc_values <- as.numeric(hc_values)
  
  ra_values <- ra_values[is.finite(ra_values)]
  hc_values <- hc_values[is.finite(hc_values)]
  
  if (length(ra_values) == 0 || length(hc_values) == 0) {
    stop("RA or HC values are empty after removing missing values.")
  }
  
  set.seed(seed)
  
  bootstrap_differences <- replicate(
    n = n_boot,
    expr = {
      
      boot_ra <- sample(
        ra_values,
        size = length(ra_values),
        replace = TRUE
      )
      
      boot_hc <- sample(
        hc_values,
        size = length(hc_values),
        replace = TRUE
      )
      
      median(boot_ra) - median(boot_hc)
    }
  )
  
  observed_difference <-
    median(ra_values) - median(hc_values)
  
  ci_values <- unname(
    quantile(
      bootstrap_differences,
      probs = c(0.025, 0.975),
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  list(
    observed_difference = observed_difference,
    bootstrap_mean = mean(
      bootstrap_differences,
      na.rm = TRUE
    ),
    bootstrap_median = median(
      bootstrap_differences,
      na.rm = TRUE
    ),
    ci_lower = ci_values[1],
    ci_upper = ci_values[2]
  )
}

# -----------------------------------------------------------------------------
# Function to run one diversity analysis
# -----------------------------------------------------------------------------

run_alpha_test <- function(
    metric_name,
    metric_values,
    group_values,
    seed = 123) {
  
  test_data <- data.frame(
    Value = as.numeric(metric_values),
    Group = factor(
      group_values,
      levels = c("HC", "RA")
    )
  ) %>%
    dplyr::filter(
      is.finite(Value),
      !is.na(Group)
    )
  
  hc_values <- test_data$Value[
    test_data$Group == "HC"
  ]
  
  ra_values <- test_data$Value[
    test_data$Group == "RA"
  ]
  
  wilcoxon_result <- wilcox.test(
    ra_values,
    hc_values,
    alternative = "two.sided",
    exact = FALSE,
    correct = TRUE
  )
  
  bootstrap_result <- bootstrap_median_difference(
    ra_values = ra_values,
    hc_values = hc_values,
    n_boot = 2000,
    seed = seed
  )
  
  data.frame(
    Diversity_Metric = metric_name,
    
    HC_n = length(hc_values),
    RA_n = length(ra_values),
    
    HC_Mean = mean(hc_values),
    HC_SD = sd(hc_values),
    HC_Median = median(hc_values),
    HC_Q1 = unname(quantile(hc_values, 0.25)),
    HC_Q3 = unname(quantile(hc_values, 0.75)),
    
    RA_Mean = mean(ra_values),
    RA_SD = sd(ra_values),
    RA_Median = median(ra_values),
    RA_Q1 = unname(quantile(ra_values, 0.25)),
    RA_Q3 = unname(quantile(ra_values, 0.75)),
    
    Median_Difference_RA_minus_HC =
      bootstrap_result$observed_difference,
    
    Median_Difference_CI_Lower =
      bootstrap_result$ci_lower,
    
    Median_Difference_CI_Upper =
      bootstrap_result$ci_upper,
    
    Bootstrap_Mean_Difference =
      bootstrap_result$bootstrap_mean,
    
    Bootstrap_Median_Difference =
      bootstrap_result$bootstrap_median,
    
    Wilcoxon_W =
      unname(wilcoxon_result$statistic),
    
    P_Value =
      wilcoxon_result$p.value,
    
    Cliffs_Delta_RA_vs_HC =
      calculate_cliffs_delta(
        ra_values = ra_values,
        hc_values = hc_values
      )
  )
}

# -----------------------------------------------------------------------------
# Run Shannon and richness analyses
# -----------------------------------------------------------------------------

alpha_test_shannon <- run_alpha_test(
  metric_name = "Shannon diversity",
  metric_values = alpha_df$Shannon,
  group_values = alpha_df$Group,
  seed = 123
)

alpha_test_richness <- run_alpha_test(
  metric_name = "Observed ASV richness",
  metric_values = alpha_df$Richness,
  group_values = alpha_df$Group,
  seed = 456
)

alpha_test_results <- dplyr::bind_rows(
  alpha_test_shannon,
  alpha_test_richness
)

print(alpha_test_results)

write.csv(
  alpha_test_results,
  file.path(
    save_dir,
    "Table_Alpha_Diversity_Tests.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Formatted table for the manuscript
# -----------------------------------------------------------------------------

alpha_test_formatted <- alpha_test_results %>%
  dplyr::mutate(
    
    HC_Median_IQR = sprintf(
      "%.3f (%.3f–%.3f)",
      HC_Median,
      HC_Q1,
      HC_Q3
    ),
    
    RA_Median_IQR = sprintf(
      "%.3f (%.3f–%.3f)",
      RA_Median,
      RA_Q1,
      RA_Q3
    ),
    
    Median_Difference_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      Median_Difference_RA_minus_HC,
      Median_Difference_CI_Lower,
      Median_Difference_CI_Upper
    ),
    
    Wilcoxon_W = round(
      Wilcoxon_W,
      1
    ),
    
    P_Value_Formatted = dplyr::case_when(
      P_Value < 0.001 ~ "<0.001",
      P_Value < 0.01 ~ sprintf("%.3f", P_Value),
      TRUE ~ sprintf("%.4f", P_Value)
    ),
    
    Cliffs_Delta_RA_vs_HC = round(
      Cliffs_Delta_RA_vs_HC,
      3
    )
  ) %>%
  dplyr::select(
    Diversity_Metric,
    HC_n,
    HC_Median_IQR,
    RA_n,
    RA_Median_IQR,
    Median_Difference_95CI,
    Wilcoxon_W,
    P_Value_Formatted,
    Cliffs_Delta_RA_vs_HC
  )

print(alpha_test_formatted)

write.csv(
  alpha_test_formatted,
  file.path(
    save_dir,
    "Table_Alpha_Diversity_Formatted_For_Manuscript.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 6) Genus-level summary and Top 20 genera plot
# =============================================================================
tax_abund <- tax_info_df %>%
  left_join(asv_profile_df, by = "ASV") %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  mutate(Genus = ifelse(is.na(Genus), "g__Unknown", Genus))

genus_abund <- tax_abund %>%
  select(Genus, all_of(sample_cols)) %>%
  group_by(Genus) %>%
  summarise(across(all_of(sample_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

genus_rel <- genus_abund
genus_rel[, sample_cols] <- sweep(
  genus_rel[, sample_cols],
  2,
  colSums(genus_rel[, sample_cols]),
  "/"
)

genus_summary <- genus_rel %>%
  mutate(
    mean_HC = rowMeans(across(all_of(hc_cols)), na.rm = TRUE),
    mean_RA = rowMeans(across(all_of(ra_cols)), na.rm = TRUE),
    diff_RA_minus_HC = mean_RA - mean_HC
  ) %>%
  arrange(desc(abs(diff_RA_minus_HC)))

write.csv(
  genus_summary,
  file.path(save_dir, "Table_Genus_Summary.csv"),
  row.names = FALSE
)

genus_bar <- genus_rel %>%
  mutate(
    mean_HC = rowMeans(across(all_of(hc_cols)), na.rm = TRUE),
    mean_RA = rowMeans(across(all_of(ra_cols)), na.rm = TRUE)
  ) %>%
  select(Genus, mean_HC, mean_RA) %>%
  pivot_longer(
    cols = c(mean_HC, mean_RA),
    names_to = "Group",
    values_to = "MeanRelAbund"
  ) %>%
  mutate(Group = recode(Group, mean_HC = "HC", mean_RA = "RA"))

top_genera <- genus_bar %>%
  group_by(Genus) %>%
  summarise(overall_mean = mean(MeanRelAbund), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  slice(1:20) %>%
  pull(Genus)

genus_bar_top <- genus_bar %>%
  filter(Genus %in% top_genera) %>%
  mutate(Genus = factor(Genus, levels = rev(top_genera)))

p_top20 <- ggplot(genus_bar_top, aes(x = Genus, y = MeanRelAbund, fill = Group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 20 Genera (Mean Relative Abundance): HC vs RA",
    x = "Genus",
    y = "Mean Relative Abundance"
  )

ggsave(
  filename = file.path(save_dir, "Figure_Top20_Genera.png"),
  plot = p_top20, width = 10, height = 7, dpi = 300
)


# =============================================================================
# 7) Shared vs HC-only vs RA-only genera
#   
# =============================================================================

library(dplyr)
library(stringr)
library(ggplot2)
library(ggVennDiagram)

# Install once if needed:
# install.packages("ggVennDiagram")

# -----------------------------------------------------------------------------
# 7.1 Identify ASVs present in each group
# -----------------------------------------------------------------------------

hc_asv_ids <- asv_profile_df %>%
  select(ASV, all_of(hc_cols)) %>%
  mutate(
    any_HC = rowSums(
      across(all_of(hc_cols), ~ .x > 0),
      na.rm = TRUE
    ) > 0
  ) %>%
  filter(any_HC) %>%
  pull(ASV)

ra_asv_ids <- asv_profile_df %>%
  select(ASV, all_of(ra_cols)) %>%
  mutate(
    any_RA = rowSums(
      across(all_of(ra_cols), ~ .x > 0),
      na.rm = TRUE
    ) > 0
  ) %>%
  filter(any_RA) %>%
  pull(ASV)

# -----------------------------------------------------------------------------
# 7.2 Create ASV-to-genus taxonomy map
# -----------------------------------------------------------------------------

asv_genus_map <- tax_info_df %>%
  transmute(
    ASV,
    Genus = str_extract(Taxon, "g__[^;]+")
  ) %>%
  filter(
    !is.na(Genus),
    Genus != "g__",
    Genus != ""
  ) %>%
  distinct(ASV, Genus)

# -----------------------------------------------------------------------------
# 7.3 Identify genera present in HC and RA
# -----------------------------------------------------------------------------

HC_genera <- asv_genus_map %>%
  filter(ASV %in% hc_asv_ids) %>%
  distinct(Genus) %>%
  pull(Genus)

RA_genera <- asv_genus_map %>%
  filter(ASV %in% ra_asv_ids) %>%
  distinct(Genus) %>%
  pull(Genus)

shared_genera <- intersect(HC_genera, RA_genera)
hc_only_genera <- setdiff(HC_genera, RA_genera)
ra_only_genera <- setdiff(RA_genera, HC_genera)

# -----------------------------------------------------------------------------
# 7.4 Save genus counts
# -----------------------------------------------------------------------------

summary_tbl <- tibble(
  Category = c(
    "Shared (HC and RA)",
    "HC only",
    "RA only"
  ),
  Count = c(
    length(shared_genera),
    length(hc_only_genera),
    length(ra_only_genera)
  )
)

print(summary_tbl)

write.csv(
  summary_tbl,
  file.path(save_dir, "Table_Shared_HConly_RAonly_Genera.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7.5 Save complete genus-category list
# -----------------------------------------------------------------------------

genus_category_tbl <- bind_rows(
  tibble(
    Genus = shared_genera,
    Category = "Shared"
  ),
  tibble(
    Genus = hc_only_genera,
    Category = "HC only"
  ),
  tibble(
    Genus = ra_only_genera,
    Category = "RA only"
  )
) %>%
  mutate(
    Genus_Label = str_remove(Genus, "^g__")
  ) %>%
  arrange(Category, Genus_Label)

write.csv(
  genus_category_tbl,
  file.path(save_dir, "Table_All_Shared_HConly_RAonly_Genera.csv"),
  row.names = FALSE
)

#
# -----------------------------------------------------------------------------
# 7.6 Create Venn diagram
# -----------------------------------------------------------------------------
install.packages("ggvenn")
library(ggvenn)

venn_input <- list(
  "Healthy controls (HC)" = HC_genera,
  "Rheumatoid arthritis (RA)" = RA_genera
)

p_shared_venn <- ggvenn(
  venn_input,
  
  # Circle colors: blue and light orange
  fill_color = c(
    "#4A90E2",  # Healthy controls
    "#F5A66F"   # Rheumatoid arthritis
  ),
  
  # Transparency produces a purple overlap
  fill_alpha = 0.55,
  
  stroke_color = c(
    "#2166AC",
    "#D95F02"
  ),
  stroke_size = 1.1,
  
  # Group names
  set_name_color = c(
    "#154C8A",
    "#B34700"
  ),
  set_name_size = 5,
  
  # Counts and percentages
  text_color = "black",
  text_size = 5,
  show_percentage = TRUE,
  digits = 0
) +
  labs(
    title = "Shared and Group-Specific Genera",
    subtitle = paste0(
      "Shared: ", length(shared_genera),
      " | HC only: ", length(hc_only_genera),
      " | RA only: ", length(ra_only_genera)
    )
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    plot.margin = margin(
      t = 15,
      r = 40,
      b = 15,
      l = 40
    )
  )

print(p_shared_venn)

ggsave(
  filename = file.path(
    save_dir,
    "Figure_Genus_Shared_HC_RA_Venn.png"
  ),
  plot = p_shared_venn,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    save_dir,
    "Figure_Genus_Shared_HC_RA_Venn.pdf"
  ),
  plot = p_shared_venn,
  width = 10,
  height = 8,
  bg = "white"
)
# -----------------------------------------------------------------------------
# 7.7 Aggregate ASV counts to genus level
# -----------------------------------------------------------------------------

genus_abundance_df <- asv_profile_df %>%
  inner_join(
    asv_genus_map,
    by = "ASV"
  ) %>%
  select(
    Genus,
    all_of(hc_cols),
    all_of(ra_cols)
  ) %>%
  group_by(Genus) %>%
  summarise(
    across(
      everything(),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 7.8 Calculate abundance and prevalence by group
# -----------------------------------------------------------------------------

hc_summary <- genus_abundance_df %>%
  transmute(
    Genus,
    Group = "HC",
    Total_Abundance = rowSums(
      across(all_of(hc_cols)),
      na.rm = TRUE
    ),
    Mean_Abundance = rowMeans(
      across(all_of(hc_cols)),
      na.rm = TRUE
    ),
    Prevalence_Percent = rowMeans(
      across(all_of(hc_cols), ~ .x > 0),
      na.rm = TRUE
    ) * 100
  )

ra_summary <- genus_abundance_df %>%
  transmute(
    Genus,
    Group = "RA",
    Total_Abundance = rowSums(
      across(all_of(ra_cols)),
      na.rm = TRUE
    ),
    Mean_Abundance = rowMeans(
      across(all_of(ra_cols)),
      na.rm = TRUE
    ),
    Prevalence_Percent = rowMeans(
      across(all_of(ra_cols), ~ .x > 0),
      na.rm = TRUE
    ) * 100
  )

genus_group_summary <- bind_rows(
  hc_summary,
  ra_summary
) %>%
  left_join(
    genus_category_tbl %>%
      select(Genus, Category),
    by = "Genus"
  ) %>%
  group_by(Group) %>%
  mutate(
    Relative_Abundance_Percent =
      100 * Total_Abundance /
      sum(Total_Abundance, na.rm = TRUE)
  ) %>%
  ungroup()

write.csv(
  genus_group_summary,
  file.path(
    save_dir,
    "Table_Genus_Abundance_Prevalence_By_Group.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7.9 Select top genera in each category
# -----------------------------------------------------------------------------

top_shared <- genus_group_summary %>%
  filter(Category == "Shared") %>%
  group_by(Genus, Category) %>%
  summarise(
    Combined_Total_Abundance = sum(
      Total_Abundance,
      na.rm = TRUE
    ),
    Mean_Relative_Abundance_Percent = mean(
      Relative_Abundance_Percent,
      na.rm = TRUE
    ),
    Mean_Prevalence_Percent = mean(
      Prevalence_Percent,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(Mean_Relative_Abundance_Percent),
    desc(Mean_Prevalence_Percent)
  ) %>%
  slice_head(n = 10) %>%
  mutate(Group = "Both")

top_hc_only <- genus_group_summary %>%
  filter(
    Category == "HC only",
    Group == "HC"
  ) %>%
  arrange(
    desc(Relative_Abundance_Percent),
    desc(Prevalence_Percent)
  ) %>%
  slice_head(n = 10) %>%
  transmute(
    Genus,
    Category,
    Group,
    Combined_Total_Abundance = Total_Abundance,
    Mean_Relative_Abundance_Percent =
      Relative_Abundance_Percent,
    Mean_Prevalence_Percent =
      Prevalence_Percent
  )

top_ra_only <- genus_group_summary %>%
  filter(
    Category == "RA only",
    Group == "RA"
  ) %>%
  arrange(
    desc(Relative_Abundance_Percent),
    desc(Prevalence_Percent)
  ) %>%
  slice_head(n = 10) %>%
  transmute(
    Genus,
    Category,
    Group,
    Combined_Total_Abundance = Total_Abundance,
    Mean_Relative_Abundance_Percent =
      Relative_Abundance_Percent,
    Mean_Prevalence_Percent =
      Prevalence_Percent
  )

top_category_genera <- bind_rows(
  top_shared,
  top_hc_only,
  top_ra_only
) %>%
  mutate(
    Genus = str_remove(Genus, "^g__")
  )

print(top_category_genera)

write.csv(
  top_category_genera,
  file.path(
    save_dir,
    "Table_Top_Shared_HConly_RAonly_Genera.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7.9 Most abundant and prevalent genera within each category
# -----------------------------------------------------------------------------

# Shared genera: combine HC and RA information
shared_summary <- genus_group_summary %>%
  filter(Category == "Shared") %>%
  group_by(Genus, Category) %>%
  summarise(
    Mean_Relative_Abundance_Percent =
      mean(Relative_Abundance_Percent, na.rm = TRUE),
    
    Mean_Prevalence_Percent =
      mean(Prevalence_Percent, na.rm = TRUE),
    
    .groups = "drop"
  )

# Top 10 shared genera by relative abundance
top_shared_abundant <- shared_summary %>%
  arrange(desc(Mean_Relative_Abundance_Percent)) %>%
  slice_head(n = 10) %>%
  mutate(
    Ranking_Type = "Most abundant",
    Group = "Shared"
  )

# Top 10 shared genera by prevalence
top_shared_prevalent <- shared_summary %>%
  arrange(desc(Mean_Prevalence_Percent)) %>%
  slice_head(n = 10) %>%
  mutate(
    Ranking_Type = "Most prevalent",
    Group = "Shared"
  )

# Top 10 HC-only genera by relative abundance
top_hc_abundant <- genus_group_summary %>%
  filter(
    Category == "HC only",
    Group == "HC"
  ) %>%
  arrange(desc(Relative_Abundance_Percent)) %>%
  slice_head(n = 10) %>%
  transmute(
    Genus,
    Category,
    Group = "HC only",
    Mean_Relative_Abundance_Percent =
      Relative_Abundance_Percent,
    Mean_Prevalence_Percent =
      Prevalence_Percent,
    Ranking_Type = "Most abundant"
  )

# Top 10 HC-only genera by prevalence
top_hc_prevalent <- genus_group_summary %>%
  filter(
    Category == "HC only",
    Group == "HC"
  ) %>%
  arrange(desc(Prevalence_Percent)) %>%
  slice_head(n = 10) %>%
  transmute(
    Genus,
    Category,
    Group = "HC only",
    Mean_Relative_Abundance_Percent =
      Relative_Abundance_Percent,
    Mean_Prevalence_Percent =
      Prevalence_Percent,
    Ranking_Type = "Most prevalent"
  )

# Top 10 RA-only genera by relative abundance
top_ra_abundant <- genus_group_summary %>%
  filter(
    Category == "RA only",
    Group == "RA"
  ) %>%
  arrange(desc(Relative_Abundance_Percent)) %>%
  slice_head(n = 10) %>%
  transmute(
    Genus,
    Category,
    Group = "RA only",
    Mean_Relative_Abundance_Percent =
      Relative_Abundance_Percent,
    Mean_Prevalence_Percent =
      Prevalence_Percent,
    Ranking_Type = "Most abundant"
  )

# Top 10 RA-only genera by prevalence
top_ra_prevalent <- genus_group_summary %>%
  filter(
    Category == "RA only",
    Group == "RA"
  ) %>%
  arrange(desc(Prevalence_Percent)) %>%
  slice_head(n = 10) %>%
  transmute(
    Genus,
    Category,
    Group = "RA only",
    Mean_Relative_Abundance_Percent =
      Relative_Abundance_Percent,
    Mean_Prevalence_Percent =
      Prevalence_Percent,
    Ranking_Type = "Most prevalent"
  )

# Combine all results
top_category_genera <- bind_rows(
  top_shared_abundant,
  top_shared_prevalent,
  top_hc_abundant,
  top_hc_prevalent,
  top_ra_abundant,
  top_ra_prevalent
) %>%
  mutate(
    Genus = str_remove(Genus, "^g__"),
    Mean_Relative_Abundance_Percent =
      round(Mean_Relative_Abundance_Percent, 4),
    Mean_Prevalence_Percent =
      round(Mean_Prevalence_Percent, 2)
  ) %>%
  select(
    Group,
    Ranking_Type,
    Genus,
    Mean_Relative_Abundance_Percent,
    Mean_Prevalence_Percent
  )

print(top_category_genera, n = 60)

write.csv(
  top_category_genera,
  file.path(
    save_dir,
    "Table_Top_Abundant_Prevalent_Genera_By_Category.csv"
  ),
  row.names = FALSE
)
# -----------------------------------------------------------------------------
# 7.10 Final count check
# -----------------------------------------------------------------------------

cat(
  "\nShared genera:", length(shared_genera),
  "\nHC-only genera:", length(hc_only_genera),
  "\nRA-only genera:", length(ra_only_genera),
  "\n"
)
# =============================================================================
# 8) Community-level ecology: ASV filtering + Bray-Curtis + PCoA + PERMANOVA
# =============================================================================
asv_mat <- t(as.matrix(asv_profile))   # rows = samples, cols = ASVs

meta <- data.frame(
  Sample = rownames(asv_mat),
  Group  = ifelse(
    rownames(asv_mat) %in% hc_cols, "HC",
    ifelse(rownames(asv_mat) %in% ra_cols, "RA", NA)
  )
) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("HC", "RA")))

asv_mat <- asv_mat[meta$Sample, , drop = FALSE]

n_samples <- nrow(asv_mat)
prev <- colSums(asv_mat > 0)
keep_asvs <- names(prev[prev >= 0.05 * n_samples])
asv_mat_filt <- asv_mat[, keep_asvs, drop = FALSE]

# =============================================================================
# 8A) Data dimensionality and sparsity before and after prevalence filtering
# =============================================================================

calculate_sparsity <- function(mat, stage_name) {
  
  total_cells <- length(mat)
  zero_cells  <- sum(mat == 0, na.rm = TRUE)
  
  data.frame(
    Stage = stage_name,
    Samples = nrow(mat),
    Features_ASVs = ncol(mat),
    Total_Values = total_cells,
    Zero_Values = zero_cells,
    Nonzero_Values = total_cells - zero_cells,
    Zero_Percentage = round(100 * zero_cells / total_cells, 3),
    Nonzero_Percentage = round(
      100 * (total_cells - zero_cells) / total_cells,
      3
    )
  )
}

sparsity_before <- calculate_sparsity(
  asv_mat,
  "Original ASV matrix"
)

sparsity_after <- calculate_sparsity(
  asv_mat_filt,
  "After 5% prevalence filtering"
)

sparsity_summary <- bind_rows(
  sparsity_before,
  sparsity_after
)

print(sparsity_summary)

write.csv(
  sparsity_summary,
  file.path(save_dir, "Table_Data_Sparsity_Summary.csv"),
  row.names = FALSE
)

# Save the prevalence of every ASV
asv_prevalence_table <- data.frame(
  ASV = colnames(asv_mat),
  Number_of_Samples_Present = colSums(asv_mat > 0, na.rm = TRUE),
  Prevalence_Percentage = round(
    100 * colSums(asv_mat > 0, na.rm = TRUE) / nrow(asv_mat),
    3
  ),
  Retained_After_5pct_Filter = colnames(asv_mat) %in% keep_asvs
)

write.csv(
  asv_prevalence_table,
  file.path(save_dir, "Table_ASV_Prevalence_and_Filtering.csv"),
  row.names = FALSE
)

cat("\nASV/SPARSITY SUMMARY\n")
cat("ASVs before filtering:", ncol(asv_mat), "\n")
cat("ASVs after filtering:", ncol(asv_mat_filt), "\n")
cat(
  "Zero percentage before filtering:",
  sparsity_before$Zero_Percentage,
  "%\n"
)
cat(
  "Zero percentage after filtering:",
  sparsity_after$Zero_Percentage,
  "%\n\n"
)

#########################################################################
install.packages("vegan")   # only once, if not installed
library(vegan)
asv_rel <- sweep(asv_mat_filt, 1, rowSums(asv_mat_filt), "/")
asv_rel[is.na(asv_rel)] <- 0

dist_bc <- vegdist(asv_rel, method = "bray")

set.seed(123)
perm <- adonis2(dist_bc ~ Group, data = meta, permutations = 999)
print(perm)

# =============================================================================
# Complete PERMANOVA reporting
# =============================================================================

perm_table <- as.data.frame(perm)
perm_table$Term <- rownames(perm_table)
rownames(perm_table) <- NULL

perm_table <- perm_table %>%
  relocate(Term)

print(perm_table)

write.csv(
  perm_table,
  file.path(save_dir, "Table_PERMANOVA_Results.csv"),
  row.names = FALSE
)

# =============================================================================
# Complete PERMDISP reporting
# betadisper() is the R implementation; the statistical method is PERMDISP.
# =============================================================================

bd <- vegan::betadisper(
  dist_bc,
  group = meta$Group,
  bias.adjust = TRUE
)

bd_anova <- anova(bd)

set.seed(123)
bd_perm <- vegan::permutest(
  bd,
  permutations = 999
)

# Standard ANOVA table
permdisp_anova_table <- as.data.frame(bd_anova)
permdisp_anova_table$Term <- rownames(permdisp_anova_table)
rownames(permdisp_anova_table) <- NULL

permdisp_anova_table <- permdisp_anova_table %>%
  relocate(Term)

print(permdisp_anova_table)

write.csv(
  permdisp_anova_table,
  file.path(save_dir, "Table_PERMDISP_ANOVA_Results.csv"),
  row.names = FALSE
)

# Permutation-test table
permdisp_permutation_table <- as.data.frame(bd_perm$tab)
permdisp_permutation_table$Term <-
  rownames(permdisp_permutation_table)
rownames(permdisp_permutation_table) <- NULL

permdisp_permutation_table <-
  permdisp_permutation_table %>%
  relocate(Term)

print(permdisp_permutation_table)

write.csv(
  permdisp_permutation_table,
  file.path(
    save_dir,
    "Table_PERMDISP_Permutation_Results.csv"
  ),
  row.names = FALSE
)

# Distance from each sample to its group centroid
distance_to_centroid_df <- data.frame(
  Sample = names(bd$distances),
  Distance_to_Centroid = as.numeric(bd$distances)
)

distance_to_centroid_df <- distance_to_centroid_df %>%
  left_join(
    meta %>% select(Sample, Group),
    by = "Sample"
  )

write.csv(
  distance_to_centroid_df,
  file.path(
    save_dir,
    "Table_PERMDISP_Distance_to_Centroid_By_Sample.csv"
  ),
  row.names = FALSE
)

# Group-level summary of dispersion
dispersion_group_summary <- distance_to_centroid_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    Mean_Distance_to_Centroid =
      mean(Distance_to_Centroid, na.rm = TRUE),
    SD_Distance_to_Centroid =
      sd(Distance_to_Centroid, na.rm = TRUE),
    Median_Distance_to_Centroid =
      median(Distance_to_Centroid, na.rm = TRUE),
    Q1_Distance_to_Centroid =
      quantile(Distance_to_Centroid, 0.25, na.rm = TRUE),
    Q3_Distance_to_Centroid =
      quantile(Distance_to_Centroid, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(dispersion_group_summary)
more_dispersed_group <- dispersion_group_summary %>%
  slice_max(
    order_by = Mean_Distance_to_Centroid,
    n = 1,
    with_ties = FALSE
  ) %>%
  pull(Group)

cat(
  "Group with greater multivariate dispersion:",
  as.character(more_dispersed_group),
  "\n"
)
write.csv(
  dispersion_group_summary,
  file.path(
    save_dir,
    "Table_PERMDISP_Distance_to_Centroid_By_Group.csv"
  ),
  row.names = FALSE
)

# Professional dispersion figure
p_dispersion <- ggplot(
  distance_to_centroid_df,
  aes(
    x = Group,
    y = Distance_to_Centroid,
    fill = Group
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 1
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Multivariate Dispersion by Disease Group",
    x = "Disease group",
    y = "Bray–Curtis distance to group centroid"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(
    save_dir,
    "Figure_PERMDISP_Distance_to_Centroid.png"
  ),
  plot = p_dispersion,
  width = 7,
  height = 5,
  dpi = 300
)

pcoa_res <- cmdscale(dist_bc, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  Sample = rownames(asv_rel),
  PC1    = pcoa_res$points[, 1],
  PC2    = pcoa_res$points[, 2]
) %>%
  left_join(meta, by = "Sample")

eig_vals <- pcoa_res$eig
var_explained <- eig_vals / sum(eig_vals[eig_vals > 0]) * 100

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    x = paste0("PCoA1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PCoA2 (", round(var_explained[2], 2), "%)"),
    color = "Group",
    title = "PCoA of ASV-level Gut Microbiome (Bray-Curtis)"
  )

ggsave(
  filename = file.path(save_dir, "Figure_PCoA_BrayCurtis.png"),
  plot = p_pcoa, width = 7, height = 6, dpi = 300
)

# =============================================================================
# 9) Differential abundance: DESeq2 (ASV level)
# =============================================================================
count_mat <- as.matrix(asv_profile)
count_mat <- count_mat[, meta$Sample, drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ Group
)

dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

res <- results(dds, contrast = c("Group", "RA", "HC"))

res_df <- as.data.frame(res) %>%
  mutate(ASV = rownames(.)) %>%
  arrange(padj)

tax_df <- as.data.frame(tax_info) %>%
  rownames_to_column(var = "ASV")

res_annot <- res_df %>%
  left_join(tax_df, by = "ASV")

final_deseq2 <- res_annot %>%
  filter(!is.na(padj)) %>%
  mutate(
    Direction = case_when(
      log2FoldChange > 0 ~ "Higher in RA",
      log2FoldChange < 0 ~ "Higher in HC",
      TRUE ~ "No change"
    ),
    Genus = str_extract(Taxon, "g__[^;]+"),
    Species = str_extract(Taxon, "s__[^;]+")
  ) %>%
  arrange(padj)

sig_final <- final_deseq2 %>% filter(padj < 0.05)

write.csv(
  final_deseq2,
  file.path(save_dir, "DESeq2_All_Annotated.csv"),
  row.names = FALSE
)
write.csv(
  sig_final,
  file.path(save_dir, "DESeq2_Significant_padj0.05.csv"),
  row.names = FALSE
)

top_deseq2 <- sig_final %>%
  select(ASV, Genus, Species, log2FoldChange, padj, Direction) %>%
  slice(1:20)

write.csv(
  top_deseq2,
  file.path(save_dir, "Table_Top20_DESeq2_Significant_ASVs.csv"),
  row.names = FALSE
)

volcano_df <- final_deseq2 %>%
  mutate(
    neglog10_padj = -log10(padj),
    category = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Higher in RA",
      !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Higher in HC",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = neglog10_padj, color = category)
) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(
    "Higher in RA" = "#D55E00",
    "Higher in HC" = "#0072B2",
    "Not significant" = "gray70"
  )) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Differentially Abundant ASVs (RA vs HC)",
    x = "Log2 Fold Change",
    y = expression(-log[10]("adjusted p-value")),
    color = NULL
  )

ggsave(
  filename = file.path(save_dir, "Figure_DESeq2_Volcano.png"),
  plot = p_volcano, width = 8, height = 6, dpi = 300
)

# =============================================================================

# Add group-specific abundance and prevalence to significant DESeq2 ASVs
# =============================================================================
library(dplyr)
library(tibble)
library(stringr)

final_deseq2 <- read.csv(
  file.path(save_dir, "DESeq2_All_Annotated.csv"),
  check.names = FALSE
)

sig_final <- read.csv(
  file.path(save_dir, "DESeq2_Significant_padj0.05.csv"),
  check.names = FALSE
)
# Sample groups
hc_samples <- meta$Sample[meta$Group == "HC"]
ra_samples <- meta$Sample[meta$Group == "RA"]

# Calculate abundance and prevalence for each ASV by group
asv_group_summary <- tibble(
  ASV = rownames(count_mat),
  
  Mean_Abundance_HC = rowMeans(
    count_mat[, hc_samples, drop = FALSE],
    na.rm = TRUE
  ),
  
  Mean_Abundance_RA = rowMeans(
    count_mat[, ra_samples, drop = FALSE],
    na.rm = TRUE
  ),
  
  Median_Abundance_HC = apply(
    count_mat[, hc_samples, drop = FALSE],
    1,
    median,
    na.rm = TRUE
  ),
  
  Median_Abundance_RA = apply(
    count_mat[, ra_samples, drop = FALSE],
    1,
    median,
    na.rm = TRUE
  ),
  
  Prevalence_HC_Percent = rowMeans(
    count_mat[, hc_samples, drop = FALSE] > 0,
    na.rm = TRUE
  ) * 100,
  
  Prevalence_RA_Percent = rowMeans(
    count_mat[, ra_samples, drop = FALSE] > 0,
    na.rm = TRUE
  ) * 100
)

# Create complete supplementary table of significant ASVs
supplementary_deseq2_table <- sig_final %>%
  left_join(
    asv_group_summary,
    by = "ASV"
  ) %>%
  transmute(
    ASV,
    Taxonomic_Annotation = Taxon,
    Genus,
    Species,
    Base_Mean = baseMean,
    Log2_Fold_Change = log2FoldChange,
    Standard_Error = lfcSE,
    Wald_Statistic = stat,
    Raw_P_Value = pvalue,
    Adjusted_P_Value_FDR = padj,
    Direction,
    Mean_Abundance_HC,
    Mean_Abundance_RA,
    Median_Abundance_HC,
    Median_Abundance_RA,
    Prevalence_HC_Percent,
    Prevalence_RA_Percent
  ) %>%
  arrange(Adjusted_P_Value_FDR)
colnames(supplementary_deseq2_table) <- c(
  "ASV",
  "Taxonomic annotation",
  "Genus",
  "Species",
  "Base mean",
  "Log2 fold change",
  "Standard error",
  "Wald statistic",
  "Raw P value",
  "Adjusted P value (FDR)",
  "Direction",
  "Mean abundance (HC)",
  "Mean abundance (RA)",
  "Median abundance (HC)",
  "Median abundance (RA)",
  "Prevalence (%) HC",
  "Prevalence (%) RA"
)
print(supplementary_deseq2_table)

write.csv(
  supplementary_deseq2_table,
  file.path(
    save_dir,
    "Supplementary_Table_DESeq2_Significant_ASVs.csv"
  ),
  row.names = FALSE
)
# =============================================================================
# 9A) Heatmap of top differentially abundant ASVs
# =============================================================================
top_n_heat <- 20

top_heat_asvs <- sig_final %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:top_n_heat) %>%
  pull(ASV)

heat_mat <- asv_mat[, top_heat_asvs, drop = FALSE]
heat_mat_log <- log10(heat_mat + 1)
heat_mat_plot <- t(heat_mat_log)

heat_tax <- tax_info_df %>%
  filter(ASV %in% rownames(heat_mat_plot)) %>%
  mutate(
    Genus = str_extract(Taxon, "g__[^;]+"),
    Genus = ifelse(is.na(Genus), "Unknown", str_remove(Genus, "g__")),
    short_id = str_sub(ASV, 1, 8),
    heat_label = paste0(Genus, " (", short_id, ")")
  ) %>%
  select(ASV, heat_label)

row_labels <- heat_tax$heat_label
names(row_labels) <- heat_tax$ASV
rownames(heat_mat_plot) <- row_labels[rownames(heat_mat_plot)]

annotation_col <- data.frame(Group = meta$Group)
rownames(annotation_col) <- meta$Sample

common_samples <- base::intersect(colnames(heat_mat_plot), rownames(annotation_col))
heat_mat_plot <- heat_mat_plot[, common_samples, drop = FALSE]
annotation_col <- annotation_col[common_samples, , drop = FALSE]

png(
  filename = file.path(save_dir, "Figure_Heatmap_Top20_DA_ASVs_full.png"),
  width = 1800, height = 1200, res = 200
)

pheatmap(
  heat_mat_plot,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize_row = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of the Top 20 Differentially Abundant ASVs Between RA and HC"
)

dev.off()

set.seed(123)

meta_sub <- meta %>%
  group_by(Group) %>%
  sample_n(size = min(100, n())) %>%
  ungroup()

samples_sub <- meta_sub$Sample

heat_mat_plot_sub <- heat_mat_plot[, samples_sub, drop = FALSE]

annotation_col_sub <- data.frame(Group = meta_sub$Group)
rownames(annotation_col_sub) <- samples_sub

png(
  filename = file.path(save_dir, "Figure_Heatmap_Top20_DA_ASVs_sampled.png"),
  width = 1800, height = 1200, res = 200
)

pheatmap(
  heat_mat_plot_sub,
  scale = "row",
  annotation_col = annotation_col_sub,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize_row = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of the Top 20 Differentially Abundant ASVs Between RA and HC"
)

dev.off()

# # =============================================================================
# 10) ASV-level Machine Learning: Baseline model comparison
#     Models: LASSO, Random Forest, SVM-RBF, XGBoost
# =============================================================================
library(tidymodels)

asv_df <- as.data.frame(asv_mat_filt) %>%
  mutate(Sample = rownames(.)) %>%
  relocate(Sample)

ml_df <- asv_df %>%
  left_join(meta, by = "Sample") %>%
  filter(!is.na(Group))

ml_df$Group <- factor(ml_df$Group, levels = c("HC", "RA"))

set.seed(123)
split <- initial_split(ml_df, prop = 0.8, strata = Group)
train_data <- training(split)
test_data  <- testing(split)

# Shared preprocessing recipe
rec_asv <- recipe(Group ~ ., data = train_data) %>%
  update_role(Sample, new_role = "id") %>%
  step_zv(all_predictors()) %>%
  step_log(all_predictors(), offset = 1) %>%
  step_normalize(all_predictors())

set.seed(123)
folds_asv <- vfold_cv(train_data, v = 5, strata = Group)

metric_fn <- yardstick::metric_set(
  yardstick::roc_auc,
  yardstick::accuracy,
  yardstick::sens,
  yardstick::spec
)


# -----------------------------------------------------------------------------
# Model specs
# -----------------------------------------------------------------------------

# 1) LASSO logistic regression
lasso_spec <- logistic_reg(
  penalty = tune(),
  mixture = 1
) %>%
  set_engine("glmnet")

# 2) Random Forest
rf_spec <- rand_forest(
  mtry = tune(),
  min_n = tune(),
  trees = 500
) %>%
  set_engine("ranger") %>%
  set_mode("classification")

# 3) SVM-RBF
svm_spec <- svm_rbf(
  cost = tune(),
  rbf_sigma = tune()
) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# 4) XGBoost
xgb_spec <- boost_tree(
  trees = 1000,
  learn_rate = tune(),
  mtry = tune(),
  tree_depth = tune(),
  min_n = tune(),
  loss_reduction = tune(),
  sample_size = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# -----------------------------------------------------------------------------
# Workflows
# -----------------------------------------------------------------------------
wf_lasso <- workflow() %>% add_recipe(rec_asv) %>% add_model(lasso_spec)
wf_rf    <- workflow() %>% add_recipe(rec_asv) %>% add_model(rf_spec)
wf_svm   <- workflow() %>% add_recipe(rec_asv) %>% add_model(svm_spec)
wf_xgb   <- workflow() %>% add_recipe(rec_asv) %>% add_model(xgb_spec)

# -----------------------------------------------------------------------------
# Tuning grids
# -----------------------------------------------------------------------------
lasso_grid <- grid_regular(
  penalty(range = c(-5, 0)),
  levels = 25
)

rf_grid <- grid_regular(
  mtry(range = c(5L, min(150L, ncol(train_data) - 2L))),
  min_n(range = c(2L, 30L)),
  levels = 5
)

svm_grid <- grid_space_filling(
  cost(),
  rbf_sigma(),
  size = 20
)

xgb_grid <- grid_latin_hypercube(
  learn_rate(range = c(-5, -1)),
  mtry(range = c(5L, min(300L, ncol(train_data) - 2L))),
  tree_depth(range = c(2L, 10L)),
  min_n(range = c(2L, 30L)),
  loss_reduction(range = c(-5, 1)),
  sample_size = sample_prop(range = c(0.5, 1.0)),
  size = 30
  
)

# -----------------------------------------------------------------------------
# Tune models
# -----------------------------------------------------------------------------
install.packages("glmnet")
library(glmnet) 

set.seed(123)
tuned_lasso <- tune_grid(
  wf_lasso,
  resamples = folds_asv,
  grid = lasso_grid,
  metrics = metric_fn
)

install.packages("ranger")
library(ranger) 

set.seed(123)
tuned_rf <- tune_grid(
  wf_rf,
  resamples = folds_asv,
  grid = rf_grid,
  metrics = metric_fn
)

install.packages("kernlab")
library(kernlab)
set.seed(123)
tuned_svm <- tune_grid(
  wf_svm,
  resamples = folds_asv,
  grid = svm_grid,
  metrics = metric_fn
)


library(xgboost)

set.seed(123)
tuned_xgb <- tune_grid(
  wf_xgb,
  resamples = folds_asv,
  grid = xgb_grid,
  metrics = metric_fn
)

# -----------------------------------------------------------------------------
# Select best hyperparameters
# -----------------------------------------------------------------------------
# for Lasso 
lasso_metrics <- collect_metrics(tuned_lasso)

best_lasso <- lasso_metrics %>%
  dplyr::filter(.metric == "roc_auc") %>%
  dplyr::arrange(dplyr::desc(mean)) %>%
  dplyr::slice(1) %>%
  dplyr::select(penalty)

best_lasso

#for Random Forest
rf_metrics <- collect_metrics(tuned_rf)

best_rf <- rf_metrics %>%
  dplyr::filter(.metric == "roc_auc") %>%
  dplyr::arrange(dplyr::desc(mean)) %>%
  dplyr::slice(1) %>%
  dplyr::select(mtry, min_n)

best_rf

#for SVM 
svm_metrics <- collect_metrics(tuned_svm)

best_svm <- svm_metrics %>%
  dplyr::filter(.metric == "roc_auc") %>%
  dplyr::arrange(dplyr::desc(mean)) %>%
  dplyr::slice(1) %>%
  dplyr::select(cost, rbf_sigma)

best_svm

# for XGBoost
xgb_metrics <- collect_metrics(tuned_xgb)

best_xgb <- xgb_metrics %>%
  dplyr::filter(.metric == "roc_auc") %>%
  dplyr::arrange(dplyr::desc(mean)) %>%
  dplyr::slice(1) %>%
  dplyr::select(learn_rate, mtry, tree_depth, min_n, loss_reduction, sample_size)

best_xgb


print(best_lasso)
print(best_rf)
print(best_svm)
print(best_xgb)

# -----------------------------------------------------------------------------
# Finalize and fit models
# -----------------------------------------------------------------------------
final_wf_lasso <- finalize_workflow(wf_lasso, best_lasso)
final_wf_rf    <- finalize_workflow(wf_rf, best_rf)
final_wf_svm   <- finalize_workflow(wf_svm, best_svm)
final_wf_xgb   <- finalize_workflow(wf_xgb, best_xgb)

final_fit_lasso <- fit(final_wf_lasso, data = train_data)
final_fit_rf    <- fit(final_wf_rf, data = train_data)
final_fit_svm   <- fit(final_wf_svm, data = train_data)
final_fit_xgb   <- fit(final_wf_xgb, data = train_data)


# =============================================================================
# Evaluation and 95% confidence intervals for all four models
# =============================================================================

get_model_predictions <- function(
    final_fit,
    test_data,
    model_name) {
  
  test_data %>%
    select(Sample, Group) %>%
    bind_cols(
      predict(
        final_fit,
        new_data = test_data,
        type = "prob"
      ),
      predict(
        final_fit,
        new_data = test_data,
        type = "class"
      )
    ) %>%
    mutate(
      Group = factor(Group, levels = c("HC", "RA")),
      Model = model_name
    )
}

pred_lasso <- get_model_predictions(
  final_fit_lasso,
  test_data,
  "LASSO Logistic"
)

pred_rf <- get_model_predictions(
  final_fit_rf,
  test_data,
  "Random Forest"
)

pred_svm <- get_model_predictions(
  final_fit_svm,
  test_data,
  "SVM-RBF"
)

pred_xgb <- get_model_predictions(
  final_fit_xgb,
  test_data,
  "XGBoost"
)

# -------------------------------------------------------------------------
# Calculate point estimates
# -------------------------------------------------------------------------

calculate_model_metrics <- function(pred_df) {
  
  data.frame(
    Accuracy = yardstick::accuracy(
      pred_df,
      truth = Group,
      estimate = .pred_class
    )$.estimate,
    
    Sensitivity = yardstick::sens(
      pred_df,
      truth = Group,
      estimate = .pred_class,
      event_level = "second"
    )$.estimate,
    
    Specificity = yardstick::spec(
      pred_df,
      truth = Group,
      estimate = .pred_class,
      event_level = "second"
    )$.estimate,
    
    ROC_AUC = yardstick::roc_auc(
      pred_df,
      truth = Group,
      .pred_RA,
      event_level = "second"
    )$.estimate
  )
}

# -------------------------------------------------------------------------
# Stratified bootstrap:
# Samples HC and RA separately so both classes occur in every bootstrap sample.
# -------------------------------------------------------------------------

bootstrap_model_metrics <- function(
    pred_df,
    n_boot = 2000,
    seed = 123) {
  
  set.seed(seed)
  
  hc_df <- pred_df %>% filter(Group == "HC")
  ra_df <- pred_df %>% filter(Group == "RA")
  
  boot_results <- vector(
    mode = "list",
    length = n_boot
  )
  
  for (b in seq_len(n_boot)) {
    
    boot_hc <- hc_df[
      sample(
        seq_len(nrow(hc_df)),
        size = nrow(hc_df),
        replace = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    boot_ra <- ra_df[
      sample(
        seq_len(nrow(ra_df)),
        size = nrow(ra_df),
        replace = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    boot_sample <- bind_rows(
      boot_hc,
      boot_ra
    )
    
    boot_results[[b]] <- calculate_model_metrics(
      boot_sample
    )
  }
  
  boot_df <- bind_rows(boot_results)
  
  data.frame(
    Accuracy_Lower =
      quantile(boot_df$Accuracy, 0.025, na.rm = TRUE),
    Accuracy_Upper =
      quantile(boot_df$Accuracy, 0.975, na.rm = TRUE),
    
    Sensitivity_Lower =
      quantile(boot_df$Sensitivity, 0.025, na.rm = TRUE),
    Sensitivity_Upper =
      quantile(boot_df$Sensitivity, 0.975, na.rm = TRUE),
    
    Specificity_Lower =
      quantile(boot_df$Specificity, 0.025, na.rm = TRUE),
    Specificity_Upper =
      quantile(boot_df$Specificity, 0.975, na.rm = TRUE),
    
    ROC_AUC_Lower =
      quantile(boot_df$ROC_AUC, 0.025, na.rm = TRUE),
    ROC_AUC_Upper =
      quantile(boot_df$ROC_AUC, 0.975, na.rm = TRUE)
  )
}

summarize_model_with_ci <- function(
    pred_df,
    model_name,
    seed) {
  
  point_estimates <- calculate_model_metrics(pred_df)
  
  confidence_intervals <- bootstrap_model_metrics(
    pred_df,
    n_boot = 2000,
    seed = seed
  )
  
  bind_cols(
    data.frame(Model = model_name),
    point_estimates,
    confidence_intervals
  )
}

results_lasso_ci <- summarize_model_with_ci(
  pred_lasso,
  "LASSO Logistic",
  seed = 101
)

results_rf_ci <- summarize_model_with_ci(
  pred_rf,
  "Random Forest",
  seed = 202
)

results_svm_ci <- summarize_model_with_ci(
  pred_svm,
  "SVM-RBF",
  seed = 303
)

results_xgb_ci <- summarize_model_with_ci(
  pred_xgb,
  "XGBoost",
  seed = 404
)

model_performance_ci <- bind_rows(
  results_lasso_ci,
  results_rf_ci,
  results_svm_ci,
  results_xgb_ci
) %>%
  arrange(desc(ROC_AUC))

print(model_performance_ci)

write.csv(
  model_performance_ci,
  file.path(
    save_dir,
    "Table_Model_Performance_With_95CI.csv"
  ),
  row.names = FALSE
)

# A formatted table for direct manuscript use
model_performance_formatted <- model_performance_ci %>%
  mutate(
    Accuracy_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      Accuracy,
      Accuracy_Lower,
      Accuracy_Upper
    ),
    Sensitivity_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      Sensitivity,
      Sensitivity_Lower,
      Sensitivity_Upper
    ),
    Specificity_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      Specificity,
      Specificity_Lower,
      Specificity_Upper
    ),
    ROC_AUC_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      ROC_AUC,
      ROC_AUC_Lower,
      ROC_AUC_Upper
    )
  ) %>%
  select(
    Model,
    Accuracy_95CI,
    Sensitivity_95CI,
    Specificity_95CI,
    ROC_AUC_95CI
  )

print(model_performance_formatted)

write.csv(
  model_performance_formatted,
  file.path(
    save_dir,
    "Table_Model_Performance_Formatted_For_Manuscript.csv"
  ),
  row.names = FALSE
)

# Retain baseline_table name so later code does not break
baseline_table <- model_performance_ci %>%
  select(
    Model,
    Accuracy,
    Sensitivity,
    Specificity,
    ROC_AUC
  )

write.csv(
  baseline_table,
  file.path(
    save_dir,
    "Table_Baseline_Model_Comparison_TestSet.csv"
  ),
  row.names = FALSE
)



# -----------------------------------------------------------------------------
# ROC curves for all models
# -----------------------------------------------------------------------------
get_roc_df <- function(final_fit, test_data, model_name) {
  preds <- test_data %>%
    select(Sample, Group) %>%
    bind_cols(
      predict(final_fit, test_data, type = "prob"),
      predict(final_fit, test_data, type = "class")
    ) %>%
    mutate(Group = factor(Group, levels = c("HC", "RA")))
  
  roc_curve(preds, truth = Group, .pred_RA, event_level = "second") %>%
    mutate(Model = model_name)
}

roc_lasso <- get_roc_df(final_fit_lasso, test_data, "LASSO Logistic")
roc_rf    <- get_roc_df(final_fit_rf, test_data, "Random Forest")
roc_svm   <- get_roc_df(final_fit_svm, test_data, "SVM-RBF")
roc_xgb   <- get_roc_df(final_fit_xgb, test_data, "XGBoost")

roc_all <- bind_rows(roc_lasso, roc_rf, roc_svm, roc_xgb)

p_roc_all <- ggplot(roc_all, aes(x = 1 - specificity, y = sensitivity, color = Model)) +
  geom_path(linewidth = 1.1) +
  geom_abline(linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "ROC Curves for ASV-level RA Classification Models",
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

ggsave(
  filename = file.path(save_dir, "Figure_ROC_Model_Comparison.png"),
  plot = p_roc_all, width = 8, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# Keep XGBoost as the primary final model for interpretation
# -----------------------------------------------------------------------------
final_fit_xgb_asv <- final_fit_xgb
best_xgb_asv <- best_xgb
tuned_xgb_asv <- tuned_xgb


# =============================================================================
# 11) XGBoost feature importance
# =============================================================================
xgb_fit_obj <- extract_fit_parsnip(final_fit_xgb_asv)$fit
imp <- xgboost::xgb.importance(model = xgb_fit_obj)

write.csv(
  imp,
  file.path(save_dir, "XGB_ASV_Feature_Importance.csv"),
  row.names = FALSE
)

top_n_imp <- 20
imp_top <- imp %>%
  slice_max(order_by = Gain, n = top_n_imp) %>%
  left_join(tax_info_df, by = c("Feature" = "ASV")) %>%
  mutate(
    Genus = str_extract(Taxon, "g__[^;]+"),
    Genus = ifelse(is.na(Genus), "Unknown", str_remove(Genus, "g__")),
    short_id = str_sub(Feature, 1, 8),
    plot_label = paste0(Genus, " (", short_id, ")")
  )

p_imp <- ggplot(imp_top, aes(x = reorder(plot_label, Gain), y = Gain)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0("Top ", top_n_imp, " ASVs by XGBoost Importance"),
    x = NULL,
    y = "Gain"
  )

ggsave(
  filename = file.path(save_dir, "Figure_XGBoost_Top20_ASV_Importance.png"),
  plot = p_imp, width = 10, height = 7, dpi = 300
)

top_asv_tax <- tibble(ASV = imp_top$Feature) %>%
  left_join(tax_info_df, by = "ASV") %>%
  mutate(
    Genus = str_extract(Taxon, "g__[^;]+"),
    Species = str_extract(Taxon, "s__[^;]+")
  )

write.csv(
  top_asv_tax,
  file.path(save_dir, "Table_XGB_Top_ASVs_With_Taxonomy.csv"),
  row.names = FALSE
)

# =============================================================================
# 11A) Merge XGBoost importance with DESeq2 results for biological interpretation
# =============================================================================
xgb_deseq2_merged <- imp %>%
  left_join(
    final_deseq2 %>%
      select(
        ASV, log2FoldChange, padj, Direction, Genus, Species, Taxon
      ),
    by = c("Feature" = "ASV")
  ) %>%
  arrange(desc(Gain))

write.csv(
  xgb_deseq2_merged,
  file.path(save_dir, "Table_XGB_Importance_Merged_with_DESeq2.csv"),
  row.names = FALSE
)

xgb_deseq2_top20 <- xgb_deseq2_merged %>%
  slice_max(order_by = Gain, n = 20)

write.csv(
  xgb_deseq2_top20,
  file.path(save_dir, "Table_XGB_Top20_Importance_Merged_with_DESeq2.csv"),
  row.names = FALSE
)

# =============================================================================
# 12) Final console summary
# =============================================================================

cat("\n================ FINAL SUMMARY ================\n")

cat("HC samples:", length(hc_cols), "\n")
cat("RA samples:", length(ra_cols), "\n")
cat("Total samples:", nrow(asv_mat), "\n\n")

cat("Original ASVs:", ncol(asv_mat), "\n")
cat("ASVs retained after 5% prevalence filtering:",
    ncol(asv_mat_filt), "\n")

cat("Zero percentage before filtering:",
    sparsity_before$Zero_Percentage, "%\n")

cat("Zero percentage after filtering:",
    sparsity_after$Zero_Percentage, "%\n\n")

cat("Primary ecological resolution: ASV level\n")
cat("Primary internal ML resolution: ASV level\n")
cat("Prevalence threshold: ASV present in at least 5% of samples\n\n")

cat("Results saved in:\n")
cat(save_dir, "\n")

cat("===============================================\n")
# =============================================================================
# 13) Resume repeated validation from completed repetitions
# =============================================================================

# =============================================================================
# Helper function for repeated validation
# =============================================================================

evaluate_repeated_model <- function(
    finalized_workflow,
    train_rep,
    test_rep,
    model_name,
    seed_value) {
  
  fitted_model <- fit(
    finalized_workflow,
    data = train_rep
  )
  
  pred_rep <- test_rep %>%
    dplyr::select(Sample, Group) %>%
    dplyr::bind_cols(
      predict(
        fitted_model,
        new_data = test_rep,
        type = "prob"
      ),
      predict(
        fitted_model,
        new_data = test_rep,
        type = "class"
      )
    ) %>%
    dplyr::mutate(
      Group = factor(Group, levels = c("HC", "RA"))
    )
  
  data.frame(
    Seed = seed_value,
    Model = model_name,
    
    Accuracy = yardstick::accuracy(
      pred_rep,
      truth = Group,
      estimate = .pred_class
    )$.estimate,
    
    Sensitivity = yardstick::sens(
      pred_rep,
      truth = Group,
      estimate = .pred_class,
      event_level = "second"
    )$.estimate,
    
    Specificity = yardstick::spec(
      pred_rep,
      truth = Group,
      estimate = .pred_class,
      event_level = "second"
    )$.estimate,
    
    ROC_AUC = yardstick::roc_auc(
      pred_rep,
      truth = Group,
      .pred_RA,
      event_level = "second"
    )$.estimate
  )
}


N_REPEATS <- 30

set.seed(123)

repeat_seeds <- sample(
  1000:999999,
  size = N_REPEATS,
  replace = FALSE
)

progress_file <- file.path(
  save_dir,
  "Repeated_Validation_All_Models_Progress.csv"
)

# Load completed results if the progress file exists
if (file.exists(progress_file)) {
  
  progress_results <- read.csv(
    progress_file,
    stringsAsFactors = FALSE
  )
  
  completed_seeds <- unique(progress_results$Seed)
  
  cat(
    "Found",
    length(completed_seeds),
    "completed repetitions.\n"
  )
  
} else {
  
  progress_results <- data.frame()
  completed_seeds <- numeric(0)
  
  cat("No previous completed repetitions found.\n")
}

# Identify seeds still needing analysis
remaining_seeds <- repeat_seeds[
  !repeat_seeds %in% completed_seeds
]

cat(
  "Remaining repetitions:",
  length(remaining_seeds),
  "\n"
)

for (current_seed in remaining_seeds) {
  
  repetition_number <- match(
    current_seed,
    repeat_seeds
  )
  
  cat(
    "\nStarting repeated validation",
    repetition_number,
    "of",
    N_REPEATS,
    "- seed:",
    current_seed,
    "\n"
  )
  
  set.seed(current_seed)
  
  split_rep <- rsample::initial_split(
    ml_df,
    prop = 0.8,
    strata = Group
  )
  
  train_rep <- rsample::training(split_rep)
  test_rep  <- rsample::testing(split_rep)
  
  # -------------------------------------------------------------------------
  # LASSO
  # -------------------------------------------------------------------------
  
  result_lasso <- evaluate_repeated_model(
    finalized_workflow = final_wf_lasso,
    train_rep = train_rep,
    test_rep = test_rep,
    model_name = "LASSO Logistic",
    seed_value = current_seed
  )
  
  cat("  LASSO completed\n")
  
  # -------------------------------------------------------------------------
  # Random Forest
  # -------------------------------------------------------------------------
  
  result_rf <- evaluate_repeated_model(
    finalized_workflow = final_wf_rf,
    train_rep = train_rep,
    test_rep = test_rep,
    model_name = "Random Forest",
    seed_value = current_seed
  )
  
  cat("  Random Forest completed\n")
  
  # -------------------------------------------------------------------------
  # SVM-RBF
  # -------------------------------------------------------------------------
  
  result_svm <- evaluate_repeated_model(
    finalized_workflow = final_wf_svm,
    train_rep = train_rep,
    test_rep = test_rep,
    model_name = "SVM-RBF",
    seed_value = current_seed
  )
  
  cat("  SVM-RBF completed\n")
  
  # -------------------------------------------------------------------------
  # XGBoost
  # -------------------------------------------------------------------------
  
  result_xgb <- evaluate_repeated_model(
    finalized_workflow = final_wf_xgb,
    train_rep = train_rep,
    test_rep = test_rep,
    model_name = "XGBoost",
    seed_value = current_seed
  )
  
  cat("  XGBoost completed\n")
  
  current_results <- dplyr::bind_rows(
    result_lasso,
    result_rf,
    result_svm,
    result_xgb
  )
  
  # Append current repetition to previously completed results
  progress_results <- dplyr::bind_rows(
    progress_results,
    current_results
  ) %>%
    dplyr::distinct(
      Seed,
      Model,
      .keep_all = TRUE
    )
  
  # Save only after all four models complete
  write.csv(
    progress_results,
    progress_file,
    row.names = FALSE
  )
  
  cat(
    "  Repetition",
    repetition_number,
    "saved successfully.\n"
  )
  
  # Clean temporary fitted objects between repetitions
  rm(
    result_lasso,
    result_rf,
    result_svm,
    result_xgb,
    current_results,
    train_rep,
    test_rep,
    split_rep
  )
  
  invisible(gc())
}

repeated_results_all_models <- progress_results %>%
  dplyr::arrange(
    match(Seed, repeat_seeds),
    Model
  )

write.csv(
  repeated_results_all_models,
  file.path(
    save_dir,
    "Repeated_Validation_All_Models_All_Runs.csv"
  ),
  row.names = FALSE
)

print(repeated_results_all_models)

# Keep only seeds from the current 30-repeat run
repeated_results_all_models <- repeated_results_all_models %>%
  dplyr::filter(Seed %in% repeat_seeds) %>%
  dplyr::distinct(
    Seed,
    Model,
    .keep_all = TRUE
  ) %>%
  dplyr::arrange(
    match(Seed, repeat_seeds),
    factor(
      Model,
      levels = c(
        "LASSO Logistic",
        "Random Forest",
        "SVM-RBF",
        "XGBoost"
      )
    )
  )

# Verify expected dimensions
cat(
  "Rows:",
  nrow(repeated_results_all_models),
  "\n"
)

cat(
  "Unique seeds:",
  dplyr::n_distinct(repeated_results_all_models$Seed),
  "\n"
)

print(
  table(repeated_results_all_models$Model)
)

write.csv(
  repeated_results_all_models,
  file.path(
    save_dir,
    "Repeated_Validation_All_Models_All_Runs.csv"
  ),
  row.names = FALSE
)

write.csv(
  repeated_results_all_models,
  file.path(
    save_dir,
    "Repeated_Validation_All_Models_Progress.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Summary of repeated validation
# =============================================================================

repeated_results_long <- repeated_results_all_models %>%
  pivot_longer(
    cols = c(
      Accuracy,
      Sensitivity,
      Specificity,
      ROC_AUC
    ),
    names_to = "Metric",
    values_to = "Value"
  )

repeated_summary_all_models <- repeated_results_long %>%
  group_by(Model, Metric) %>%
  summarise(
    Repetitions = sum(!is.na(Value)),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    Minimum = min(Value, na.rm = TRUE),
    Maximum = max(Value, na.rm = TRUE),
    
    # 95% CI for mean performance across repeated splits
    CI_Lower = Mean -
      qt(0.975, df = Repetitions - 1) *
      SD / sqrt(Repetitions),
    
    CI_Upper = Mean +
      qt(0.975, df = Repetitions - 1) *
      SD / sqrt(Repetitions),
    
    .groups = "drop"
  )

print(repeated_summary_all_models)

write.csv(
  repeated_summary_all_models,
  file.path(
    save_dir,
    "Repeated_Validation_All_Models_Summary.csv"
  ),
  row.names = FALSE
)

# Formatted summary table
repeated_summary_formatted <- repeated_summary_all_models %>%
  mutate(
    Mean_SD = sprintf(
      "%.3f ± %.3f",
      Mean,
      SD
    ),
    Mean_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      Mean,
      CI_Lower,
      CI_Upper
    )
  ) %>%
  select(
    Model,
    Metric,
    Repetitions,
    Mean_SD,
    Mean_95CI,
    Median,
    Q1,
    Q3,
    Minimum,
    Maximum
  )

write.csv(
  repeated_summary_formatted,
  file.path(
    save_dir,
    "Repeated_Validation_Formatted_For_Manuscript.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# Figure: repeated performance for all models
# =============================================================================

p_repeat_all_models <- ggplot(
  repeated_results_long,
  aes(
    x = Model,
    y = Value,
    fill = Model
  )
) +
  geom_boxplot(
    outlier.alpha = 0.5
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0(
      "Performance Across ",
      N_REPEATS,
      " Repeated Stratified Train-Test Splits"
    ),
    x = "Machine-learning model",
    y = "Performance"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    save_dir,
    "Figure_Repeated_Validation_All_Models.png"
  ),
  plot = p_repeat_all_models,
  width = 12,
  height = 8,
  dpi = 300
)

# =============================================================================
# Repeated ROC-AUC ranking
# =============================================================================

repeated_auc_ranking <- repeated_summary_all_models %>%
  filter(Metric == "ROC_AUC") %>%
  arrange(desc(Mean))

print(repeated_auc_ranking)

write.csv(
  repeated_auc_ranking,
  file.path(
    save_dir,
    "Table_Repeated_ROC_AUC_Model_Ranking.csv"
  ),
  row.names = FALSE
)
