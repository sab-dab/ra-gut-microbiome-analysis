library(dplyr)
library(tibble)
library(stringr)
library(DESeq2)
library(ggplot2)

project_dir <- getwd()

output_path <- file.path(
  project_dir,
  "external_validation_results"
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

# -----------------------------------------------------------------------------
# 1. Load external genus count table
# -----------------------------------------------------------------------------

external_genus_counts <- read.csv(
  file.path(
    output_path,
    "External_Genus_Count_Table_Combined.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Rows = samples; columns = genera
rownames(external_genus_counts) <- external_genus_counts$Sample
external_genus_counts$Sample <- NULL

cat(
  "External genus table:",
  nrow(external_genus_counts),
  "samples x",
  ncol(external_genus_counts),
  "genera\n"
)

# -----------------------------------------------------------------------------
# 2. Create RA/HC metadata
# -----------------------------------------------------------------------------

metadata_all <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

external_meta <- metadata_all %>%
  transmute(
    Sample = Run,
    Group = case_when(
      grepl("^RA_", `Sample Name`) ~ "RA",
      grepl("^GUT_", `Sample Name`) ~ "HC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    Sample %in% rownames(external_genus_counts),
    !is.na(Group)
  )

external_meta <- external_meta %>%
  distinct(Sample, .keep_all = TRUE) %>%
  arrange(match(Sample, rownames(external_genus_counts)))

external_meta$Group <- factor(
  external_meta$Group,
  levels = c("HC", "RA")
)

rownames(external_meta) <- external_meta$Sample

# Match count table to metadata order
external_genus_counts <- external_genus_counts[
  external_meta$Sample,
  ,
  drop = FALSE
]

cat(
  "RA samples:",
  sum(external_meta$Group == "RA"),
  "\nHC samples:",
  sum(external_meta$Group == "HC"),
  "\n"
)

stopifnot(
  identical(
    rownames(external_genus_counts),
    rownames(external_meta)
  )
)

# -----------------------------------------------------------------------------
# 3. Filter very rare genera
# Retain genera present in at least 5% of external samples
# -----------------------------------------------------------------------------

minimum_samples <- ceiling(
  0.05 * nrow(external_genus_counts)
)

keep_genera <- colSums(
  external_genus_counts > 0,
  na.rm = TRUE
) >= minimum_samples

external_genus_filtered <- external_genus_counts[
  ,
  keep_genera,
  drop = FALSE
]

cat(
  "Genera tested after 5% prevalence filtering:",
  ncol(external_genus_filtered),
  "\n"
)

# -----------------------------------------------------------------------------
# 4. Prepare genus-by-sample count matrix for DESeq2
# -----------------------------------------------------------------------------

external_count_mat <- t(
  as.matrix(external_genus_filtered)
)

storage.mode(external_count_mat) <- "integer"

# -----------------------------------------------------------------------------
# 5. Run genus-level DESeq2 analysis
# -----------------------------------------------------------------------------

dds_external <- DESeqDataSetFromMatrix(
  countData = external_count_mat,
  colData = external_meta,
  design = ~ Group
)

dds_external <- estimateSizeFactors(
  dds_external,
  type = "poscounts"
)

dds_external <- DESeq(
  dds_external
)

external_res <- results(
  dds_external,
  contrast = c("Group", "RA", "HC")
)

external_da <- as.data.frame(external_res) %>%
  rownames_to_column("Genus") %>%
  mutate(
    Direction = case_when(
      log2FoldChange > 0 ~ "Higher in RA",
      log2FoldChange < 0 ~ "Higher in HC",
      TRUE ~ "No change"
    )
  ) %>%
  arrange(padj)

# -----------------------------------------------------------------------------
# 6. Add abundance and prevalence summaries
# -----------------------------------------------------------------------------

hc_samples <- external_meta$Sample[
  external_meta$Group == "HC"
]

ra_samples <- external_meta$Sample[
  external_meta$Group == "RA"
]

external_group_summary <- tibble(
  Genus = colnames(external_genus_filtered),
  
  Mean_Abundance_HC = colMeans(
    external_genus_filtered[
      hc_samples,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  ),
  
  Mean_Abundance_RA = colMeans(
    external_genus_filtered[
      ra_samples,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  ),
  
  Prevalence_HC_Percent = colMeans(
    external_genus_filtered[
      hc_samples,
      ,
      drop = FALSE
    ] > 0,
    na.rm = TRUE
  ) * 100,
  
  Prevalence_RA_Percent = colMeans(
    external_genus_filtered[
      ra_samples,
      ,
      drop = FALSE
    ] > 0,
    na.rm = TRUE
  ) * 100
)

external_da_complete <- external_da %>%
  left_join(
    external_group_summary,
    by = "Genus"
  )

external_da_significant <- external_da_complete %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

cat(
  "FDR-significant genera:",
  nrow(external_da_significant),
  "\n"
)

print(external_da_significant)


# =============================================================================
# Summarize significant external differential-abundance findings
# =============================================================================

n_external_sig <- nrow(external_da_significant)

n_external_ra <- sum(
  external_da_significant$Direction == "Higher in RA"
)

n_external_hc <- sum(
  external_da_significant$Direction == "Higher in HC"
)

cat(
  "\nExternal DA summary:\n",
  "Total FDR-significant genera:", n_external_sig, "\n",
  "Higher in RA:", n_external_ra, "\n",
  "Higher in HC:", n_external_hc, "\n"
)
# -----------------------------------------------------------------------------
# 7. Save complete and significant results
# -----------------------------------------------------------------------------

write.csv(
  external_da_complete,
  file.path(
    output_path,
    "External_Genus_DESeq2_All_Results.csv"
  ),
  row.names = FALSE
)

write.csv(
  external_da_significant,
  file.path(
    output_path,
    "External_Genus_DESeq2_Significant_FDR0.05.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# External differential abundance figure
# Only FDR-significant genera are displayed
# =============================================================================

external_da_plot <- external_da_significant %>%
  mutate(
    Genus_Label = stringr::str_replace_all(
      Genus,
      "_",
      " "
    ),
    Lower_CI = log2FoldChange - 1.96 * lfcSE,
    Upper_CI = log2FoldChange + 1.96 * lfcSE,
    FDR_Label = paste0(
      "FDR = ",
      format(
        padj,
        scientific = TRUE,
        digits = 2
      )
    )
  )

p_external_da <- ggplot(
  external_da_plot,
  aes(
    x = log2FoldChange,
    y = Genus_Label
  )
) +
  geom_segment(
    aes(
      x = Lower_CI,
      xend = Upper_CI,
      y = Genus_Label,
      yend = Genus_Label
    ),
    linewidth = 1
  ) +
  geom_point(
    size = 4,
    colour = "#0072B2"
  ) +
  geom_text(
    aes(
      x = log2FoldChange,
      label = FDR_Label
    ),
    nudge_y = 0.12,
    size = 4,
    fontface = "bold"
  ) +
  scale_x_continuous(
    limits = c(
      floor(min(external_da_plot$Lower_CI)) - 1,
      ceiling(max(external_da_plot$Upper_CI)) + 1
    )
  ) +
  labs(
    title = "External Cohort Differential Abundance",
    subtitle = "Exploratory genus-level DESeq2 analysis",
    x = "Log2 fold change (RA versus HC)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 12
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      10,
      20,
      10,
      20
    )
  )

print(p_external_da)

ggsave(
  filename = file.path(
    output_path,
    "Figure_External_Genus_DESeq2_Significant.png"
  ),
  plot = p_external_da,
  width = 7,
  height = 3.5,
  dpi = 600,
  bg = "white"
)

# =============================================================================
# Supplementary Table 
# External genus-level differential abundance results
# =============================================================================

supplementary_external_da <- external_da_complete %>%
  select(
    Genus,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,
    Direction,
    Mean_Abundance_HC,
    Mean_Abundance_RA,
    Prevalence_HC_Percent,
    Prevalence_RA_Percent
  ) %>%
  arrange(padj)

colnames(supplementary_external_da) <- c(
  "Genus",
  "Base mean",
  "Log2 fold change",
  "Standard error",
  "Wald statistic",
  "Raw P value",
  "Adjusted P value (FDR)",
  "Direction",
  "Mean abundance (HC)",
  "Mean abundance (RA)",
  "Prevalence (%) HC",
  "Prevalence (%) RA"
)

write.csv(
  supplementary_external_da,
  file.path(
    output_path,
    "Supplementary_Table_External_Genus_Differential_Abundance.csv"
  ),
  row.names = FALSE
)

print(supplementary_external_da)

# =============================================================================
# Compare external significant genera with discovery ANCOM-BC2 findings
# =============================================================================

discovery_da_file <- file.path(
  "C:/Users/sabir/Downloads/gut microbiota",
  "results_final",
  "ANCOMBC2_Genus_RA_vs_HC_Robust_Significant.csv"
)

discovery_da <- read.csv(
  discovery_da_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Normalize genus names for matching
normalize_genus_name <- function(x) {
  x <- tolower(x)
  x <- gsub("^g__", "", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("[^a-z0-9_]", "", x)
  trimws(x)
}

external_sig_clean <- normalize_genus_name(
  external_da_significant$Genus
)

discovery_sig_clean <- normalize_genus_name(
  discovery_da$taxon
)

overlap_da <- base::intersect(
  external_sig_clean,
  discovery_sig_clean
)

cat(
  "\nCross-cohort differential-abundance overlap:\n",
  "External significant genera:",
  length(external_sig_clean),
  "\n",
  "Discovery robust significant genera:",
  length(discovery_sig_clean),
  "\n",
  "Significant genera shared across cohorts:",
  length(overlap_da),
  "\n"
)

if (length(overlap_da) > 0) {
  cat(
    "Shared significant genera:",
    paste(overlap_da, collapse = ", "),
    "\n"
  )
} else {
  cat(
    "No significant genera were shared across cohorts.\n"
  )
}

# Identify the external significant genus/genus names
cat(
  "\nExternal FDR-significant genus/genera:\n",
  paste(
    external_da_significant$Genus,
    collapse = ", "
  ),
  "\n"
)

# Save summary
cross_cohort_da_summary <- data.frame(
  Metric = c(
    "External FDR-significant genera",
    "Discovery robust ANCOM-BC2 significant genera",
    "Significant genera shared across cohorts"
  ),
  Number = c(
    length(external_sig_clean),
    length(discovery_sig_clean),
    length(overlap_da)
  )
)

write.csv(
  cross_cohort_da_summary,
  file.path(
    output_path,
    "Cross_Cohort_DA_Overlap_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Genus = overlap_da
  ),
  file.path(
    output_path,
    "Cross_Cohort_DA_Shared_Genera.csv"
  ),
  row.names = FALSE
)
