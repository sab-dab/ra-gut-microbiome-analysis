

save_dir <- file.path(
  getwd(),
  "results_final"
)

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# =============================================================================
# 1) Packages
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
  BiocManager::install("ANCOMBC", update = FALSE, ask = FALSE)
}


BiocManager::install("phyloseq", ask = FALSE, update = FALSE)

library(ANCOMBC)
library(dplyr)
library(tibble)
library(stringr)
library(phyloseq)

# =============================================================================
# 2) Load discovery data
# =============================================================================

asv_profile <- readRDS("1.ASV.profile.rds")
tax_info    <- readRDS("1.taxonomy.info.rds")

# Samples as rows, ASVs as columns
asv_counts <- t(as.matrix(asv_profile))

metadata <- data.frame(
  Sample = rownames(asv_counts),
  Group = ifelse(
    grepl("^HC", rownames(asv_counts)),
    "HC",
    "RA"
  ),
  row.names = rownames(asv_counts)
)

metadata$Group <- factor(
  metadata$Group,
  levels = c("HC", "RA")
)

cat("Samples:", nrow(asv_counts), "\n")
cat("ASVs:", ncol(asv_counts), "\n")
cat("HC:", sum(metadata$Group == "HC"), "\n")
cat("RA:", sum(metadata$Group == "RA"), "\n")

# =============================================================================
# 3) Aggregate raw ASV counts to genus level
# =============================================================================

tax_df <- as.data.frame(tax_info, stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("ASV")

tax_df$Genus <- stringr::str_extract(
  tax_df$Taxon,
  "g__[^;]+"
)

tax_df$Genus <- ifelse(
  is.na(tax_df$Genus),
  "Unknown",
  stringr::str_remove(tax_df$Genus, "^g__")
)

# Match taxonomy rows to ASV-count columns
tax_df <- tax_df[
  match(colnames(asv_counts), tax_df$ASV),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    tax_df$ASV,
    colnames(asv_counts)
  )
)

# Sum raw ASV counts belonging to the same genus
genus_counts <- rowsum(
  t(asv_counts),
  group = tax_df$Genus,
  reorder = FALSE
)

# Convert back to samples × genera
genus_counts <- t(genus_counts)

# Remove unclassified genus before DA analysis
if ("Unknown" %in% colnames(genus_counts)) {
  genus_counts <- genus_counts[
    ,
    colnames(genus_counts) != "Unknown",
    drop = FALSE
  ]
}

cat("Samples in genus table:", nrow(genus_counts), "\n")
cat("Genera before prevalence filtering:", ncol(genus_counts), "\n")

# =============================================================================
# 4) Genus prevalence filtering
# Retain genera present in at least 5% of samples
# =============================================================================

genus_prevalence <- colSums(genus_counts > 0)

keep_genera <- genus_prevalence >=
  ceiling(0.05 * nrow(genus_counts))

genus_counts_filt <- genus_counts[
  ,
  keep_genera,
  drop = FALSE
]

cat(
  "Genera after 5% prevalence filtering:",
  ncol(genus_counts_filt),
  "\n"
)

write.csv(
  data.frame(
    Genus = colnames(genus_counts),
    Samples_Present = genus_prevalence,
    Prevalence_Percentage = round(
      100 * genus_prevalence / nrow(genus_counts),
      3
    ),
    Retained = keep_genera
  ),
  file.path(
    save_dir,
    "Table_ANCOMBC2_Genus_Prevalence_Filter.csv"
  ),
  row.names = FALSE
)

cat("Genera before prevalence filtering:", ncol(genus_counts), "\n")
cat("Genera after 5% prevalence filtering:", ncol(genus_counts_filt), "\n")

# =============================================================================
# 5) Run ANCOM-BC2
# =============================================================================

ancom_input <- t(genus_counts_filt)

ancom_res <- ANCOMBC::ancombc2(
  data = ancom_input,
  meta_data = metadata,
  fix_formula = "Group",
  p_adj_method = "BH",
  prv_cut = 0,
  lib_cut = 0,
  group = "Group",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1
)

names(ancom_res)
ancom_primary <- ancom_res$res

str(ancom_primary)
head(ancom_primary)
colnames(ancom_primary)

ancombc2_results <- ancom_primary %>%
  dplyr::select(
    taxon,
    lfc_GroupRA,
    se_GroupRA,
    W_GroupRA,
    p_GroupRA,
    q_GroupRA,
    diff_GroupRA,
    passed_ss_GroupRA,
    diff_robust_GroupRA
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      lfc_GroupRA > 0 ~ "Higher in RA",
      lfc_GroupRA < 0 ~ "Higher in HC",
      TRUE ~ "No difference"
    )
  ) %>%
  dplyr::arrange(q_GroupRA)

write.csv(
  ancombc2_results,
  file.path(
    save_dir,
    "ANCOMBC2_Genus_RA_vs_HC_All.csv"
  ),
  row.names = FALSE
)

ancombc2_significant <- ancombc2_results %>%
  dplyr::filter(
    diff_robust_GroupRA == TRUE
  )

print(ancombc2_significant)

write.csv(
  ancombc2_significant,
  file.path(
    save_dir,
    "ANCOMBC2_Genus_RA_vs_HC_Robust_Significant.csv"
  ),
  row.names = FALSE
)
# =============================================================================
# Figure 3: Robust genus-level differential abundance by ANCOM-BC2
# =============================================================================

library(ggplot2)
library(dplyr)
library(stringr)

figure3_data <- ancombc2_significant %>%
  dplyr::mutate(
    Genus = taxon,
    Direction = factor(
      Direction,
      levels = c("Higher in HC", "Higher in RA")
    ),
    Genus = stringr::str_replace_all(Genus, "_", " ")
  ) %>%
  dplyr::arrange(lfc_GroupRA) %>%
  dplyr::mutate(
    Genus = factor(Genus, levels = Genus)
  )

p_figure3 <- ggplot(
  figure3_data,
  aes(
    x = lfc_GroupRA,
    y = Genus,
    fill = Direction
  )
) +
  geom_col(width = 0.75) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Higher in HC" = "#3B82B8",
      "Higher in RA" = "#D96B5F"
    )
  ) +
  labs(
    title = "Robust Genus-Level Differential Abundance by ANCOM-BC2",
    subtitle = "Positive coefficients indicate higher abundance in RA; negative coefficients indicate higher abundance in HC",
    x = "Bias-corrected coefficient for RA versus HC",
    y = "Genus",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text.y = element_text(size = 8),
    legend.position = "top"
  )

print(p_figure3)

ggsave(
  filename = file.path(
    save_dir,
    "Figure3_ANCOMBC2_Robust_Differential_Abundance.png"
  ),
  plot = p_figure3,
  width = 9,
  height = 10,
  dpi = 600
)

ggsave(
  filename = file.path(
    save_dir,
    "Figure3_ANCOMBC2_Robust_Differential_Abundance.pdf"
  ),
  plot = p_figure3,
  width = 9,
  height = 10
)
