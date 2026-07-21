# ============================================================
# SCRIPT 1: DADA2 PROCESSING ONLY
# ============================================================

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

library(dada2)
library(dplyr)
library(tibble)

fnFs <- sort(list.files(external_path, pattern = "_1.fastq.gz$", full.names = TRUE))
sample_names <- sub("_1.fastq.gz$", "", basename(fnFs))
fnRs <- file.path(external_path, paste0(sample_names, "_2.fastq.gz"))

cat("Forward files:", length(fnFs), "\n")
cat("Reverse files:", length(fnRs), "\n")
cat("Missing reverse files:", sum(!file.exists(fnRs)), "\n")

if (length(fnFs) == 0) stop("No forward FASTQ files found.")
if (any(!file.exists(fnRs))) stop("One or more reverse FASTQ files are missing.")

filtered_path <- file.path(external_path, "filtered")
if (dir.exists(filtered_path)) unlink(filtered_path, recursive = TRUE)
dir.create(filtered_path, recursive = TRUE)

filtFs <- file.path(filtered_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_path, paste0(sample_names, "_R_filt.fastq.gz"))

filter_output <- dada2::filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(280, 250),
  maxN = 0,
  maxEE = c(5, 5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
)

write.csv(as.data.frame(filter_output),
          file.path(output_path, "External_DADA2_Filtering_Input_Output.csv"),
          row.names = TRUE)

errF <- dada2::learnErrors(filtFs, nbases = 1e7, multithread = FALSE)
errR <- dada2::learnErrors(filtRs, nbases = 1e7, multithread = FALSE)
saveRDS(errF, file.path(output_path, "external_errF.rds"))
saveRDS(errR, file.path(output_path, "external_errR.rds"))

derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

dadaFs <- dada2::dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada2::dada(derepRs, err = errR, multithread = FALSE)

mergers <- dada2::mergePairs(
  dadaFs, derepFs, dadaRs, derepRs,
  minOverlap = 8,
  maxMismatch = 2,
  verbose = TRUE
)

seqtab <- dada2::makeSequenceTable(mergers)
saveRDS(seqtab, file.path(output_path, "external_seqtab.rds"))
###############################################################################
#this part was run 
seqtab_nochim <- dada2::removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)
saveRDS(seqtab_nochim, file.path(output_path, "external_seqtab_nochim.rds"))

getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(
  input = filter_output[, 1],
  filtered = filter_output[, 2],
  denoisedF = vapply(dadaFs, getN, numeric(1)),
  denoisedR = vapply(dadaRs, getN, numeric(1)),
  merged = vapply(mergers, getN, numeric(1)),
  nonchim = rowSums(seqtab_nochim)
)
rownames(track) <- sample_names
saveRDS(track, file.path(output_path, "external_read_tracking.rds"))

track_df <- as.data.frame(track) |>
  tibble::rownames_to_column("Sample") |>
  dplyr::mutate(
    Filtered_Percentage = 100 * filtered / input,
    Merged_Percentage = 100 * merged / input,
    Nonchim_Percentage = 100 * nonchim / input
  )
write.csv(track_df,
          file.path(output_path, "Table_External_Read_Retention_By_Sample.csv"),
          row.names = FALSE)

qc_summary <- data.frame(
  Metric = c(
    "FASTQ pairs",
    "Samples in sequence table",
    "Samples with nonzero merged reads",
    "Samples with nonzero nonchimeric reads",
    "ASVs before chimera removal",
    "ASVs after chimera removal",
    "Zero-depth samples after chimera removal"
  ),
  Value = c(
    length(sample_names),
    nrow(seqtab),
    sum(rowSums(seqtab) > 0),
    sum(rowSums(seqtab_nochim) > 0),
    ncol(seqtab),
    ncol(seqtab_nochim),
    sum(rowSums(seqtab_nochim) == 0)
  )
)
write.csv(qc_summary,
          file.path(output_path, "Table_External_DADA2_QC_Summary.csv"),
          row.names = FALSE)

print(track)
print(qc_summary)
cat("\nSCRIPT 1 COMPLETE\n")
cat("Results saved in:", output_path, "\n")

############################################################################


seqtab <- readRDS(
  file.path(output_path, "external_seqtab.rds")
)

dim(seqtab)

nrow(seqtab)
ncol(seqtab)
summary(rowSums(seqtab))
summary(colSums(seqtab))

seq_lengths <- nchar(colnames(seqtab))

summary(seq_lengths)
table(seq_lengths)

asv_totals <- colSums(seqtab)

table(asv_totals <= 1)
table(asv_totals <= 2)

summary(seq_lengths)
sort(table(seq_lengths), decreasing = TRUE)[1:20]

lengths <- nchar(colnames(seqtab))

dominant_length_by_sample <- apply(
  seqtab,
  1,
  function(x) {
    totals_by_length <- tapply(x, lengths, sum)
    as.numeric(names(which.max(totals_by_length)))
  }
)

data.frame(
  Sample = rownames(seqtab),
  Dominant_Length = dominant_length_by_sample,
  Total_Reads = rowSums(seqtab)
)

table(dominant_length_by_sample)

write.csv(
  data.frame(
    Sample = rownames(seqtab),
    Dominant_Length = dominant_length_by_sample,
    Total_Reads = rowSums(seqtab)
  ),
  file.path(
    output_path,
    "External_Dominant_Sequence_Length_By_Sample.csv"
  ),
  row.names = FALSE
)

seqtab_456 <- seqtab[
  dominant_length_by_sample == 456,
  ,
  drop = FALSE
]

seqtab_476 <- seqtab[
  dominant_length_by_sample == 476,
  ,
  drop = FALSE
]

seqtab_481 <- seqtab[
  dominant_length_by_sample == 481,
  ,
  drop = FALSE
]

seqtab_456 <- seqtab_456[
  ,
  nchar(colnames(seqtab_456)) >= 455 &
    nchar(colnames(seqtab_456)) <= 459,
  drop = FALSE
]

seqtab_476 <- seqtab_476[
  ,
  nchar(colnames(seqtab_476)) >= 474 &
    nchar(colnames(seqtab_476)) <= 478,
  drop = FALSE
]

seqtab_481 <- seqtab_481[
  ,
  nchar(colnames(seqtab_481)) >= 480 &
    nchar(colnames(seqtab_481)) <= 482,
  drop = FALSE
]

keep_asv_456 <- colSums(seqtab_456) >= 10

seqtab_456_filt <- seqtab_456[
  ,
  keep_asv_456,
  drop = FALSE
]

dim(seqtab_456_filt)
nochim_456 <- dada2::removeBimeraDenovo(
  seqtab_456_filt,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

saveRDS(
  nochim_456,
  file.path(
    output_path,
    "external_seqtab_456_nochim.rds"
  )
)

qc_456 <- data.frame(
  Samples = nrow(nochim_456),
  ASVs_before_chimera = ncol(seqtab_456_filt),
  ASVs_after_chimera = ncol(nochim_456),
  Chimeras_removed = ncol(seqtab_456_filt) - ncol(nochim_456),
  Samples_with_reads = sum(rowSums(nochim_456) > 0),
  Zero_depth_samples = sum(rowSums(nochim_456) == 0)
)

print(qc_456)

write.csv(
  qc_456,
  file.path(output_path, "Table_External_456bp_QC_Summary.csv"),
  row.names = FALSE
)

keep_asv_476 <- colSums(seqtab_476) >= 10

seqtab_476_filt <- seqtab_476[
  ,
  keep_asv_476,
  drop = FALSE
]

dim(seqtab_476_filt)

nochim_476 <- dada2::removeBimeraDenovo(
  seqtab_476_filt,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

saveRDS(
  nochim_476,
  file.path(
    output_path,
    "external_seqtab_476_nochim.rds"
  )
)

keep_asv_481 <- colSums(seqtab_481) >= 10

seqtab_481_filt <- seqtab_481[
  ,
  keep_asv_481,
  drop = FALSE
]

dim(seqtab_481_filt)

nochim_481 <- dada2::removeBimeraDenovo(
  seqtab_481_filt,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

saveRDS(
  nochim_481,
  file.path(
    output_path,
    "external_seqtab_481_nochim.rds"
  )
)

dim(nochim_456)
dim(nochim_476)
dim(nochim_481)

sum(rowSums(nochim_456) > 0)
sum(rowSums(nochim_476) > 0)
sum(rowSums(nochim_481) > 0)


length_group_qc <- data.frame(
  Length_Group = c("456 bp", "476 bp", "481 bp"),
  Samples = c(
    nrow(nochim_456),
    nrow(nochim_476),
    nrow(nochim_481)
  ),
  ASVs_After_Chimera_Removal = c(
    ncol(nochim_456),
    ncol(nochim_476),
    ncol(nochim_481)
  ),
  Samples_With_Reads = c(
    sum(rowSums(nochim_456) > 0),
    sum(rowSums(nochim_476) > 0),
    sum(rowSums(nochim_481) > 0)
  )
)

print(length_group_qc)

write.csv(
  length_group_qc,
  file.path(output_path, "Table_External_Length_Group_QC.csv"),
  row.names = FALSE
)

#assign taxonomy to each length group
silva_path <- "C:/silva_nr99_v138.1_train_set.fa.gz"

taxa_456 <- dada2::assignTaxonomy(
  nochim_456,
  silva_path,
  multithread = FALSE,
  tryRC = TRUE
)

taxa_476 <- dada2::assignTaxonomy(
  nochim_476,
  silva_path,
  multithread = FALSE,
  tryRC = TRUE
)

taxa_481 <- dada2::assignTaxonomy(
  nochim_481,
  silva_path,
  multithread = FALSE,
  tryRC = TRUE
) 

saveRDS(
  taxa_456,
  file.path(output_path, "external_taxonomy_456.rds")
)

saveRDS(
  taxa_476,
  file.path(output_path, "external_taxonomy_476.rds")
)

saveRDS(
  taxa_481,
  file.path(output_path, "external_taxonomy_481.rds")
) 

dim(taxa_456)
dim(taxa_476)
dim(taxa_481)
#############################################################################
#aggregate ASV to genus
normalize_genus <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unknown"
  x <- trimws(x)
  x <- sub("^g__", "", x)
  x <- gsub("[\\[\\]]", "", x)
  x <- gsub("[[:space:]]+", "_", x)
  x
}

aggregate_to_genus <- function(seqtab_group, taxa_group) {
  taxa_df <- as.data.frame(taxa_group, stringsAsFactors = FALSE)
  taxa_df <- taxa_df[colnames(seqtab_group), , drop = FALSE]
  
  genus_names <- normalize_genus(taxa_df$Genus)
  
  genus_mat <- rowsum(
    t(seqtab_group),
    group = genus_names,
    reorder = FALSE
  )
  
  as.data.frame(
    t(genus_mat),
    check.names = FALSE
  )
}

genus_456 <- aggregate_to_genus(nochim_456, taxa_456)
genus_476 <- aggregate_to_genus(nochim_476, taxa_476)
genus_481 <- aggregate_to_genus(nochim_481, taxa_481)

#combine
all_genera <- Reduce(
  union,
  list(
    colnames(genus_456),
    colnames(genus_476),
    colnames(genus_481)
  )
)

add_missing_genera <- function(df, all_genera) {
  missing <- base::setdiff(all_genera, colnames(df))
  
  for (g in missing) {
    df[[g]] <- 0
  }
  
  df[, all_genera, drop = FALSE]
}

genus_456_aligned <- add_missing_genera(genus_456, all_genera)
genus_476_aligned <- add_missing_genera(genus_476, all_genera)
genus_481_aligned <- add_missing_genera(genus_481, all_genera)

external_genus_counts <- rbind(
  genus_456_aligned,
  genus_476_aligned,
  genus_481_aligned
)

dim(external_genus_counts)
sum(rowSums(external_genus_counts) > 0)

saveRDS(
  external_genus_counts,
  file.path(
    output_path,
    "external_genus_counts_combined.rds"
  )
)

write.csv(
  cbind(
    Sample = rownames(external_genus_counts),
    external_genus_counts
  ),
  file.path(
    output_path,
    "External_Genus_Count_Table_Combined.csv"
  ),
  row.names = FALSE
)


discovery_features <- readRDS(
  file.path(
    output_path,
    "discovery_genus_training_features.rds"
  )
)

shared_genera <- intersect(
  discovery_features,
  colnames(external_genus_counts)
)

missing_genera <- base::setdiff(
  discovery_features,
  colnames(external_genus_counts)
)

cat("Discovery genera:", length(discovery_features), "\n")
cat("External genera:", ncol(external_genus_counts), "\n")
cat("Shared genera:", length(shared_genera), "\n")
cat("Missing discovery genera:", length(missing_genera), "\n")
cat(
  "Overlap percentage:",
  round(
    100 * length(shared_genera) / length(discovery_features),
    2
  ),
  "%\n"
)

head(discovery_features, 20)
head(colnames(external_genus_counts), 20)

normalize_genus_name <- function(x) {
  x <- as.character(x)
  x <- sub("^g__", "", x)
  x <- gsub("^X", "", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "", x)
  x <- tolower(x)
  trimws(x)
}

discovery_features_clean <- normalize_genus_name(
  discovery_features
)

external_features_clean <- normalize_genus_name(
  colnames(external_genus_counts)
)

shared_clean <- intersect(
  discovery_features_clean,
  external_features_clean
)

cat("Shared genera after normalization:", length(shared_clean), "\n")
head(shared_clean, 30)

str(discovery_features)

discovery_features <- colnames(
  readRDS(
    file.path(
      output_path,
      "discovery_genus_training_data.rds"
    )
  )
)

discovery_features <- base::setdiff(
  discovery_features,
  c("Sample", "Group")
)

discovery_features_clean <- normalize_genus_name(discovery_features)

external_features_original <- colnames(external_genus_counts)
external_features_clean <- normalize_genus_name(external_features_original)

shared_genera_clean <- intersect(
  discovery_features_clean,
  external_features_clean
)

missing_genera_clean <- base::setdiff(
  discovery_features_clean,
  external_features_clean
)

cat("Discovery genera:", length(discovery_features_clean), "\n")
cat("External genera:", length(external_features_clean), "\n")
cat("Shared genera:", length(shared_genera_clean), "\n")
cat(
  "Overlap percentage:",
  round(
    100 * length(shared_genera_clean) /
      length(discovery_features_clean),
    2
  ),
  "%\n"
)

shared_genera_report <- base::setdiff(
  shared_genera_clean,
  "unknown"
)


feature_overlap <- data.frame(
  Discovery_Genera = length(discovery_features_clean),
  External_Genera = length(external_features_clean),
  Shared_Genera_Including_Unknown = length(shared_genera_clean),
  Shared_Named_Genera = length(shared_genera_report),
  Overlap_Percentage = round(
    100 * length(shared_genera_clean) /
      length(discovery_features_clean),
    2
  )
)

write.csv(
  feature_overlap,
  file.path(
    output_path,
    "Table_Corrected_Genus_Feature_Overlap.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(Shared_Genus = shared_genera_report),
  file.path(
    output_path,
    "Table_Corrected_Shared_Genera.csv"
  ),
  row.names = FALSE
)
