

library(dplyr)
library(tidymodels)
library(xgboost) 

external_genus_counts <- readRDS(
  file.path(
    output_path,
    "external_genus_counts_combined.rds"
  )
)

dim(external_genus_counts)

head(colnames(external_genus_counts))

discovery_training <- readRDS(
  file.path(
    output_path,
    "discovery_genus_training_data.rds"
  )
)

discovery_features <- setdiff(
  colnames(discovery_training),
  c("Sample","Group")
)

normalize_genus_name <- function(x){
  
  x <- tolower(x)
  x <- gsub("^g__", "", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("[^a-z0-9_]", "", x)
  trimws(x)
  
}

disc_clean <- normalize_genus_name(discovery_features)

ext_clean <- normalize_genus_name(
  colnames(external_genus_counts)
)

colnames(external_genus_counts) <- ext_clean
names(disc_clean) <- discovery_features


#step 4
external_prediction <- matrix(
  0,
  nrow=nrow(external_genus_counts),
  ncol=length(discovery_features)
)

colnames(external_prediction) <- discovery_features
rownames(external_prediction) <- rownames(external_genus_counts)

#step 5
shared <- intersect(
  disc_clean,
  ext_clean
)

length(shared)
for(g in shared){
  
  discovery_name <- names(disc_clean)[disc_clean==g]
  
  external_prediction[,discovery_name] <-
    external_genus_counts[,g]
  
}

#step 6
external_prediction <- as.data.frame(
  external_prediction
)

external_prediction$Sample <- rownames(external_prediction)

external_prediction <- external_prediction[
  ,
  c("Sample", discovery_features)
]

#step 7
workflow <- readRDS(
  file.path(
    output_path,
    "finalized_genus_xgboost_workflow.rds"
  )
)



workflow <- readRDS(
  file.path(
    output_path,
    "finalized_genus_xgboost_workflow.rds"
  )
)

discovery_training <- readRDS(
  file.path(
    output_path,
    "discovery_genus_training_data.rds"
  )
)

set.seed(123)

fitted_workflow <- fit(
  workflow,
  data = discovery_training
)

pred_prob <- predict(
  fitted_workflow,
  new_data = external_prediction,
  type = "prob"
)

pred_class <- predict(
  fitted_workflow,
  new_data = external_prediction,
  type = "class"
)

saveRDS(
  fitted_workflow,
  file.path(
    output_path,
    "fitted_genus_xgboost_workflow.rds"
  )
)

fitted_workflow <- readRDS(
  file.path(
    output_path,
    "fitted_genus_xgboost_workflow.rds"
  )
)


pred_prob <- predict(
  fitted_workflow,
  new_data = external_prediction,
  type = "prob"
)

pred_class <- predict(
  fitted_workflow,
  new_data = external_prediction,
  type = "class"
) 

external_results <- cbind(
  external_prediction["Sample"],
  pred_prob,
  pred_class
)

head(external_results)

write.csv(
  external_results,
  file.path(
    output_path,
    "External_Genus_Predictions.csv"
  ),
  row.names = FALSE
)
