#!/usr/bin/env Rscript


# NON-IMMUNE CELLS ANALYSIS VERIFICATION SCRIPT

# This script verifies the non-immune cells analysis 

cat("=== NON-IMMUNE CELLS ANALYSIS VERIFICATION ===\n")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})


# 1. LOAD DATA
 

cat("Loading data...\n")

# Load results
scmap_results <- readRDS("scmap_projection_pcafeatures_with_cell_ids_20250808_014756.RDS")
pca_results <- readRDS("Jan-PCA/clonetracer_malignant_predicted_validation.RDS")
query_seurat <- readRDS("39014174_diet.RDS")

# 2. EXTRACT PREDICTIONS AND GROUND TRUTH

cat("Extracting predictions and ground truth...\n")

# Extract predictions
scmap_predictions <- scmap_results$predictions
pca_predictions <- tools::toTitleCase(pca_results@meta.data$predicted_status)
ground_truth <- query_seurat@meta.data$status
cell_types <- query_seurat@meta.data$ct_simple

# Filter to labeled cells only
labeled_mask <- query_seurat@meta.data$status %in% c('healthy', 'leukemic')

# 3. VERIFY CELL TYPE DISTRIBUTION

cat("\n=== CELL TYPE DISTRIBUTION ===\n")

# Get labeled cells only
labeled_data <- data.frame(
  cell_id = colnames(query_seurat)[labeled_mask],
  scmap_pred = scmap_predictions[labeled_mask],
  pca_pred = pca_predictions[labeled_mask],
  ground_truth = ground_truth[labeled_mask],
  cell_type = cell_types[labeled_mask],
  stringsAsFactors = FALSE
)

cat("Total labeled cells:", nrow(labeled_data), "\n")

# Cell type distribution
cell_distribution <- labeled_data %>%
  group_by(cell_type) %>%
  summarise(
    total_cells = n(),
    healthy_cells = sum(ground_truth == "healthy"),
    leukemic_cells = sum(ground_truth == "leukemic"),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_cells))

cat("Cell type distribution:\n")
print(cell_distribution)

# 4. VERIFY IMMUNE CELL EXCLUSION

cat("\n=== IMMUNE CELL EXCLUSION ===\n")

# Define immune cells to exclude
immune_cells <- c("T cells", "B cells", "NK cells")

# Count immune cells
immune_counts <- cell_distribution[cell_distribution$cell_type %in% immune_cells, ]
cat("Immune cells to be excluded:\n")
print(immune_counts)

total_immune <- sum(immune_counts$total_cells)
cat("Total immune cells to exclude:", total_immune, "\n")


# 5. VERIFY NON-IMMUNE CELL DISTRIBUTION


cat("\n=== NON-IMMUNE CELL DISTRIBUTION ===\n")

# Filter out immune cells
non_immune_data <- labeled_data[!labeled_data$cell_type %in% immune_cells, ]

cat("Non-immune labeled cells:", nrow(non_immune_data), "\n")

# Non-immune cell type distribution
non_immune_distribution <- non_immune_data %>%
  group_by(cell_type) %>%
  summarise(
    total_cells = n(),
    healthy_cells = sum(ground_truth == "healthy"),
    leukemic_cells = sum(ground_truth == "leukemic"),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_cells))

cat("\nNon-immune cell type distribution:\n")
print(non_immune_distribution)

# 6. VERIFY MALIGNANT CELL PROPORTION

cat("\n=== MALIGNANT CELL PROPORTION ===\n")

total_non_immune <- nrow(non_immune_data)
malignant_non_immune <- sum(non_immune_data$ground_truth == "leukemic")
healthy_non_immune <- sum(non_immune_data$ground_truth == "healthy")

malignant_proportion <- malignant_non_immune / total_non_immune * 100
healthy_proportion <- healthy_non_immune / total_non_immune * 100

cat("Non-immune cell breakdown:\n")
cat("Total cells:", total_non_immune, "\n")
cat("Malignant cells:", malignant_non_immune, "\n")
cat("Healthy cells:", healthy_non_immune, "\n")
cat("Malignant proportion:", malignant_proportion, "%\n")
cat("Healthy proportion:", healthy_proportion, "%\n")

# 7. VERIFY SCMAP PERFORMANCE ON NON-IMMUNE CELLS


cat("\n=== SCMAP NON-IMMUNE PERFORMANCE ===\n")

# Calculate SCMap confusion matrix for non-immune cells
scmap_non_immune_conf <- table(
  Predicted = non_immune_data$scmap_pred,
  Actual = non_immune_data$ground_truth
)

cat("SCMap confusion matrix (non-immune cells):\n")
print(scmap_non_immune_conf)

# Extract values
scmap_ni_tp <- scmap_non_immune_conf["Malignant", "leukemic"]
scmap_ni_tn <- scmap_non_immune_conf["Healthy", "healthy"]
scmap_ni_fp <- scmap_non_immune_conf["Malignant", "healthy"]
scmap_ni_fn <- scmap_non_immune_conf["Healthy", "leukemic"]

# Calculate metrics
scmap_ni_accuracy <- (scmap_ni_tp + scmap_ni_tn) / total_non_immune * 100
scmap_ni_sensitivity <- scmap_ni_tp / (scmap_ni_tp + scmap_ni_fn) * 100
scmap_ni_specificity <- scmap_ni_tn / (scmap_ni_tn + scmap_ni_fp) * 100

cat("\nSCMap non-immune metrics:\n")
cat("Accuracy:", scmap_ni_accuracy, "%\n")
cat("Sensitivity:", scmap_ni_sensitivity, "%\n")
cat("Specificity:", scmap_ni_specificity, "%\n")

# 8. VERIFY PCA PERFORMANCE ON NON-IMMUNE CELLS

cat("\n=== PCA NON-IMMUNE PERFORMANCE ===\n")

# Calculate PCA confusion matrix for non-immune cells
pca_non_immune_conf <- table(
  Predicted = non_immune_data$pca_pred,
  Actual = non_immune_data$ground_truth
)

cat("PCA confusion matrix (non-immune cells):\n")
print(pca_non_immune_conf)

# Extract values
pca_ni_tp <- pca_non_immune_conf["Leukemic", "leukemic"]
pca_ni_tn <- pca_non_immune_conf["Healthy", "healthy"]
pca_ni_fp <- pca_non_immune_conf["Leukemic", "healthy"]
pca_ni_fn <- pca_non_immune_conf["Healthy", "leukemic"]

# Calculate metrics
pca_ni_accuracy <- (pca_ni_tp + pca_ni_tn) / total_non_immune * 100
pca_ni_sensitivity <- pca_ni_tp / (pca_ni_tp + pca_ni_fn) * 100
pca_ni_specificity <- pca_ni_tn / (pca_ni_tn + pca_ni_fp) * 100

cat("\nPCA non-immune metrics:\n")
cat("Accuracy:", pca_ni_accuracy, "%\n")
cat("Sensitivity:", pca_ni_sensitivity, "%\n")
cat("Specificity:", pca_ni_specificity, "%\n")

# 9. SUMMARY


cat("=== VERIFICATION SUMMARY ===\n")

cat("Final Results:\n")
cat("Total non-immune cells:", total_non_immune, "\n")
cat("Malignant proportion:", malignant_proportion, "%\n")
cat("SCMap non-immune accuracy:", scmap_ni_accuracy, "%\n")
cat("PCA non-immune accuracy:", pca_ni_accuracy, "%\n")
cat("Accuracy difference (SCMap - PCA):", scmap_ni_accuracy - pca_ni_accuracy, "%\n")

cat("\nVerification completed.\n") 

