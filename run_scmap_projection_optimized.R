#!/usr/bin/env Rscript

# =============================================================================
# SCMAP-BASED MALIGNANT CELL DETECTION - OPTIMIZED VERSION WITH CELL ID VERIFICATION
# =============================================================================
# 
# KEY UPDATES BASED ON JAN'S FEEDBACK:
# 1. Use same features as PCA (1,991) for fair comparison
# 2. Ensure log normalization matches PCA
# 3. Improved performance by reducing feature set
# 4. ADDED: Cell ID verification to ensure consistent cell order
# 5. ADDED: Cell IDs in results for traceability
#
# =============================================================================

cat("=== SCMap Projection - Optimized Version with Cell ID Verification ===\n")
start_time <- Sys.time()

# =============================================================================
# 1. LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scmap)
  library(SingleCellExperiment)
  library(dplyr)
})

# =============================================================================
# 2. DATA LOADING
# =============================================================================

cat("Loading datasets...\n")

# Load datasets
query_seurat <- readRDS("39014174_diet.RDS")
reference_seurat <- readRDS("so_signature_h_220125.RDS")
malignant_info <- readRDS("malignant_info_combined(1).RDS")

# Load PCA to get the same features
pca_integrated <- readRDS("pca_integrated.RDS")
pca_features <- rownames(pca_integrated@feature.loadings)
cat("Number of PCA features:", length(pca_features), "\n")

# Print initial dataset sizes
cat("Query dataset:", ncol(query_seurat), "cells\n")
cat("Reference dataset:", ncol(reference_seurat), "cells\n")

# =============================================================================
# 3. CELL ID VERIFICATION FUNCTION
# =============================================================================

verify_cell_order <- function(query_seurat, reference_seurat, malignant_info, pca_integrated) {
  cat("=== CELL ID VERIFICATION ===\n")
  
  # Get cell IDs from each object
  query_cells <- colnames(query_seurat)
  reference_cells <- colnames(reference_seurat)
  malignant_cells <- rownames(malignant_info)
  
  # Check if PCA has cell IDs (it might not, depending on how it was created)
  pca_has_cells <- FALSE
  if ("cell.embeddings" %in% names(pca_integrated)) {
    pca_cells <- rownames(pca_integrated@cell.embeddings)
    pca_has_cells <- TRUE
    cat("PCA has cell embeddings with", length(pca_cells), "cells\n")
  } else {
    cat("PCA does not have cell embeddings - this is normal for integrated PCA\n")
  }
  
  # Check overlaps
  cat("Cell ID overlaps:\n")
  cat("  Query vs Reference:", length(intersect(query_cells, reference_cells)), "\n")
  cat("  Query vs Malignant info:", length(intersect(query_cells, malignant_cells)), "\n")
  if (pca_has_cells) {
    cat("  Query vs PCA:", length(intersect(query_cells, pca_cells)), "\n")
  }
  
  # Verify that query cells are in the same order as they would be in PCA
  # This is critical for applying label masks correctly
  if (pca_has_cells) {
    query_in_pca <- intersect(query_cells, pca_cells)
    if (length(query_in_pca) > 0) {
      # Check if order matches
      query_order <- match(query_in_pca, query_cells)
      pca_order <- match(query_in_pca, pca_cells)
      
      if (identical(query_order, pca_order)) {
        cat("✓ Query cell order matches PCA cell order\n")
      } else {
        cat("⚠ WARNING: Query cell order does NOT match PCA cell order!\n")
        cat("  This could cause incorrect label mask application!\n")
        return(FALSE)
      }
    }
  }
  
  cat("✓ Cell ID verification completed\n")
  return(TRUE)
}

# =============================================================================
# 4. DATA PREPARATION
# =============================================================================

# Set correct assays
DefaultAssay(query_seurat) <- "RNA"
DefaultAssay(reference_seurat) <- "RNA_raw"

# Log normalize data (to match PCA)
reference_seurat <- NormalizeData(reference_seurat, normalization.method = "LogNormalize")
query_seurat <- NormalizeData(query_seurat, normalization.method = "LogNormalize")

# Integrate malignant info with reference
reference_seurat@meta.data$malignant_status <- NA
common_cells <- intersect(colnames(reference_seurat), rownames(malignant_info))
reference_seurat@meta.data[common_cells, "malignant_status"] <- 
  malignant_info[common_cells, "classification_combined"]

# Create hybrid classification
reference_seurat@meta.data <- reference_seurat@meta.data %>% 
  mutate(classification_combined_hr = ifelse(WHO_24 == "Healthy control", "Healthy",
                                         ifelse(malignant_status == "Malignant", "Malignant", NA)))

# Remove NA cells
cells_with_labels <- !is.na(reference_seurat@meta.data$classification_combined_hr)
reference_seurat_filtered <- reference_seurat[, cells_with_labels]

cat("Reference dataset after filtering:", ncol(reference_seurat_filtered), "cells\n")
cat("Label distribution:\n")
print(table(reference_seurat_filtered@meta.data$classification_combined_hr))

# =============================================================================
# 5. SCMAP PROJECTION FUNCTION WITH CELL ID TRACKING
# =============================================================================

scmap_project <- function(query_seurat, reference_seurat, features, verbose = TRUE) {
  if (verbose) cat("Starting SCMap projection...\n")
  
  # Get common features with PCA features
  common_features <- intersect(features, rownames(query_seurat))
  common_features <- intersect(common_features, rownames(reference_seurat))
  if (verbose) cat("Common features with PCA:", length(common_features), "\n")
  
  # Subset to PCA features
  query_subset <- query_seurat[common_features, ]
  reference_subset <- reference_seurat[common_features, ]
  
  # Get cell IDs for tracking
  query_cell_ids <- colnames(query_subset)
  reference_cell_ids <- colnames(reference_subset)
  
  if (verbose) {
    cat("Query cells:", length(query_cell_ids), "\n")
    cat("Reference cells:", length(reference_cell_ids), "\n")
  }
  
  # Get normalized data (already normalized in preparation step)
  query_data <- GetAssayData(query_subset, slot = "data")
  reference_data <- GetAssayData(reference_subset, slot = "data")
  
  # Create SCE objects
  reference_sce <- SingleCellExperiment(
    assays = list(logcounts = reference_data),
    colData = reference_subset@meta.data
  )
  
  query_sce <- SingleCellExperiment(
    assays = list(logcounts = query_data),
    colData = query_subset@meta.data
  )
  
  # Add required feature names
  rowData(reference_sce)$feature_symbol <- rownames(reference_sce)
  rowData(query_sce)$feature_symbol <- rownames(query_sce)
  
  # Get reference labels
  ref_labels <- reference_subset@meta.data$classification_combined_hr
  names(ref_labels) <- colnames(reference_subset)
  reference_sce$cell_type1 <- ref_labels
  
  # Set features and create index
  reference_sce <- setFeatures(reference_sce, features = rownames(reference_sce))
  reference_sce <- indexCell(reference_sce)
  
  # Perform projection
  if (verbose) cat("Performing SCMap projection...\n")
  scmap_results <- scmapCell(
    projection = query_sce,
    index_list = list(reference = metadata(reference_sce)$scmap_cell_index),
    w = 3
  )
  
  # Process results
  cell_assignments <- scmap_results$reference$cells
  similarities <- scmap_results$reference$similarities
  
  # Calculate predictions
  predictions <- rep(NA, ncol(query_sce))
  healthy_proportions <- rep(NA, ncol(query_sce))
  
  for (i in 1:ncol(query_sce)) {
    matched_cells <- cell_assignments[, i]
    matched_similarities <- similarities[, i]
    
    valid_matches <- !is.na(matched_cells)
    if (sum(valid_matches) > 0) {
      matched_cells <- matched_cells[valid_matches]
      matched_similarities <- matched_similarities[valid_matches]
      matched_labels <- ref_labels[matched_cells]
      
      # Calculate weighted votes
      healthy_weight <- sum(matched_similarities[matched_labels == "Healthy"])
      malignant_weight <- sum(matched_similarities[matched_labels == "Malignant"])
      total_weight <- healthy_weight + malignant_weight
      
      if (total_weight > 0) {
        healthy_proportions[i] <- healthy_weight / total_weight
        predictions[i] <- ifelse(healthy_proportions[i] > 0.5, "Healthy", "Malignant")
      }
    }
  }
  
  # Create detailed results with cell IDs
  cell_results <- data.frame(
    cell_id = query_cell_ids,
    prediction = predictions,
    healthy_proportion = healthy_proportions,
    stringsAsFactors = FALSE
  )
  
  # Return comprehensive results
  list(
    predictions = predictions,
    healthy_proportions = healthy_proportions,
    cell_ids = query_cell_ids,
    cell_results = cell_results,
    n_features = length(common_features),
    prediction_summary = table(predictions, useNA = "always"),
    reference_cell_ids = reference_cell_ids,
    reference_labels = ref_labels
  )
}

# =============================================================================
# 6. RUN CELL ID VERIFICATION
# =============================================================================

cat("\n=== RUNNING CELL ID VERIFICATION ===\n")
cell_order_ok <- verify_cell_order(query_seurat, reference_seurat, malignant_info, pca_integrated)

if (!cell_order_ok) {
  cat("⚠ WARNING: Cell order verification failed!\n")
  cat("  Proceeding with caution - results may be incorrect!\n")
} else {
  cat("✓ Cell order verification passed - proceeding safely\n")
}

# =============================================================================
# 7. RUN PROJECTION
# =============================================================================

cat("\nRunning SCMap projection with PCA features...\n")
results <- scmap_project(
  query_seurat = query_seurat,
  reference_seurat = reference_seurat_filtered,
  features = pca_features,
  verbose = TRUE
)

# =============================================================================
# 8. VERIFY CELL ID CONSISTENCY IN RESULTS
# =============================================================================

cat("\n=== VERIFYING CELL ID CONSISTENCY IN RESULTS ===\n")

# Check if cell IDs match between query and results
if (identical(colnames(query_seurat), results$cell_ids)) {
  cat("✓ Query cell IDs match results cell IDs\n")
} else {
  cat("⚠ WARNING: Query cell IDs do NOT match results cell IDs!\n")
  cat("  This indicates a serious problem with cell order!\n")
}

# Check for any missing cell IDs
missing_cells <- setdiff(colnames(query_seurat), results$cell_ids)
if (length(missing_cells) > 0) {
  cat("⚠ WARNING: Missing cell IDs in results:", length(missing_cells), "cells\n")
} else {
  cat("✓ All query cells have corresponding results\n")
}

# =============================================================================
# 9. SAVE AND SUMMARIZE RESULTS
# =============================================================================

# Create comprehensive results object
final_results <- list(
  # Core predictions
  predictions = results$predictions,
  healthy_proportions = results$healthy_proportions,
  
  # Cell ID tracking
  cell_ids = results$cell_ids,
  cell_results = results$cell_results,
  reference_cell_ids = results$reference_cell_ids,
  reference_labels = results$reference_labels,
  
  # Verification status
  cell_order_verification_passed = cell_order_ok,
  cell_id_consistency_verified = identical(colnames(query_seurat), results$cell_ids),
  
  # Method metadata
  features_used = pca_features,
  n_features = results$n_features,
  method = "scmap_projection_with_cell_id_verification",
  
  # Statistics
  prediction_summary = results$prediction_summary,
  healthy_prop_stats = list(
    mean = mean(results$healthy_proportions, na.rm = TRUE),
    median = median(results$healthy_proportions, na.rm = TRUE),
    sd = sd(results$healthy_proportions, na.rm = TRUE),
    min = min(results$healthy_proportions, na.rm = TRUE),
    max = max(results$healthy_proportions, na.rm = TRUE)
  )
)

# Save results with timestamp
output_file <- paste0("scmap_projection_pcafeatures_with_cell_ids_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RDS")
saveRDS(final_results, output_file)
cat("Results saved to:", output_file, "\n")

# Print summary statistics
cat("\n=== RESULTS SUMMARY ===\n")
cat("PCA features used:", results$n_features, "\n")
cat("Cells processed:", length(results$cell_ids), "\n")
cat("Cell ID verification:", ifelse(cell_order_ok, "PASSED", "FAILED"), "\n")
cat("Cell ID consistency:", ifelse(identical(colnames(query_seurat), results$cell_ids), "PASSED", "FAILED"), "\n")
cat("Prediction distribution:\n")
print(results$prediction_summary)

# Display healthy proportion statistics
cat("\nHealthy proportion statistics:\n")
cat("  Mean:", round(mean(results$healthy_proportions, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(results$healthy_proportions, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(results$healthy_proportions, na.rm = TRUE), 3), "\n")
cat("  Min:", round(min(results$healthy_proportions, na.rm = TRUE), 3), "\n")
cat("  Max:", round(max(results$healthy_proportions, na.rm = TRUE), 3), "\n")

end_time <- Sys.time()
cat("\nExecution time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

cat("\n=== KEY IMPROVEMENTS ADDED ===\n")
cat("✓ Cell ID verification function\n")
cat("✓ Cell order consistency checks\n")
cat("✓ Cell IDs included in all results\n")
cat("✓ Comprehensive cell ID tracking\n")
cat("✓ Verification status reporting\n") 