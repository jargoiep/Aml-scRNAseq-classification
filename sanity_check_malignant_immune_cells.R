#!/usr/bin/env Rscript

# =============================================================================
# SANITY CHECK: MALIGNANT IMMUNE CELLS ANALYSIS
# =============================================================================
# 
# Quick 5-10 minute check to see if there are malignant T-cells and B-cells
# This addresses mentor's question about why sensitivity changed by 3% 
# when evaluating on non-immune cells only
# =============================================================================

cat("=== Sanity Check: Malignant Immune Cells Analysis ===\n")
start_time <- Sys.time()

# =============================================================================
# 1. LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# =============================================================================
# 2. DATA LOADING
# =============================================================================

cat("Loading datasets...\n")

# Load query dataset (clone tracer data)
query_seurat <- readRDS("39014174_diet.RDS")

cat("Query dataset:", ncol(query_seurat), "cells\n")

# =============================================================================
# 3. DATA PREPARATION
# =============================================================================

# Set correct assay
DefaultAssay(query_seurat) <- "RNA"

# Filter to labeled cells only (excluding unsure and NA)
labeled_mask <- query_seurat@meta.data$status %in% c('healthy', 'leukemic')
query_labeled <- query_seurat[, labeled_mask]

cat("Labeled cells (healthy + leukemic):", ncol(query_labeled), "cells\n")

# =============================================================================
# 4. BASIC DATA EXPLORATION
# =============================================================================

cat("\n=== BASIC DATA EXPLORATION ===\n")

# Check available metadata columns
cat("Available metadata columns:\n")
print(colnames(query_labeled@meta.data))

# Check for cell type column
cell_type_col <- NULL
if ("ct_simple" %in% colnames(query_labeled@meta.data)) {
  cell_type_col <- "ct_simple"
  cat("Using 'ct_simple' for cell type information\n")
} else if ("cell_type" %in% colnames(query_labeled@meta.data)) {
  cell_type_col <- "cell_type"
  cat("Using 'cell_type' for cell type information\n")
} else {
  cat("ERROR: No cell type column found (ct_simple or cell_type)\n")
  quit(status = 1)
}

# Get cell types and status
cell_types <- query_labeled@meta.data[[cell_type_col]]
status <- query_labeled@meta.data$status

cat("\nAll unique cell types in dataset:\n")
print(sort(unique(cell_types)))

# =============================================================================
# 5. OVERALL DISTRIBUTION
# =============================================================================

cat("\n=== OVERALL DISTRIBUTION ===\n")

# Overall status distribution
cat("Overall status distribution:\n")
overall_table <- table(status)
print(overall_table)
cat("Total cells:", sum(overall_table), "\n")
cat("Malignant percentage:", round((overall_table["leukemic"]/sum(overall_table))*100, 1), "%\n")

# Cell type distribution
cat("\nCell type distribution:\n")
celltype_table <- table(cell_types)
print(celltype_table)

# =============================================================================
# 6. FOCUS ON IMMUNE CELLS
# =============================================================================

cat("\n=== IMMUNE CELLS ANALYSIS ===\n")

# Define immune cell types (same as in LASSO analysis)
immune_types <- c("T cells", "B cells", "NK cells")

cat("Defined immune cell types:", paste(immune_types, collapse = ", "), "\n")

# Check which immune types are actually present
present_immune_types <- immune_types[immune_types %in% unique(cell_types)]
cat("Present immune cell types:", paste(present_immune_types, collapse = ", "), "\n")

# Identify immune cells
immune_mask <- cell_types %in% immune_types
n_immune_total <- sum(immune_mask)

cat("Total immune cells (T, B, NK):", n_immune_total, "\n")

if (n_immune_total > 0) {
  # Immune cell breakdown by type and status
  immune_cell_types <- cell_types[immune_mask]
  immune_status <- status[immune_mask]
  
  cat("\nImmune cells cross-tabulation (Cell Type vs Status):\n")
  immune_cross_table <- table(immune_cell_types, immune_status)
  print(immune_cross_table)
  
  # =============================================================================
  # 7. MALIGNANT IMMUNE CELLS IDENTIFICATION
  # =============================================================================
  
  cat("\n=== MALIGNANT IMMUNE CELLS IDENTIFICATION ===\n")
  
  # Find malignant immune cells
  malignant_immune_mask <- immune_mask & (status == "leukemic")
  n_malignant_immune <- sum(malignant_immune_mask)
  
  cat("Total malignant immune cells found:", n_malignant_immune, "\n")
  
  if (n_malignant_immune > 0) {
    cat("\n⚠️  MALIGNANT IMMUNE CELLS DETECTED! ⚠️\n")
    
    # Breakdown by cell type
    malignant_immune_types <- cell_types[malignant_immune_mask]
    malignant_immune_breakdown <- table(malignant_immune_types)
    
    cat("\nMalignant immune cells by type:\n")
    print(malignant_immune_breakdown)
    
    # Calculate percentages for each immune type
    cat("\n=== DETAILED PERCENTAGES ===\n")
    for (immune_type in present_immune_types) {
      n_total_type <- sum(cell_types == immune_type)
      n_healthy_type <- sum(cell_types == immune_type & status == "healthy")
      n_malignant_type <- sum(cell_types == immune_type & status == "leukemic")
      
      if (n_total_type > 0) {
        pct_malignant <- (n_malignant_type / n_total_type) * 100
        pct_healthy <- (n_healthy_type / n_total_type) * 100
        
        cat(sprintf("📊 %s:\n", immune_type))
        cat(sprintf("   - Total: %d cells\n", n_total_type))
        cat(sprintf("   - Healthy: %d cells (%.1f%%)\n", n_healthy_type, pct_healthy))
        cat(sprintf("   - Malignant: %d cells (%.1f%%)\n", n_malignant_type, pct_malignant))
        cat("\n")
      }
    }
    
    # Show some example cell IDs
    malignant_immune_cells <- colnames(query_labeled)[malignant_immune_mask]
    cat("Example malignant immune cell IDs (first 10):\n")
    print(head(malignant_immune_cells, 10))
    
  } else {
    cat("\n✅ NO MALIGNANT IMMUNE CELLS FOUND\n")
    cat("All immune cells are classified as healthy\n")
  }
  
  # =============================================================================
  # 8. SENSITIVITY CHANGE EXPLANATION
  # =============================================================================
  
  cat("\n=== SENSITIVITY CHANGE EXPLANATION ===\n")
  
  # Calculate what happens when we remove immune cells
  non_immune_mask <- !immune_mask
  
  # Original dataset
  original_healthy <- sum(status == "healthy")
  original_malignant <- sum(status == "leukemic")
  
  # Non-immune dataset
  nonimmune_healthy <- sum(status[non_immune_mask] == "healthy")
  nonimmune_malignant <- sum(status[non_immune_mask] == "leukemic")
  
  # Removed cells
  removed_healthy <- original_healthy - nonimmune_healthy
  removed_malignant <- original_malignant - nonimmune_malignant
  
  cat("Impact of removing immune cells:\n")
  cat(sprintf("- Healthy cells removed: %d (%.1f%% of original healthy)\n", 
             removed_healthy, (removed_healthy/original_healthy)*100))
  cat(sprintf("- Malignant cells removed: %d (%.1f%% of original malignant)\n", 
             removed_malignant, (removed_malignant/original_malignant)*100))
  
  if (removed_malignant > 0) {
    cat("\n🔍 EXPLANATION FOR SENSITIVITY CHANGE:\n")
    cat("Sensitivity changed because we removed", removed_malignant, "malignant immune cells\n")
    cat("These were likely misclassified as healthy in the full dataset model,\n")
    cat("so removing them slightly improved the sensitivity calculation.\n")
  } else {
    cat("\n🤔 SENSITIVITY CHANGE UNEXPLAINED:\n")
    cat("No malignant immune cells were removed, so the 3% sensitivity change\n")
    cat("might be due to other factors or data processing differences.\n")
  }
  
} else {
  cat("No immune cells found in the dataset.\n")
}

# =============================================================================
# 9. SUMMARY STATISTICS
# =============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")

# Create summary table
summary_stats <- data.frame(
  Category = c("Total Cells", "Healthy Cells", "Malignant Cells", 
               "Immune Cells", "Malignant Immune Cells", "Non-Immune Cells"),
  Count = c(
    length(status),
    sum(status == "healthy"),
    sum(status == "leukemic"),
    sum(immune_mask),
    if(exists("malignant_immune_mask")) sum(malignant_immune_mask) else 0,
    sum(!immune_mask)
  ),
  stringsAsFactors = FALSE
)

summary_stats$Percentage <- round((summary_stats$Count / length(status)) * 100, 1)

print(summary_stats)

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

# Create results list
results <- list(
  timestamp = Sys.time(),
  total_cells = length(status),
  immune_cells = sum(immune_mask),
  malignant_immune_cells = if(exists("malignant_immune_mask")) sum(malignant_immune_mask) else 0,
  immune_cross_table = if(exists("immune_cross_table")) immune_cross_table else NULL,
  malignant_immune_breakdown = if(exists("malignant_immune_breakdown")) malignant_immune_breakdown else NULL,
  summary_stats = summary_stats,
  immune_types_analyzed = present_immune_types
)

# Save results
output_file <- paste0("malignant_immune_sanity_check_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RDS")
saveRDS(results, output_file)
cat("\nResults saved to:", output_file, "\n")

# =============================================================================
# 11. FINAL SUMMARY FOR MENTOR
# =============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n=== FINAL SUMMARY FOR MENTOR ===\n")
cat("Analysis completed in", round(runtime, 2), "minutes\n")

if (exists("malignant_immune_mask") && sum(malignant_immune_mask) > 0) {
  cat("🔴 FINDING: Malignant immune cells detected in dataset\n")
  cat("   This explains the 3% sensitivity change when removing immune cells\n")
  
  # Quick summary for mentor
  for (immune_type in present_immune_types) {
    n_total <- sum(cell_types == immune_type)
    n_malignant <- sum(cell_types == immune_type & status == "leukemic")
    if (n_total > 0 && n_malignant > 0) {
      pct <- round((n_malignant/n_total)*100, 1)
      cat(sprintf("   - %s: %d/%d (%.1f%%) are malignant\n", immune_type, n_malignant, n_total, pct))
    }
  }
} else {
  cat("🟢 FINDING: No malignant immune cells found\n")
  cat("   Sensitivity change may be due to other factors\n")
}

cat("\n✅ Sanity check completed successfully!\n")
