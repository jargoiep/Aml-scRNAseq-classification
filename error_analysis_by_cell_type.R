#!/usr/bin/env Rscript

# =============================================================================
# ERROR DISTRIBUTION ANALYSIS BY CELL TYPES
# =============================================================================
# 
# Analyzes how errors (specificity and sensitivity) are distributed across
# different cell types as requested by mentor:
# - HSCs (Immature cells in our data)
# - Monocytes
# - Erythroid cells
# - T-cells, B-cells, NK-cells
# - Early myeloid, Dendritic, Other
#
# Also visualizes cell numbers for each cell type to identify types with
# too few cells for meaningful metrics
# =============================================================================

cat("=== Error Distribution Analysis by Cell Type ===\n")
start_time <- Sys.time()

# =============================================================================
# 1. LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(glmnet)
  library(dplyr)
  library(Matrix)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
})

# =============================================================================
# 2. LOAD DATA AND LASSO RESULTS
# =============================================================================

cat("Loading datasets and LASSO results...\n")

# Load datasets
query_seurat <- readRDS("39014174_diet.RDS")

# Check for existing LASSO results
lasso_files <- list.files(pattern = "lasso_regression_results_.*\\.RDS")
if (length(lasso_files) > 0) {
  latest_lasso_file <- lasso_files[order(file.mtime(lasso_files), decreasing = TRUE)][1]
  cat("Loading latest LASSO results from:", latest_lasso_file, "\n")
  lasso_results <- readRDS(latest_lasso_file)
} else {
  cat("No existing LASSO results found. Please run LASSO analysis first.\n")
  cat("Example: Rscript LASSO_analysis.R\n")
  quit(status = 1)
}

# =============================================================================
# 3. DATA PREPARATION
# =============================================================================

# Prepare query data
DefaultAssay(query_seurat) <- "RNA"

# Filter to labeled cells only
labeled_mask <- query_seurat@meta.data$status %in% c('healthy', 'leukemic')
query_labeled <- query_seurat[, labeled_mask]

cat("Labeled cells for analysis:", ncol(query_labeled), "cells\n")

# Get cell types and true labels
cell_types <- query_labeled@meta.data$ct_simple
true_labels <- ifelse(query_labeled@meta.data$status == "leukemic", 1, 0)

# Check available cell types
cat("\nAvailable cell types:\n")
celltype_counts <- table(cell_types)
print(celltype_counts)

# =============================================================================
# 4. FUNCTION DEFINITIONS
# =============================================================================

#' Calculate performance metrics for a specific cell type
calculate_celltype_metrics <- function(predictions, true_labels, probabilities, 
                                     cell_type_mask, cell_type_name) {
  
  # Extract subset for this cell type
  subset_pred <- predictions[cell_type_mask]
  subset_true <- true_labels[cell_type_mask]
  subset_prob <- probabilities[cell_type_mask]
  
  n_cells <- sum(cell_type_mask)
  n_healthy <- sum(subset_true == 0)
  n_malignant <- sum(subset_true == 1)
  
  if (n_cells == 0) {
    return(list(
      cell_type = cell_type_name,
      n_cells = 0,
      n_healthy = 0,
      n_malignant = 0,
      malignant_percentage = 0,
      accuracy = NA,
      sensitivity = NA,
      specificity = NA,
      auc = NA,
      confusion_matrix = NULL,
      true_positive = 0,
      true_negative = 0,
      false_positive = 0,
      false_negative = 0
    ))
  }
  
  # Calculate confusion matrix
  confusion_matrix <- table(Predicted = subset_pred, Actual = subset_true)
  
  # Handle missing categories
  if (!"0" %in% colnames(confusion_matrix)) {
    confusion_matrix <- cbind(confusion_matrix, "0" = 0)
  }
  if (!"1" %in% colnames(confusion_matrix)) {
    confusion_matrix <- cbind(confusion_matrix, "1" = 0)
  }
  if (!"0" %in% rownames(confusion_matrix)) {
    confusion_matrix <- rbind(confusion_matrix, "0" = 0)
  }
  if (!"1" %in% rownames(confusion_matrix)) {
    confusion_matrix <- rbind(confusion_matrix, "1" = 0)
  }
  
  # Reorder to standard format
  confusion_matrix <- confusion_matrix[c("0", "1"), c("0", "1")]
  
  # Calculate metrics
  TP <- confusion_matrix["1", "1"]
  TN <- confusion_matrix["0", "0"] 
  FP <- confusion_matrix["1", "0"]
  FN <- confusion_matrix["0", "1"]
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  
  # Calculate malignant percentage
  malignant_percentage <- (n_malignant / n_cells) * 100
  
  # Calculate AUC if we have both classes
  auc <- NA
  if (length(unique(subset_true)) > 1 && length(subset_prob) > 1) {
    tryCatch({
      auc <- as.numeric(roc(subset_true, subset_prob, quiet = TRUE)$auc)
    }, error = function(e) {
      auc <<- NA
    })
  }
  
  return(list(
    cell_type = cell_type_name,
    n_cells = n_cells,
    n_healthy = n_healthy,
    n_malignant = n_malignant,
    malignant_percentage = malignant_percentage,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    auc = auc,
    confusion_matrix = confusion_matrix,
    true_positive = TP,
    true_negative = TN,
    false_positive = FP,
    false_negative = FN
  ))
}

#' Create comprehensive visualizations
create_celltype_visualizations <- function(results_df) {
  
  # Remove cell types with too few cells for plotting
  results_filtered <- results_df[results_df$n_cells >= 10, ]
  
  # 1. Cell count visualization
  p1 <- ggplot(results_df, aes(x = reorder(cell_type, n_cells), y = n_cells)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 50, color = "red", linetype = "dashed", linewidth = 1) +
    coord_flip() +
    labs(title = "Cell Count by Type", 
         subtitle = "Red line indicates minimum threshold (50 cells)",
         x = "Cell Type", y = "Number of Cells") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  if (nrow(results_filtered) > 0) {
    # 2. Accuracy by cell type
    p2 <- ggplot(results_filtered, aes(x = reorder(cell_type, accuracy), y = accuracy)) +
      geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
      coord_flip() +
      ylim(0, 1) +
      labs(title = "Accuracy by Cell Type", 
           subtitle = "Only cell types with ≥10 cells",
           x = "Cell Type", y = "Accuracy") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    
    # 3. Sensitivity vs Specificity
    p3 <- ggplot(results_filtered, aes(x = specificity, y = sensitivity)) +
      geom_point(aes(size = n_cells, color = cell_type), alpha = 0.7) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = "Sensitivity vs Specificity by Cell Type",
           subtitle = "Point size indicates cell count",
           x = "Specificity", y = "Sensitivity",
           size = "Cell Count", color = "Cell Type") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "bottom")
    
    # 4. Malignant percentage by cell type
    p4 <- ggplot(results_df, aes(x = reorder(cell_type, malignant_percentage), y = malignant_percentage)) +
      geom_bar(stat = "identity", fill = "salmon", alpha = 0.7) +
      coord_flip() +
      labs(title = "Malignant Cell Percentage by Type",
           x = "Cell Type", y = "Malignant Percentage (%)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(list(cell_counts = p1, accuracy = p2, sens_spec = p3, malignant_pct = p4))
  } else {
    return(list(cell_counts = p1))
  }
}

# =============================================================================
# 5. ANALYZE LASSO MODELS
# =============================================================================

cat("\n=== Analyzing Error Distribution Across Cell Types ===\n")

# Get available models
model_names <- names(lasso_results)
cat("Available models:\n")
for (i in seq_along(model_names)) {
  cat(sprintf("%d. %s\n", i, model_names[i]))
}

# Select best performing models to analyze
models_to_analyze <- c("full_pca_weighted", "full_raw_weighted")
if (!all(models_to_analyze %in% model_names)) {
  cat("Some requested models not found. Using available models:\n")
  models_to_analyze <- intersect(models_to_analyze, model_names)
  if (length(models_to_analyze) == 0) {
    models_to_analyze <- model_names[1]  # Use first available
  }
}

cat("Analyzing models:", paste(models_to_analyze, collapse = ", "), "\n")

# Store results for each model
all_model_results <- list()

for (model_name in models_to_analyze) {
  
  cat("\n--- Analyzing Model:", model_name, "---\n")
  
  model_result <- lasso_results[[model_name]]
  
  # Get predictions from model
  if (!is.null(model_result$predictions)) {
    predictions <- model_result$predictions$predictions
    probabilities <- model_result$predictions$probabilities
    cell_ids <- model_result$predictions$cell_ids
    
    # Match predictions to cell types
    query_cell_names <- colnames(query_labeled)
    pred_indices <- match(cell_ids, query_cell_names)
    
    if (any(is.na(pred_indices))) {
      cat("Warning: Some prediction cell IDs not found in query data\n")
      valid_mask <- !is.na(pred_indices)
      pred_indices <- pred_indices[valid_mask]
      predictions <- predictions[valid_mask]
      probabilities <- probabilities[valid_mask]
      cell_ids <- cell_ids[valid_mask]
    }
    
    # Get corresponding cell types and true labels
    pred_cell_types <- cell_types[pred_indices]
    pred_true_labels <- true_labels[pred_indices]
    
    # Analyze each cell type
    unique_celltypes <- unique(pred_cell_types)
    celltype_results <- list()
    
    for (ct in unique_celltypes) {
      ct_mask <- pred_cell_types == ct
      ct_result <- calculate_celltype_metrics(
        predictions = predictions,
        true_labels = pred_true_labels,
        probabilities = probabilities,
        cell_type_mask = ct_mask,
        cell_type_name = ct
      )
      celltype_results[[ct]] <- ct_result
    }
    
    # Convert to data frame
    results_df <- do.call(rbind, lapply(celltype_results, function(x) {
      data.frame(
        cell_type = x$cell_type,
        n_cells = x$n_cells,
        n_healthy = x$n_healthy,
        n_malignant = x$n_malignant,
        malignant_percentage = x$malignant_percentage,
        accuracy = x$accuracy,
        sensitivity = x$sensitivity,
        specificity = x$specificity,
        auc = x$auc,
        true_positive = x$true_positive,
        true_negative = x$true_negative,
        false_positive = x$false_positive,
        false_negative = x$false_negative,
        stringsAsFactors = FALSE
      )
    }))
    
    # Sort by cell count (descending)
    results_df <- results_df[order(results_df$n_cells, decreasing = TRUE), ]
    
    # Print detailed results
    cat("\n=== DETAILED RESULTS FOR", toupper(model_name), "===\n")
    print(results_df)
    
    # =============================================================================
    # 6. CELL TYPE ANALYSIS SUMMARY
    # =============================================================================
    
    cat("\n=== CELL TYPE ANALYSIS SUMMARY ===\n")
    
    # Identify cell types with insufficient data
    cat("Cell types with < 50 cells (metrics may be unreliable):\n")
    low_count_types <- results_df[results_df$n_cells < 50, ]
    if (nrow(low_count_types) > 0) {
      print(low_count_types[, c("cell_type", "n_cells")])
    } else {
      cat("✅ All cell types have sufficient sample size (≥50 cells)\n")
    }
    
    # Performance analysis
    reliable_types <- results_df[results_df$n_cells >= 50, ]
    
    if (nrow(reliable_types) > 0) {
      cat("\n=== PERFORMANCE ANALYSIS (≥50 cells only) ===\n")
      
      cat("🏆 Best performing cell types (accuracy):\n")
      best_accuracy <- reliable_types[order(reliable_types$accuracy, decreasing = TRUE), ]
      print(head(best_accuracy[, c("cell_type", "n_cells", "accuracy", "sensitivity", "specificity")], 3))
      
      cat("\n⚠️ Worst performing cell types (accuracy):\n")
      worst_accuracy <- reliable_types[order(reliable_types$accuracy), ]
      print(head(worst_accuracy[, c("cell_type", "n_cells", "accuracy", "sensitivity", "specificity")], 3))
      
      cat("\n🔍 Cell types with lowest sensitivity (missing malignant cells):\n")
      worst_sensitivity <- reliable_types[order(reliable_types$sensitivity), ]
      print(head(worst_sensitivity[, c("cell_type", "n_cells", "sensitivity", "false_negative", "n_malignant")], 3))
      
      cat("\n🚨 Cell types with lowest specificity (false positives):\n")
      worst_specificity <- reliable_types[order(reliable_types$specificity), ]
      print(head(worst_specificity[, c("cell_type", "n_cells", "specificity", "false_positive", "n_healthy")], 3))
      
      # Immune vs non-immune comparison
      immune_types <- c("T cells", "B cells", "NK cells")
      immune_results <- reliable_types[reliable_types$cell_type %in% immune_types, ]
      nonimmune_results <- reliable_types[!reliable_types$cell_type %in% immune_types, ]
      
      if (nrow(immune_results) > 0 && nrow(nonimmune_results) > 0) {
        cat("\n=== IMMUNE vs NON-IMMUNE COMPARISON ===\n")
        cat("Immune cells average accuracy:", round(mean(immune_results$accuracy, na.rm = TRUE), 3), "\n")
        cat("Non-immune cells average accuracy:", round(mean(nonimmune_results$accuracy, na.rm = TRUE), 3), "\n")
        
        cat("Immune cells average sensitivity:", round(mean(immune_results$sensitivity, na.rm = TRUE), 3), "\n")
        cat("Non-immune cells average sensitivity:", round(mean(nonimmune_results$sensitivity, na.rm = TRUE), 3), "\n")
        
        cat("Immune cells average specificity:", round(mean(immune_results$specificity, na.rm = TRUE), 3), "\n")
        cat("Non-immune cells average specificity:", round(mean(nonimmune_results$specificity, na.rm = TRUE), 3), "\n")
      }
    }
    
    # Store results
    all_model_results[[model_name]] <- list(
      results_df = results_df,
      celltype_results = celltype_results,
      reliable_types = reliable_types,
      low_count_types = low_count_types
    )
    
    # Create and save visualizations
    cat("\nCreating visualizations...\n")
    plots <- create_celltype_visualizations(results_df)
    
    if (length(plots) > 1) {
      pdf(paste0("celltype_error_analysis_", model_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"), 
          width = 14, height = 10)
      grid.arrange(plots$cell_counts, plots$accuracy, plots$sens_spec, plots$malignant_pct, 
                   ncol = 2, nrow = 2)
      dev.off()
    }
    
  } else {
    cat("Warning: No predictions found for model", model_name, "\n")
  }
}

# =============================================================================
# 7. COMPARATIVE ANALYSIS
# =============================================================================

if (length(all_model_results) > 1) {
  cat("\n=== COMPARATIVE ANALYSIS ACROSS MODELS ===\n")
  
  # Compare accuracy across models for each cell type
  comparison_data <- list()
  
  for (model_name in names(all_model_results)) {
    model_data <- all_model_results[[model_name]]$results_df
    model_data$model <- model_name
    comparison_data[[model_name]] <- model_data
  }
  
  comparison_df <- do.call(rbind, comparison_data)
  
  # Create comparison visualization
  comparison_plot <- ggplot(comparison_df[comparison_df$n_cells >= 50, ], 
                           aes(x = cell_type, y = accuracy, fill = model)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    coord_flip() +
    ylim(0, 1) +
    labs(title = "Model Comparison: Accuracy by Cell Type",
         subtitle = "Only cell types with ≥50 cells",
         x = "Cell Type", y = "Accuracy", fill = "Model") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Save comparison plot
  ggsave(paste0("celltype_model_comparison_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"), 
         comparison_plot, width = 12, height = 8)
}

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

# Save comprehensive results
output_results <- list(
  timestamp = Sys.time(),
  models_analyzed = models_to_analyze,
  celltype_analysis = all_model_results,
  summary = "Error distribution analysis by cell type completed"
)

output_file <- paste0("celltype_error_analysis_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RDS")
saveRDS(output_results, output_file)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", output_file, "\n")
cat("Visualizations saved as PDF files\n")

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("Analysis completed in", round(runtime, 2), "minutes\n")

# =============================================================================
# 9. FINAL SUMMARY FOR MENTOR
# =============================================================================

cat("\n=== KEY FINDINGS FOR MENTOR ===\n")

for (model_name in names(all_model_results)) {
  results_df <- all_model_results[[model_name]]$results_df
  reliable_results <- results_df[results_df$n_cells >= 50, ]
  
  cat("\nModel:", model_name, "\n")
  cat("- Total cell types analyzed:", nrow(results_df), "\n")
  cat("- Cell types with ≥50 cells:", nrow(reliable_results), "\n")
  cat("- Cell types with <50 cells:", sum(results_df$n_cells < 50), "\n")
  
  if (nrow(reliable_results) > 0) {
    best_type <- reliable_results[which.max(reliable_results$accuracy), ]
    worst_type <- reliable_results[which.min(reliable_results$accuracy), ]
    
    cat("- Best performing cell type:", best_type$cell_type, 
        sprintf("(%.1f%% accuracy)", best_type$accuracy * 100), "\n")
    cat("- Worst performing cell type:", worst_type$cell_type, 
        sprintf("(%.1f%% accuracy)", worst_type$accuracy * 100), "\n")
    
    # Immune cell summary
    immune_types <- c("T cells", "B cells", "NK cells")
    immune_results <- reliable_results[reliable_results$cell_type %in% immune_types, ]
    
    if (nrow(immune_results) > 0) {
      cat("- Immune cells average accuracy:", 
          sprintf("%.1f%%", mean(immune_results$accuracy, na.rm = TRUE) * 100), "\n")
    }
  }
}

cat("\n✅ Error analysis by cell type completed!\n")
cat("📊 Check the generated PDF files for detailed visualizations.\n")
