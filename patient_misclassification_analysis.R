#!/usr/bin/env Rscript

# =============================================================================
# PATIENT MISCLASSIFICATION ANALYSIS
# =============================================================================
# 
# This script analyzes model performance at the patient level to identify
# patient-specific factors affecting classification accuracy.
# 
# Key features:
# - Patient-level accuracy, sensitivity, specificity analysis
# - Misclassification pattern visualization by patient and cell type
# - False positive/negative rate analysis
# - Proper facet scaling as requested by mentor (scales = "free")
# =============================================================================

cat("=== Patient Misclassification Analysis ===\n")
start_time <- Sys.time()

# =============================================================================
# 1. LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(pROC)
})

# =============================================================================
# 2. DATA LOADING
# =============================================================================

cat("Loading datasets...\n")

# Load query dataset (clone tracer data)
query_seurat <- readRDS("39014174_diet.RDS")

# Load LASSO results
lasso_results <- readRDS("lasso_regression_results_20250902_160849.RDS")

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

# Ensure patient ID column exists
if (!"patient" %in% colnames(query_labeled@meta.data)) {
  # Try alternative patient ID columns
  patient_cols <- c("RhapsodyID", "sample_id", "donor_id")
  found_col <- NULL
  for (col in patient_cols) {
    if (col %in% colnames(query_labeled@meta.data)) {
      found_col <- col
      break
    }
  }
  if (!is.null(found_col)) {
    query_labeled@meta.data$patient <- query_labeled@meta.data[[found_col]]
  } else {
    stop("No patient ID column found. Please ensure patient information is available.")
  }
}

# Check for cell type column
cell_type_col <- NULL
if ("ct_simple" %in% colnames(query_labeled@meta.data)) {
  cell_type_col <- "ct_simple"
} else if ("cell_type" %in% colnames(query_labeled@meta.data)) {
  cell_type_col <- "cell_type"
} else {
  stop("No cell type column found (ct_simple or cell_type)")
}

# =============================================================================
# 4. PATIENT-LEVEL ANALYSIS FUNCTION
# =============================================================================

#' Analyze patient-level performance
analyze_patient_performance <- function(seurat_obj, predictions, model_name) {
  
  # Get metadata
  metadata <- seurat_obj@meta.data
  
  # Convert status to numeric (0 = healthy, 1 = leukemic)
  true_labels <- ifelse(metadata$status == "leukemic", 1, 0)
  
  # Get patient and cell type information
  patients <- metadata$patient
  cell_types <- metadata[[cell_type_col]]
  
  # Initialize results
  patient_results <- list()
  
  # Analyze each patient
  unique_patients <- unique(patients)
  cat("Analyzing", length(unique_patients), "patients...\n")
  
  for (patient in unique_patients) {
    patient_mask <- patients == patient
    
    if (sum(patient_mask) == 0) next
    
    # Extract patient data
    patient_pred <- predictions[patient_mask]
    patient_true <- true_labels[patient_mask]
    patient_cell_types <- cell_types[patient_mask]
    
    # Calculate overall patient metrics
    confusion_matrix <- table(Predicted = patient_pred, Actual = patient_true)
    
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
    
    n_cells <- sum(patient_mask)
    n_healthy <- sum(patient_true == 0)
    n_malignant <- sum(patient_true == 1)
    
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
    
    # Calculate rates
    false_positive_rate <- ifelse((TN + FP) > 0, FP / (TN + FP), 0)
    false_negative_rate <- ifelse((TP + FN) > 0, FN / (TP + FN), 0)
    
    # Store results
    patient_results[[patient]] <- list(
      patient = patient,
      n_cells = n_cells,
      n_healthy = n_healthy,
      n_malignant = n_malignant,
      malignant_percentage = (n_malignant / n_cells) * 100,
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      false_positive_rate = false_positive_rate,
      false_negative_rate = false_negative_rate,
      TP = TP, TN = TN, FP = FP, FN = FN,
      confusion_matrix = confusion_matrix
    )
  }
  
  return(patient_results)
}

#' Analyze misclassification patterns by patient and cell type
analyze_misclassification_patterns <- function(seurat_obj, predictions, model_name) {
  
  # Get metadata
  metadata <- seurat_obj@meta.data
  
  # Convert status to numeric (0 = healthy, 1 = leukemic)
  true_labels <- ifelse(metadata$status == "leukemic", 1, 0)
  
  # Get patient and cell type information
  patients <- metadata$patient
  cell_types <- metadata[[cell_type_col]]
  
  # Create detailed misclassification data
  misclass_data <- data.frame(
    patient = patients,
    cell_type = cell_types,
    true_label = true_labels,
    predicted_label = predictions,
    stringsAsFactors = FALSE
  )
  
  # Calculate misclassification types
  misclass_data$classification_type <- with(misclass_data, {
    ifelse(true_label == 0 & predicted_label == 1, "False Positive",
    ifelse(true_label == 1 & predicted_label == 0, "False Negative",
    ifelse(true_label == 0 & predicted_label == 0, "True Negative",
           "True Positive")))
  })
  
  # Focus on misclassifications only
  misclass_only <- misclass_data[misclass_data$classification_type %in% c("False Positive", "False Negative"), ]
  
  # Count misclassifications by patient and cell type
  misclass_summary <- misclass_only %>%
    group_by(patient, cell_type, classification_type) %>%
    summarise(count = n(), .groups = "drop")
  
  return(list(
    full_data = misclass_data,
    misclass_only = misclass_only,
    summary = misclass_summary
  ))
}

# =============================================================================
# 5. CREATE VISUALIZATIONS
# =============================================================================

#' Create patient-level visualizations
create_patient_visualizations <- function(patient_results, misclass_patterns, model_name) {
  
  # Convert patient results to data frame
  patient_df <- do.call(rbind, lapply(patient_results, function(x) {
    data.frame(
      patient = x$patient,
      n_cells = x$n_cells,
      n_malignant = x$n_malignant,
      malignant_percentage = x$malignant_percentage,
      accuracy = x$accuracy,
      sensitivity = x$sensitivity,
      specificity = x$specificity,
      false_positive_rate = x$false_positive_rate,
      false_negative_rate = x$false_negative_rate,
      stringsAsFactors = FALSE
    )
  }))
  
  # Sort by accuracy for better visualization
  patient_df <- patient_df[order(patient_df$accuracy, decreasing = TRUE), ]
  
  # 1. Patient accuracy distribution
  p1 <- ggplot(patient_df, aes(x = accuracy)) +
    geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7, color = "black") +
    geom_vline(xintercept = 0.9, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = "Distribution of Patient Classification Accuracy",
         subtitle = "Red line indicates 90% accuracy threshold",
         x = "Accuracy", y = "Number of Patients") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # 2. Patient cell count vs accuracy (colored by malignant cell count)
  p2 <- ggplot(patient_df, aes(x = n_cells, y = accuracy)) +
    geom_point(aes(color = n_malignant), size = 3, alpha = 0.7) +
    scale_color_gradient(low = "lightblue", high = "darkred", name = "Malignant\nCells") +
    geom_hline(yintercept = 0.9, color = "red", linetype = "dashed") +
    labs(title = "Patient Cell Count vs Classification Accuracy",
         subtitle = "Color indicates number of malignant cells",
         x = "Number of Cells", y = "Accuracy") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # 3. False positive rate by patient
  p3 <- ggplot(patient_df, aes(x = reorder(patient, false_positive_rate), y = false_positive_rate)) +
    geom_bar(stat = "identity", fill = "salmon", alpha = 0.7) +
    coord_flip() +
    labs(title = "False Positive Rate by Patient",
         subtitle = "Healthy cells mistakenly classified as malignant",
         x = "Patient", y = "False Positive Rate") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # 4. False negative rate by patient
  p4 <- ggplot(patient_df, aes(x = reorder(patient, false_negative_rate), y = false_negative_rate)) +
    geom_bar(stat = "identity", fill = "lightcoral", alpha = 0.7) +
    coord_flip() +
    labs(title = "False Negative Rate by Patient",
         subtitle = "Malignant cells mistakenly classified as healthy",
         x = "Patient", y = "False Negative Rate") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # 5. Misclassification pattern by patient and cell type
  # THIS IS THE KEY PLOT THAT MENTOR WANTED FIXED WITH PROPER SCALING
  if (nrow(misclass_patterns$summary) > 0) {
    p5 <- ggplot(misclass_patterns$summary, aes(x = cell_type, y = count, fill = classification_type)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
      facet_wrap(~patient, scales = "free") +  # THIS IS THE FIX: scales = "free"
      coord_flip() +
      labs(title = "Misclassification Pattern by Patient and Cell Type",
           subtitle = "Each patient panel scaled independently for better readability",
           x = "Cell Type", y = "Number of Misclassified Cells", 
           fill = "Misclassification\nType") +
      scale_fill_manual(values = c("False Positive" = "salmon", "False Negative" = "#9999ff")) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 8))
  } else {
    p5 <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No misclassifications found", size = 6) +
      theme_void()
  }
  
  return(list(
    accuracy_dist = p1,
    cell_count_vs_accuracy = p2,
    false_positive_by_patient = p3,
    false_negative_by_patient = p4,
    misclass_pattern = p5,
    patient_df = patient_df
  ))
}

# =============================================================================
# 6. ANALYZE MODELS
# =============================================================================

# Define models to analyze
models_to_analyze <- c("full_pca_weighted", "full_raw_weighted")

all_results <- list()

for (model_name in models_to_analyze) {
  cat("\n=== ANALYZING MODEL:", toupper(model_name), "===\n")
  
  # Check if model results exist
  if (model_name %in% names(lasso_results)) {
    model_data <- lasso_results[[model_name]]
    
    if ("predictions" %in% names(model_data) && "predictions" %in% names(model_data$predictions)) {
      predictions <- model_data$predictions$predictions
      
      # Ensure predictions match the labeled cells
      if (length(predictions) == ncol(query_labeled)) {
        
        # Analyze patient performance
        patient_results <- analyze_patient_performance(query_labeled, predictions, model_name)
        
        # Analyze misclassification patterns
        misclass_patterns <- analyze_misclassification_patterns(query_labeled, predictions, model_name)
        
        # Create visualizations
        plots <- create_patient_visualizations(patient_results, misclass_patterns, model_name)
        
        # Print summary statistics
        cat("\n=== PATIENT PERFORMANCE SUMMARY ===\n")
        patient_df <- plots$patient_df
        
        cat("Total patients analyzed:", nrow(patient_df), "\n")
        cat("Patients with >90% accuracy:", sum(patient_df$accuracy > 0.9, na.rm = TRUE), "\n")
        cat("Patients with >95% accuracy:", sum(patient_df$accuracy > 0.95, na.rm = TRUE), "\n")
        
        if (nrow(patient_df) > 0) {
          best_patient <- patient_df[which.max(patient_df$accuracy), ]
          worst_patient <- patient_df[which.min(patient_df$accuracy), ]
          
          cat("Best performing patient:", best_patient$patient, 
              sprintf("(%.1f%% accuracy)", best_patient$accuracy * 100), "\n")
          cat("Worst performing patient:", worst_patient$patient, 
              sprintf("(%.1f%% accuracy)", worst_patient$accuracy * 100), "\n")
          
          cat("Average accuracy:", sprintf("%.1f%%", mean(patient_df$accuracy, na.rm = TRUE) * 100), "\n")
          cat("Median accuracy:", sprintf("%.1f%%", median(patient_df$accuracy, na.rm = TRUE) * 100), "\n")
        }
        
        # Save plots to PDF
        pdf_file <- paste0("patient_misclassification_analysis_", model_name, "_", 
                          format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
        
        pdf(pdf_file, width = 16, height = 12)
        
        # Page 1: Overview plots
        grid.arrange(plots$accuracy_dist, plots$cell_count_vs_accuracy, 
                    plots$false_positive_by_patient, plots$false_negative_by_patient,
                    ncol = 2, nrow = 2)
        
        # Page 2: Misclassification pattern (the key plot with proper scaling)
        print(plots$misclass_pattern)
        
        dev.off()
        
        cat("Visualizations saved to:", pdf_file, "\n")
        
        # Store results
        all_results[[model_name]] <- list(
          patient_results = patient_results,
          misclass_patterns = misclass_patterns,
          plots = plots,
          summary_stats = list(
            total_patients = nrow(patient_df),
            patients_above_90 = sum(patient_df$accuracy > 0.9, na.rm = TRUE),
            patients_above_95 = sum(patient_df$accuracy > 0.95, na.rm = TRUE),
            mean_accuracy = mean(patient_df$accuracy, na.rm = TRUE),
            median_accuracy = median(patient_df$accuracy, na.rm = TRUE)
          )
        )
        
      } else {
        cat("Warning: Prediction length mismatch for model", model_name, "\n")
      }
    } else {
      cat("Warning: No predictions found for model", model_name, "\n")
    }
  } else {
    cat("Warning: Model", model_name, "not found in LASSO results\n")
  }
}

# =============================================================================
# 7. SAVE COMPREHENSIVE RESULTS
# =============================================================================

# Save all results
output_results <- list(
  timestamp = Sys.time(),
  models_analyzed = names(all_results),
  patient_analysis = all_results,
  summary = "Patient-level misclassification analysis completed with proper plot scaling"
)

output_file <- paste0("patient_misclassification_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RDS")
saveRDS(output_results, output_file)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", output_file, "\n")

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("Analysis completed in", round(runtime, 2), "minutes\n")

# =============================================================================
# 8. FINAL SUMMARY FOR MENTOR
# =============================================================================

cat("\n=== KEY FINDINGS FOR MENTOR ===\n")

for (model_name in names(all_results)) {
  results <- all_results[[model_name]]
  stats <- results$summary_stats
  
  cat("\nModel:", model_name, "\n")
  cat("- Total patients:", stats$total_patients, "\n")
  cat("- Patients with >90% accuracy:", stats$patients_above_90, 
      sprintf("(%.1f%%)", (stats$patients_above_90/stats$total_patients)*100), "\n")
  cat("- Average accuracy:", sprintf("%.1f%%", stats$mean_accuracy * 100), "\n")
  
  # Identify problematic patients
  patient_df <- results$plots$patient_df
  low_acc_patients <- patient_df[patient_df$accuracy < 0.9, ]
  
  if (nrow(low_acc_patients) > 0) {
    cat("- Patients with <90% accuracy:\n")
    for (i in 1:nrow(low_acc_patients)) {
      p <- low_acc_patients[i, ]
      cat(sprintf("  * %s: %.1f%% accuracy (%d cells, %d malignant)\n", 
                  p$patient, p$accuracy*100, p$n_cells, p$n_malignant))
    }
  }
}

cat("\n Patient misclassification analysis completed!\n")
cat("Key improvement: Misclassification pattern plots now use scales='free' for better readability\n")
cat("Each patient panel is scaled independently as requested by mentor\n")
