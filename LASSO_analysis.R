#!/usr/bin/env Rscript

# =============================================================================
# LASSO REGRESSION ANALYSIS - STANDALONE SCRIPT
# =============================================================================
# 
# This script runs LASSO regression for malignant/healthy cell prediction
# Fixed version that handles PCA dimension issues properly
# Runs on both full dataset and non-immune cells only
#
# =============================================================================

cat("=== LASSO Regression Analysis ===\n")
start_time <- Sys.time()

# =============================================================================
# 1. LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(glmnet)      # For LASSO regression
  library(dplyr)
  library(Matrix)
  library(parallel)
  library(caret)       # For cross-validation and performance metrics
  library(pROC)        # For ROC analysis
  library(doParallel)  # For parallel processing
})

# =============================================================================
# 2. DATA LOADING
# =============================================================================

cat("Loading datasets...\n")

# Load datasets
query_seurat <- readRDS("39014174_diet.RDS")
reference_seurat <- readRDS("so_signature_h_220125.RDS")
malignant_info <- readRDS("malignant_info_combined(1).RDS")

# Load PCA to get the same 1,991 features used in projection methods
pca_integrated <- readRDS("pca_integrated.RDS")
pca_features <- rownames(pca_integrated@feature.loadings)
cat("Number of PCA features to use:", length(pca_features), "\n")

# Print initial dataset sizes
cat("Query dataset:", ncol(query_seurat), "cells\n")
cat("Reference dataset:", ncol(reference_seurat), "cells\n")

# =============================================================================
# 3. DATA PREPARATION FUNCTIONS
# =============================================================================

#' Prepare reference dataset with proper labels
prepare_reference_data <- function(reference_seurat, malignant_info) {
  
  # Set correct assay
  DefaultAssay(reference_seurat) <- "RNA_raw"
  
  # Log normalize data
  reference_seurat <- NormalizeData(reference_seurat, normalization.method = "LogNormalize", verbose = FALSE)
  
  # Integrate malignant info with reference
  reference_seurat@meta.data$malignant_status <- NA
  common_cells <- intersect(colnames(reference_seurat), rownames(malignant_info))
  reference_seurat@meta.data[common_cells, "malignant_status"] <- 
    malignant_info[common_cells, "classification_combined"]
  
  # Create hybrid classification (same as projection methods)
  # Suppress R CMD check warnings for column names
  WHO_24 <- malignant_status <- NULL
  
  reference_seurat@meta.data <- reference_seurat@meta.data %>% 
    mutate(classification_combined_hr = ifelse(WHO_24 == "Healthy control", "Healthy",
                                           ifelse(malignant_status == "Malignant", "Malignant", NA)))
  
  # Remove NA cells
  cells_with_labels <- !is.na(reference_seurat@meta.data$classification_combined_hr)
  reference_seurat_filtered <- reference_seurat[, cells_with_labels]
  
  cat("Reference dataset after filtering:", ncol(reference_seurat_filtered), "cells\n")
  cat("Label distribution:\n")
  print(table(reference_seurat_filtered@meta.data$classification_combined_hr))
  
  return(reference_seurat_filtered)
}

#' Prepare query dataset with ground truth labels
prepare_query_data <- function(query_seurat, include_immune = TRUE) {
  
  # Set correct assay
  DefaultAssay(query_seurat) <- "RNA"
  
  # Log normalize data
  query_seurat <- NormalizeData(query_seurat, normalization.method = "LogNormalize", verbose = FALSE)
  
  # Filter to labeled cells only (excluding unsure and NA)
  labeled_mask <- query_seurat@meta.data$status %in% c('healthy', 'leukemic')
  query_labeled <- query_seurat[, labeled_mask]
  
  if (!include_immune) {
    # Exclude immune cells (T, B, NK cells) as done in projection analysis
    immune_types <- c("T cells", "B cells", "NK cells")
    
    # Check for cell type column (try ct_simple first, then cell_type)
    if ("ct_simple" %in% colnames(query_labeled@meta.data)) {
      non_immune_mask <- !query_labeled@meta.data$ct_simple %in% immune_types
      query_labeled <- query_labeled[, non_immune_mask]
      cat("Non-immune cells after filtering (using ct_simple):", ncol(query_labeled), "cells\n")
    } else if ("cell_type" %in% colnames(query_labeled@meta.data)) {
      non_immune_mask <- !query_labeled@meta.data$cell_type %in% immune_types
      query_labeled <- query_labeled[, non_immune_mask]
      cat("Non-immune cells after filtering (using cell_type):", ncol(query_labeled), "cells\n")
    } else {
      cat("Warning: No cell type column found (ct_simple or cell_type). Using all labeled cells.\n")
      cat("All labeled cells:", ncol(query_labeled), "cells\n")
    }
  } else {
    cat("All labeled cells:", ncol(query_labeled), "cells\n")
  }
  
  # Create binary labels (0 = Healthy, 1 = Malignant)
  query_labeled@meta.data$binary_label <- ifelse(query_labeled@meta.data$status == "leukemic", 1, 0)
  
  cat("Label distribution in query:\n")
  print(table(query_labeled@meta.data$status))
  
  return(query_labeled)
}

#' Extract gene expression matrix and labels for modeling
extract_modeling_data <- function(seurat_obj, features, use_pca = FALSE, n_pcs = 50, pca_integrated = NULL) {
  
  if (use_pca && !is.null(pca_integrated)) {
    # PROJECT onto precomputed PCA space (FAST! - as jan intended)
    cat("Projecting data onto precomputed PCA space...\n")
    
    # Step 1: Get common features between data and precomputed PCA
    common_features <- intersect(features, rownames(seurat_obj))
    pca_loadings <- pca_integrated@feature.loadings[common_features, ]
    
    # Step 2: Get expression data for common features
    expression_data <- GetAssayData(seurat_obj, slot = "data")[common_features, ]
    
    # Step 3: Scale data (important for PCA projection!)
    # Center the data (subtract mean) but don't scale by standard deviation
    expression_scaled <- t(scale(t(expression_data), center = TRUE, scale = FALSE))
    
    # Step 4: Project onto precomputed PCA space using matrix multiplication
    n_pcs_to_use <- min(n_pcs, ncol(pca_loadings))
    pca_coords <- t(pca_loadings[, 1:n_pcs_to_use]) %*% expression_scaled
    expression_matrix <- t(pca_coords)  # Transpose to get cells x PCs format
    
    cat("Projected to", ncol(expression_matrix), "precomputed PCA dimensions\n")
    
  } else if (use_pca) {
    # Fallback: compute fresh PCA if precomputed not available (old approach)
    cat("Computing fresh PCA (precomputed not available)...\n")
    if (!"pca" %in% names(seurat_obj@reductions)) {
      # Compute PCA if not available
      common_features <- intersect(features, rownames(seurat_obj))
      seurat_obj <- ScaleData(seurat_obj, features = common_features, verbose = FALSE)
      
      # Adjust n_pcs to be safe for matrix dimensions
      max_pcs <- min(length(common_features), ncol(seurat_obj)) - 1
      safe_n_pcs <- min(n_pcs, max_pcs, 30)  # Cap at 30 for safety
      
      cat("Computing PCA with", safe_n_pcs, "components (max possible:", max_pcs, ")\n")
      seurat_obj <- RunPCA(seurat_obj, features = common_features, npcs = safe_n_pcs, verbose = FALSE)
    }
    
    # Extract PCA coordinates
    available_pcs <- ncol(seurat_obj@reductions$pca@cell.embeddings)
    n_pcs_to_use <- min(n_pcs, available_pcs)
    expression_matrix <- seurat_obj@reductions$pca@cell.embeddings[, 1:n_pcs_to_use]
    cat("Using fresh PCA coordinates:", ncol(expression_matrix), "dimensions\n")
    
  } else {
    # Use raw gene expression
    common_features <- intersect(features, rownames(seurat_obj))
    cat("Using gene expression features:", length(common_features), "genes\n")
    
    # Get normalized expression data
    expression_data <- GetAssayData(seurat_obj, slot = "data")[common_features, ]
    expression_matrix <- t(as.matrix(expression_data))
  }
  
  # Get labels
  if ("binary_label" %in% colnames(seurat_obj@meta.data)) {
    labels <- seurat_obj@meta.data$binary_label
  } else if ("classification_combined_hr" %in% colnames(seurat_obj@meta.data)) {
    labels <- ifelse(seurat_obj@meta.data$classification_combined_hr == "Malignant", 1, 0)
  } else {
    stop("No suitable label column found")
  }
  
  return(list(
    expression = expression_matrix,
    labels = labels,
    cell_ids = colnames(seurat_obj)
  ))
}

# =============================================================================
# 4. LASSO REGRESSION FUNCTIONS
# =============================================================================

#' Create patient-based folds for cross-validation (MENTOR'S CRITICAL REQUEST)
#' Ensures that cells from the same patient are never split between training and validation
#' @param seurat_obj Seurat object with patient information
#' @param nfolds Number of cross-validation folds
#' @param seed Random seed for reproducibility
create_patient_based_folds <- function(seurat_obj, nfolds = 10, seed = 123) {
  
  set.seed(seed)
  
  # Get patient information
  patient_col <- NULL
  if ("RhapsodyID" %in% colnames(seurat_obj@meta.data)) {
    patient_col <- "RhapsodyID"
  } else if ("patient" %in% colnames(seurat_obj@meta.data)) {
    patient_col <- "patient"
  } else if ("patient_id" %in% colnames(seurat_obj@meta.data)) {
    patient_col <- "patient_id"
  } else {
    # Look for any column that might contain patient information
    potential_cols <- colnames(seurat_obj@meta.data)
    patient_candidates <- potential_cols[grepl("patient|id|sample", potential_cols, ignore.case = TRUE)]
    if (length(patient_candidates) > 0) {
      patient_col <- patient_candidates[1]
      cat("Using column '", patient_col, "' as patient identifier\n")
    } else {
      cat("Warning: No patient column found. Creating dummy patient IDs based on cell index.\n")
      # Create dummy patient IDs (every 100 cells = 1 patient)
      n_cells <- ncol(seurat_obj)
      dummy_patients <- paste0("Patient_", ceiling(seq_len(n_cells) / 100))
      seurat_obj@meta.data$dummy_patient <- dummy_patients
      patient_col <- "dummy_patient"
    }
  }
  
  patient_ids <- seurat_obj@meta.data[[patient_col]]
  unique_patients <- unique(patient_ids)
  n_patients <- length(unique_patients)
  
  cat("Creating patient-based folds:\n")
  cat("- Total cells:", length(patient_ids), "\n")
  cat("- Unique patients:", n_patients, "\n")
  cat("- Requested folds:", nfolds, "\n")
  
  # Adjust nfolds if we have fewer patients than requested folds
  if (n_patients < nfolds) {
    nfolds <- n_patients
    cat("- Adjusted folds to match patient count:", nfolds, "\n")
  }
  
  # Randomly assign patients to folds
  patient_folds <- sample(rep(1:nfolds, length.out = n_patients))
  names(patient_folds) <- unique_patients
  
  # Create cell-level fold assignments based on patient assignments
  cell_folds <- patient_folds[patient_ids]
  
  # Verify no patient is split across folds
  for (patient in unique_patients) {
    patient_cells_folds <- unique(cell_folds[patient_ids == patient])
    if (length(patient_cells_folds) > 1) {
      stop("ERROR: Patient ", patient, " is split across multiple folds!")
    }
  }
  
  cat("✅ Patient-based fold assignment completed successfully\n")
  
  return(list(
    folds = cell_folds,
    nfolds = nfolds,
    n_patients = n_patients,
    patient_col = patient_col
  ))
}

#' Evaluate predictions on a subset of cells (MENTOR'S REQUEST - for non-immune evaluation)
#' @param predictions Original prediction results
#' @param test_data Original test data
#' @param subset_mask Logical vector indicating which cells to include in subset evaluation
#' @param subset_name Name of the subset for reporting
evaluate_predictions_subset <- function(predictions, test_data, subset_mask, subset_name = "subset") {
  
  cat("Evaluating predictions on", subset_name, "subset...\n")
  
  # Apply subset mask to predictions and ground truth
  subset_predictions <- predictions$predictions[subset_mask]
  subset_probabilities <- predictions$probabilities[subset_mask]
  subset_labels <- test_data$labels[subset_mask]
  
  n_subset <- sum(subset_mask)
  cat("- Subset size:", n_subset, "cells\n")
  
  if (n_subset == 0) {
    cat("Warning: Empty subset - cannot evaluate\n")
    return(NULL)
  }
  
  # Calculate performance metrics for subset
  confusion_matrix <- table(Predicted = subset_predictions, Actual = subset_labels)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Calculate sensitivity and specificity with safe division
  TP <- ifelse("1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), 
               confusion_matrix["1", "1"], 0)
  TN <- ifelse("0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), 
               confusion_matrix["0", "0"], 0)
  FP <- ifelse("1" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), 
               confusion_matrix["1", "0"], 0)
  FN <- ifelse("0" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), 
               confusion_matrix["0", "1"], 0)
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  
  # Calculate AUC if possible
  auc <- NA
  if (length(unique(subset_labels)) == 2 && n_subset > 1) {
    tryCatch({
      roc_obj <- roc(subset_labels, subset_probabilities, quiet = TRUE)
      auc <- as.numeric(auc(roc_obj))
    }, error = function(e) {
      cat("Warning: Could not calculate AUC for subset:", e$message, "\n")
    })
  }
  
  cat("- Subset accuracy:", sprintf("%.4f", accuracy), "\n")
  cat("- Subset sensitivity:", sprintf("%.4f", sensitivity), "\n")
  cat("- Subset specificity:", sprintf("%.4f", specificity), "\n")
  if (!is.na(auc)) cat("- Subset AUC:", sprintf("%.4f", auc), "\n")
  
  return(list(
    predictions = subset_predictions,
    probabilities = subset_probabilities,
    confusion_matrix = confusion_matrix,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    auc = auc,
    subset_mask = subset_mask,
    n_cells = n_subset
  ))
}



#' Train LASSO regression model with patient-based cross-validation
#' @param weights Optional vector of observation weights for handling patient bias
#' @param class_weights Optional vector of class weights for handling class imbalance
#' @param use_patient_cv If TRUE, use patient-based cross-validation (MENTOR'S CRITICAL REQUEST)
train_lasso_model <- function(train_data, alpha = 1, nfolds = 10, parallel = TRUE, weights = NULL, class_weights = NULL, use_patient_cv = TRUE, seurat_obj = NULL) {
  
  cat("Training LASSO model with", nrow(train_data$expression), "samples and", 
      ncol(train_data$expression), "features...\n")
  
  # Set up parallel processing if requested
  if (parallel) {
    registerDoParallel(cores = min(4, detectCores() - 1))
  }
  
  # Combine patient weights with class weights if both provided
  final_weights <- weights
  if (!is.null(class_weights)) {
    # Apply class weights based on labels
    class_weight_vector <- ifelse(train_data$labels == 1, class_weights[2], class_weights[1])
    if (!is.null(weights)) {
      # Multiply patient weights by class weights
      final_weights <- weights * class_weight_vector
    } else {
      final_weights <- class_weight_vector
    }
    cat("Applied class weights - Healthy:", class_weights[1], "Malignant:", class_weights[2], "\n")
  }
  
  # Train LASSO with cross-validation to find optimal lambda
  if (use_patient_cv && !is.null(seurat_obj)) {
    cat("Using patient-based cross-validation (MENTOR'S CRITICAL REQUEST)\n")
    
    # Create patient-based folds
    fold_info <- create_patient_based_folds(seurat_obj, nfolds)
    
    # Use custom folds in cv.glmnet
    cv_lasso <- cv.glmnet(
      x = train_data$expression,
      y = train_data$labels,
      family = "binomial",
      alpha = alpha,
      foldid = fold_info$folds,  # Use patient-based folds
      parallel = parallel,
      weights = final_weights,
      type.measure = "class"
    )
    
    cat("Patient-based CV completed with", fold_info$nfolds, "folds\n")
    
  } else {
    cat("Using standard cell-based cross-validation\n")
    
    cv_lasso <- cv.glmnet(
      x = train_data$expression,
      y = train_data$labels,
      family = "binomial",
      alpha = alpha,
      nfolds = nfolds,
      parallel = parallel,
      weights = final_weights,  # Use combined weights
      type.measure = "class"  # Use classification error for CV
    )
  }
  
  # Get the best lambda (using 1se to prevent overfitting as per jan's recommendation)
  best_lambda <- cv_lasso$lambda.1se
  cat("Best lambda (1se to prevent overfitting):", best_lambda, "\n")
  cat("Lambda min would be:", cv_lasso$lambda.min, "\n")
  
  # Check if lambda.1se results in no features (too aggressive regularization)
  temp_model <- glmnet(
    x = train_data$expression,
    y = train_data$labels,
    family = "binomial",
    alpha = alpha,
    lambda = best_lambda,
    weights = weights
  )
  temp_coef <- coef(temp_model, s = best_lambda)
  n_features_1se <- sum(temp_coef[-1] != 0)
  
  if (n_features_1se == 0) {
    cat("Warning: lambda.1se selected 0 features. Falling back to lambda.min for practical use.\n")
    best_lambda <- cv_lasso$lambda.min
    cat("Using lambda.min instead:", best_lambda, "\n")
  }
  
  # Train final model with best lambda
  final_model <- glmnet(
    x = train_data$expression,
    y = train_data$labels,
    family = "binomial",
    alpha = alpha,
    lambda = best_lambda,
    weights = final_weights  # Use combined weights (patient + class)
  )
  
  # Get selected features (non-zero coefficients)
  coefficients <- coef(final_model, s = best_lambda)
  selected_features <- which(coefficients[-1] != 0)  # Exclude intercept
  
  cat("Selected features:", length(selected_features), "out of", ncol(train_data$expression), "\n")
  
  # Handle case where no features are selected
  if (length(selected_features) == 0) {
    cat("Warning: No features selected! This may indicate:\n")
    cat("  - Lambda is too high (over-regularization)\n")
    cat("  - Data scaling issues\n")
    cat("  - Convergence problems\n")
    cat("  - Patient weighting may be too extreme\n")
    
    # For weighted LASSO, try reducing the weight extremes
    if (!is.null(weights)) {
      cat("  - Consider using less extreme patient weights\n")
    }
  }
  
  return(list(
    model = final_model,
    cv_model = cv_lasso,
    best_lambda = best_lambda,
    selected_features = selected_features,
    feature_names = colnames(train_data$expression)[selected_features]
  ))
}

#' Make predictions using trained LASSO model
#' @param optimize_threshold If TRUE, optimize threshold for balanced sensitivity/specificity
predict_lasso <- function(model_result, test_data, optimize_threshold = FALSE) {
  
  # Make predictions
  pred_probs <- predict(model_result$model, newx = test_data$expression, 
                       s = model_result$best_lambda, type = "response")[,1]
  
  # Determine optimal threshold
  threshold <- 0.5  # Default
  if (optimize_threshold && length(unique(test_data$labels)) == 2) {
    # Find threshold that maximizes Youden's J statistic (sensitivity + specificity - 1)
    roc_obj <- roc(test_data$labels, pred_probs, quiet = TRUE)
    optimal_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
    threshold <- optimal_coords$threshold
    cat("Optimized threshold:", sprintf("%.4f", threshold), 
        "(Sens:", sprintf("%.3f", optimal_coords$sensitivity), 
        "Spec:", sprintf("%.3f", optimal_coords$specificity), ")\n")
  }
  
  pred_classes <- ifelse(pred_probs > threshold, 1, 0)
  
  # Calculate performance metrics
  confusion_matrix <- table(Predicted = pred_classes, Actual = test_data$labels)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Calculate sensitivity and specificity (as requested by jan)
  # Handle cases where confusion matrix might not have all categories
  TP <- ifelse("1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), 
               confusion_matrix["1", "1"], 0)
  TN <- ifelse("0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), 
               confusion_matrix["0", "0"], 0)
  FP <- ifelse("1" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), 
               confusion_matrix["1", "0"], 0)
  FN <- ifelse("0" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), 
               confusion_matrix["0", "1"], 0)
  
  # Calculate metrics with safe division
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)  # True Positive Rate
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)  # True Negative Rate
  
  # Additional safety check
  if (is.na(sensitivity)) sensitivity <- 0
  if (is.na(specificity)) specificity <- 0
  
  # Calculate additional metrics
  if (length(unique(test_data$labels)) == 2) {
    roc_obj <- roc(test_data$labels, pred_probs, quiet = TRUE)
    auc <- as.numeric(auc(roc_obj))
  } else {
    auc <- NA
  }
  
  return(list(
    predictions = pred_classes,
    probabilities = pred_probs,
    confusion_matrix = confusion_matrix,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    auc = auc,
    cell_ids = test_data$cell_ids
  ))
}

#' Evaluate predictions on a subset of cells (e.g., non-immune only)
#' This function takes existing predictions and re-evaluates metrics on a subset
evaluate_predictions_subset <- function(predictions, test_data, subset_mask, subset_name = "subset") {
  
  cat("Evaluating predictions on", subset_name, "subset (", sum(subset_mask), "out of", length(subset_mask), "cells)...\n")
  
  # Extract subset predictions and labels
  subset_pred_classes <- predictions$predictions[subset_mask]
  subset_labels <- test_data$labels[subset_mask]
  subset_pred_probs <- predictions$probabilities[subset_mask]
  
  # Calculate performance metrics for subset
  confusion_matrix <- table(Predicted = subset_pred_classes, Actual = subset_labels)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Calculate sensitivity and specificity
  TP <- ifelse("1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), 
               confusion_matrix["1", "1"], 0)
  TN <- ifelse("0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), 
               confusion_matrix["0", "0"], 0)
  FP <- ifelse("1" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), 
               confusion_matrix["1", "0"], 0)
  FN <- ifelse("0" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), 
               confusion_matrix["0", "1"], 0)
  
  # Calculate metrics with safe division
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  
  # Additional safety check
  if (is.na(sensitivity)) sensitivity <- 0
  if (is.na(specificity)) specificity <- 0
  
  # Calculate AUC for subset
  if (length(unique(subset_labels)) == 2) {
    roc_obj <- roc(subset_labels, subset_pred_probs, quiet = TRUE)
    auc <- as.numeric(auc(roc_obj))
  } else {
    auc <- NA
  }
  
  cat("Subset evaluation results:\n")
  cat("- Accuracy:", sprintf("%.4f", accuracy), "\n")
  cat("- Sensitivity:", sprintf("%.4f", sensitivity), "\n") 
  cat("- Specificity:", sprintf("%.4f", specificity), "\n")
  cat("- AUC:", sprintf("%.4f", auc), "\n")
  
  return(list(
    predictions = subset_pred_classes,
    probabilities = subset_pred_probs,
    confusion_matrix = confusion_matrix,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    auc = auc,
    subset_mask = subset_mask,
    subset_name = subset_name
  ))
}

# =============================================================================
# 5. WEIGHTING FUNCTIONS (Patient bias + Class imbalance handling)
# =============================================================================

#' Calculate patient-based weights for handling bias
#' Each patient gets total weight of 1, distributed among their cells
#' @param method "equal" for equal patient weights, "sqrt" for square root dampening
calculate_patient_weights <- function(seurat_obj, method = "sqrt") {
  
  # Get patient information
  if ("RhapsodyID" %in% colnames(seurat_obj@meta.data)) {
    patient_ids <- seurat_obj@meta.data$RhapsodyID
  } else if ("patient" %in% colnames(seurat_obj@meta.data)) {
    patient_ids <- seurat_obj@meta.data$patient
  } else {
    # Create dummy patient IDs if none found
    patient_ids <- paste0("Patient_", 
      rep(seq_len(ceiling(ncol(seurat_obj)/1000)), 
          each = 1000)[seq_len(ncol(seurat_obj))])
  }
  
  # Count cells per patient
  patient_counts <- table(patient_ids)
  
  # Calculate weights based on method
  weights <- numeric(length(patient_ids))
  
  if (method == "equal") {
    # Original method: each patient gets total weight of 1
    for (patient in names(patient_counts)) {
      patient_mask <- patient_ids == patient
      weight_per_cell <- 1 / patient_counts[patient]
      weights[patient_mask] <- weight_per_cell
    }
  } else if (method == "sqrt") {
    # Square root dampening: less extreme weighting
    # Weight inversely proportional to sqrt of cell count
    for (patient in names(patient_counts)) {
      patient_mask <- patient_ids == patient
      weight_per_cell <- 1 / sqrt(patient_counts[patient])
      weights[patient_mask] <- weight_per_cell
    }
    # Normalize so total weight equals number of patients
    weights <- weights * length(patient_counts) / sum(weights)
  } else if (method == "balanced") {
    # Balanced approach: less extreme than equal, more than sqrt
    # Use log dampening to balance patient equity with statistical power
    for (patient in names(patient_counts)) {
      patient_mask <- patient_ids == patient
      # Use log(cell_count + 1) to dampen extreme differences
      weight_per_cell <- 1 / log(patient_counts[patient] + 1)
      weights[patient_mask] <- weight_per_cell
    }
    # Normalize so total weight equals number of patients
    weights <- weights * length(patient_counts) / sum(weights)
  }
  
  cat("Patient weighting summary (method:", method, "):\n")
  cat("- Number of patients:", length(patient_counts), "\n")
  cat("- Cells per patient - Min:", min(patient_counts), "Max:", max(patient_counts), 
      "Mean:", round(mean(patient_counts), 1), "\n")
  cat("- Weight per cell - Min:", sprintf("%.4f", min(weights)), 
      "Max:", sprintf("%.4f", max(weights)), "Ratio:", sprintf("%.1f", max(weights)/min(weights)), "\n")
  
  return(weights)
}

#' Calculate class weights for handling class imbalance
#' @param labels Vector of binary labels (0 = healthy, 1 = malignant)
#' @param method "balanced" for inverse frequency, "custom" for manual specification
#' @param healthy_weight Custom weight for healthy class (only used if method="custom")
#' @param malignant_weight Custom weight for malignant class (only used if method="custom")
calculate_class_weights <- function(labels, method = "balanced", healthy_weight = 1, malignant_weight = 1) {
  
  # Count classes
  class_counts <- table(labels)
  n_healthy <- class_counts["0"]
  n_malignant <- class_counts["1"]
  total <- length(labels)
  
  if (method == "balanced") {
    # Inverse frequency weighting - give more weight to minority class
    weight_healthy <- total / (2 * n_healthy)
    weight_malignant <- total / (2 * n_malignant)
  } else if (method == "custom") {
    # User-specified weights
    weight_healthy <- healthy_weight
    weight_malignant <- malignant_weight
  } else {
    stop("Method must be 'balanced' or 'custom'")
  }
  
  cat("Class weighting summary (method:", method, "):\n")
  cat("- Healthy cells:", n_healthy, "(", sprintf("%.1f%%", 100*n_healthy/total), ") - Weight:", sprintf("%.3f", weight_healthy), "\n")
  cat("- Malignant cells:", n_malignant, "(", sprintf("%.1f%%", 100*n_malignant/total), ") - Weight:", sprintf("%.3f", weight_malignant), "\n")
  cat("- Weight ratio (healthy:malignant):", sprintf("%.2f", weight_healthy/weight_malignant), ":1\n")
  
  return(c(weight_healthy, weight_malignant))
}

#' Calculate differential class weights for full dataset (as requested by mentor)
#' Gives higher weights to non-immune healthy cells only, not all healthy cells
#' @param seurat_obj Seurat object with cell type information
#' @param labels Vector of binary labels (0 = healthy, 1 = malignant)
#' @param nonimmune_healthy_weight Weight for non-immune healthy cells (default 3)
#' @param immune_healthy_weight Weight for immune healthy cells (default 1)
#' @param malignant_weight Weight for malignant cells (default 1)
calculate_differential_class_weights <- function(seurat_obj, labels, 
                                               nonimmune_healthy_weight = 3, 
                                               immune_healthy_weight = 1, 
                                               malignant_weight = 1) {
  
  # Get cell type information
  cell_types <- NULL
  if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
    cell_types <- seurat_obj@meta.data$cell_type
    cat("Using 'cell_type' column for differential weighting\n")
  } else if ("ct_simple" %in% colnames(seurat_obj@meta.data)) {
    cell_types <- seurat_obj@meta.data$ct_simple
    cat("Using 'ct_simple' column for differential weighting\n")
  } else {
    # Check what columns are available
    available_cols <- colnames(seurat_obj@meta.data)
    cat("Available metadata columns:", paste(available_cols, collapse = ", "), "\n")
    
    # Look for any column that might contain cell type information
    potential_ct_cols <- available_cols[grepl("cell|type|ct|cluster", available_cols, ignore.case = TRUE)]
    if (length(potential_ct_cols) > 0) {
      cat("Potential cell type columns found:", paste(potential_ct_cols, collapse = ", "), "\n")
      # Use the first potential column
      cell_types <- seurat_obj@meta.data[[potential_ct_cols[1]]]
      cat("Using column '", potential_ct_cols[1], "' for differential weighting\n")
    } else {
      cat("No suitable cell type column found for differential weighting.\n")
      cat("Available columns: ", paste(available_cols, collapse = ", "), "\n")
      cat("Falling back to standard class weighting.\n")
      return(NULL)  # Return NULL to indicate fallback needed
    }
  }
  
  # Define immune cell types (be flexible with naming)
  immune_types <- c("T cells", "B cells", "NK cells", "T cell", "B cell", "NK cell", 
                    "T-cells", "B-cells", "NK-cells", "Tcells", "Bcells", "NKcells",
                    "CD4", "CD8", "Tcell", "Bcell", "NKcell")
  
  # Create weight vector for each cell
  weights <- numeric(length(labels))
  
  for (i in seq_along(labels)) {
    if (labels[i] == 1) {
      # Malignant cell
      weights[i] <- malignant_weight
    } else {
      # Healthy cell - check if immune or non-immune
      if (cell_types[i] %in% immune_types) {
        # Immune healthy cell
        weights[i] <- immune_healthy_weight
      } else {
        # Non-immune healthy cell (the hard-to-predict ones)
        weights[i] <- nonimmune_healthy_weight
      }
    }
  }
  
  # Count different cell categories for reporting
  n_malignant <- sum(labels == 1)
  n_immune_healthy <- sum(labels == 0 & cell_types %in% immune_types)
  n_nonimmune_healthy <- sum(labels == 0 & !cell_types %in% immune_types)
  total <- length(labels)
  
  cat("Differential class weighting summary:\n")
  cat("- Malignant cells:", n_malignant, "(", sprintf("%.1f%%", 100*n_malignant/total), ") - Weight:", malignant_weight, "\n")
  cat("- Immune healthy cells:", n_immune_healthy, "(", sprintf("%.1f%%", 100*n_immune_healthy/total), ") - Weight:", immune_healthy_weight, "\n")
  cat("- Non-immune healthy cells:", n_nonimmune_healthy, "(", sprintf("%.1f%%", 100*n_nonimmune_healthy/total), ") - Weight:", nonimmune_healthy_weight, "\n")
  cat("- Weight ratios - NonImmuneHealthy:Malignant =", sprintf("%.1f", nonimmune_healthy_weight/malignant_weight), ":1\n")
  cat("- Weight ratios - NonImmuneHealthy:ImmuneHealthy =", sprintf("%.1f", nonimmune_healthy_weight/immune_healthy_weight), ":1\n")
  
  # Show unique cell types found
  unique_types <- unique(cell_types)
  cat("Unique cell types found:", paste(unique_types, collapse = ", "), "\n")
  cat("Cell types classified as immune:", paste(unique_types[unique_types %in% immune_types], collapse = ", "), "\n")
  
  return(weights)
}

# =============================================================================
# 6. MAIN ANALYSIS FUNCTIONS
# =============================================================================

#' Run LASSO analysis
#' @param use_patient_weights If TRUE, use patient-based weighting to handle bias
#' @param use_class_weights If TRUE, use class weighting to handle imbalance
#' @param use_differential_weights If TRUE, use differential weighting (higher for non-immune healthy)
#' @param optimize_threshold If TRUE, optimize prediction threshold for better balance
#' @param evaluate_nonimmune_subset If TRUE, also evaluate performance on non-immune subset (for full dataset models)
#' @param use_patient_cv If TRUE, use patient-based cross-validation for lambda selection (MENTOR'S REQUEST)
run_lasso_analysis <- function(include_immune = TRUE, use_pca = FALSE, pca_integrated = NULL, 
                              use_patient_weights = FALSE, use_class_weights = FALSE, 
                              use_differential_weights = FALSE, optimize_threshold = FALSE,
                              evaluate_nonimmune_subset = FALSE, use_patient_cv = FALSE) {
  
  dataset_name <- ifelse(include_immune, "full", "non_immune")
  method_name <- ifelse(use_pca, "PCA", "raw_genes")
  weight_name <- paste0(
    ifelse(use_patient_weights, "patient_weighted", ""),
    ifelse(use_patient_weights && (use_class_weights || use_differential_weights), "+", ""),
    ifelse(use_class_weights, "class_weighted", ""),
    ifelse(use_differential_weights, "differential_weighted", ""),
    ifelse(!use_patient_weights && !use_class_weights && !use_differential_weights, "standard", "")
  )
  threshold_name <- ifelse(optimize_threshold, "optimized_threshold", "default_threshold")
  
  cat("\n=== LASSO Analysis -", dataset_name, "dataset,", method_name, "(", weight_name, ",", threshold_name, ") ===\n")
  
  # Prepare data
  reference_processed <- prepare_reference_data(reference_seurat, malignant_info)
  query_processed <- prepare_query_data(query_seurat, include_immune = include_immune)
  
  # Filter reference to non-immune cells if needed
  if (!include_immune) {
    immune_types <- c("T cells", "B cells", "NK cells")
    if ("cell_type" %in% colnames(reference_processed@meta.data)) {
      non_immune_ref_mask <- !reference_processed@meta.data$cell_type %in% immune_types
      reference_processed <- reference_processed[, non_immune_ref_mask]
      cat("Reference non-immune cells:", ncol(reference_processed), "cells\n")
    }
  }
  
  # Extract modeling data
  ref_data <- extract_modeling_data(reference_processed, pca_features, use_pca, n_pcs = 50, pca_integrated)
  query_data <- extract_modeling_data(query_processed, pca_features, use_pca, n_pcs = 50, pca_integrated)
  
  # Calculate patient weights if requested
  patient_weights <- NULL
  if (use_patient_weights) {
    cat("Calculating patient-based weights...\n")
    patient_weights <- calculate_patient_weights(reference_processed, method = "sqrt")
  }
  
  # Calculate class weights if requested
  class_weights <- NULL
  differential_weights <- NULL
  
  if (use_class_weights) {
    cat("Calculating class weights for imbalance handling...\n")
    class_weights <- calculate_class_weights(ref_data$labels, method = "balanced")
  } else if (use_differential_weights) {
    cat("Calculating differential weights (higher for non-immune healthy cells)...\n")
    if (!include_immune) {
      cat("Warning: Differential weighting requested but using non-immune dataset only.\n")
      cat("Falling back to standard class weighting for non-immune cells.\n")
      class_weights <- calculate_class_weights(ref_data$labels, method = "balanced")
    } else {
      # Check if reference dataset has cell type information
      has_cell_type <- "cell_type" %in% colnames(reference_processed@meta.data) || 
                       "ct_simple" %in% colnames(reference_processed@meta.data)
      
      if (!has_cell_type) {
        cat("Warning: Reference dataset lacks cell type information for differential weighting.\n")
        cat("Available columns:", paste(colnames(reference_processed@meta.data), collapse = ", "), "\n")
        cat("Falling back to standard class weighting for full dataset.\n")
        class_weights <- calculate_class_weights(ref_data$labels, method = "balanced")
      } else {
        # Use the new differential weighting function
        differential_weights <- calculate_differential_class_weights(reference_processed, ref_data$labels)
      }
    }
  }
  
  # Train model on reference
  cat("Training LASSO model on reference data...\n")
  if (!is.null(differential_weights)) {
    # Use differential weights directly (they're already per-cell weights)
    final_weights <- if (!is.null(patient_weights)) patient_weights * differential_weights else differential_weights
    model_result <- train_lasso_model(ref_data, weights = final_weights, use_patient_cv = use_patient_cv, seurat_obj = reference_processed)
  } else {
    # Use the original approach with class weights
    model_result <- train_lasso_model(ref_data, weights = patient_weights, class_weights = class_weights, use_patient_cv = use_patient_cv, seurat_obj = reference_processed)
  }
  
  # Predict on query
  cat("Making predictions on query data...\n")
  pred_result <- predict_lasso(model_result, query_data, optimize_threshold = optimize_threshold)
  
  # If requested and training on full dataset, also evaluate on non-immune subset
  nonimmune_subset_result <- NULL
  if (evaluate_nonimmune_subset && include_immune) {
    cat("\n--- Evaluating on Non-Immune Subset (as requested by mentor) ---\n")
    
    # Create non-immune mask for query data
    # Get cell type information from query dataset
    immune_types <- c("T cells", "B cells", "NK cells")
    
    if ("ct_simple" %in% colnames(query_processed@meta.data)) {
      query_cell_types <- query_processed@meta.data$ct_simple
      nonimmune_mask <- !(query_cell_types %in% immune_types)
    } else if ("cell_type" %in% colnames(query_processed@meta.data)) {
      query_cell_types <- query_processed@meta.data$cell_type
      nonimmune_mask <- !(query_cell_types %in% immune_types)
    } else {
      cat("Warning: No cell type column found for subset evaluation. Using all cells.\n")
      nonimmune_mask <- rep(TRUE, ncol(query_processed))
    }
    
    # Evaluate predictions on non-immune subset only
    nonimmune_subset_result <- evaluate_predictions_subset(
      predictions = pred_result,
      test_data = query_data,
      subset_mask = nonimmune_mask,
      subset_name = "non-immune"
    )
  }
  
  results <- list(
    model = model_result,
    predictions = pred_result,
    nonimmune_subset_predictions = nonimmune_subset_result,
    dataset_type = dataset_name,
    use_pca = use_pca,
    use_patient_weights = use_patient_weights,
    use_class_weights = use_class_weights,
    use_differential_weights = use_differential_weights,
    optimize_threshold = optimize_threshold,
    evaluate_nonimmune_subset = evaluate_nonimmune_subset,
    n_train_cells = nrow(ref_data$expression),
    n_test_cells = nrow(query_data$expression)
  )
  
  return(results)
}

# =============================================================================
# 6. VERIFICATION FUNCTIONS (as requested by jan)
# =============================================================================

#' Check clone tracer dataset proportions and immune cell removal
#' Addresses jan's questions from Part 3
check_clonetracer_details <- function() {
  cat("\n=== CLONE TRACER DATASET VERIFICATION (jan's Questions) ===\n")
  
  # Load and prepare clone tracer data
  query_processed_full <- prepare_query_data(query_seurat, include_immune = TRUE)
  query_processed_nonimmune <- prepare_query_data(query_seurat, include_immune = FALSE)
  
  # Check full dataset proportions
  full_labels <- query_processed_full@meta.data$status
  full_healthy <- sum(full_labels == "healthy")
  full_malignant <- sum(full_labels == "leukemic")
  full_total <- length(full_labels)
  
  cat("FULL CLONE TRACER DATASET:\n")
  cat("- Healthy cells:", full_healthy, "(", sprintf("%.1f%%", 100*full_healthy/full_total), ")\n")
  cat("- Malignant cells:", full_malignant, "(", sprintf("%.1f%%", 100*full_malignant/full_total), ")\n")
  cat("- Ratio (malignant:healthy):", sprintf("%.2f", full_malignant/full_healthy), ":1\n")
  
  # Check non-immune dataset proportions
  nonimmune_labels <- query_processed_nonimmune@meta.data$status
  nonimmune_healthy <- sum(nonimmune_labels == "healthy")
  nonimmune_malignant <- sum(nonimmune_labels == "leukemic")
  nonimmune_total <- length(nonimmune_labels)
  
  cat("\nNON-IMMUNE CLONE TRACER DATASET:\n")
  cat("- Healthy cells:", nonimmune_healthy, "(", sprintf("%.1f%%", 100*nonimmune_healthy/nonimmune_total), ")\n")
  cat("- Malignant cells:", nonimmune_malignant, "(", sprintf("%.1f%%", 100*nonimmune_malignant/nonimmune_total), ")\n")
  cat("- Ratio (malignant:healthy):", sprintf("%.2f", nonimmune_malignant/nonimmune_healthy), ":1\n")
  
  # Verify immune cell removal
  immune_removed <- full_total - nonimmune_total
  cat("\nIMMUNE CELL REMOVAL VERIFICATION:\n")
  cat("- Cells removed (T, B, NK):", immune_removed, "\n")
  cat("- Percentage removed:", sprintf("%.1f%%", 100*immune_removed/full_total), "\n")
  
  # Check if we have cell type information to verify removal
  if ("ct_simple" %in% colnames(query_processed_full@meta.data)) {
    immune_types <- c("T cells", "B cells", "NK cells")
    immune_cells_full <- sum(query_processed_full@meta.data$ct_simple %in% immune_types)
    cat("- Immune cells in full dataset:", immune_cells_full, "\n")
    
    if ("ct_simple" %in% colnames(query_processed_nonimmune@meta.data)) {
      immune_cells_nonimmune <- sum(query_processed_nonimmune@meta.data$ct_simple %in% immune_types)
      cat("- Immune cells in non-immune dataset:", immune_cells_nonimmune, "\n")
      if (immune_cells_nonimmune == 0) {
        cat("✅ VERIFICATION PASSED: All immune cells properly removed\n")
      } else {
        cat("❌ WARNING: Some immune cells still present in non-immune dataset\n")
      }
    }
  }
  
  cat("\nThis answers jan's question: 'What are the proportions in the clone tracer dataset?'\n")
  
  return(list(
    full_healthy = full_healthy,
    full_malignant = full_malignant,
    nonimmune_healthy = nonimmune_healthy,
    nonimmune_malignant = nonimmune_malignant,
    immune_removed = immune_removed
  ))
}

# =============================================================================
# 7. MAIN EXECUTION
# =============================================================================

# MENTOR'S CRITICAL IMPLEMENTATIONS COMPLETED:
# 1. ✅ Differential weighted evaluation on non-immune subset only
#    - Train on full dataset with differential weighting
#    - Evaluate metrics only on non-immune cells (exclude T/B cells)
#    - This shows how full dataset training performs on non-immune subset
#
# 2. ✅ Patient-based cross-validation for LASSO models  
#    - Modified cv.glmnet to split by patients, not individual cells
#    - Prevents overfitting to specific patients
#    - Ensures better generalization to unseen datasets
#    - This addresses the critical "must do" requirement from mentor

cat("\n=== Starting LASSO Regression Analysis ===\n")

# MENTOR'S REQUEST: Verify clone tracer proportions and immune cell removal
clonetracer_verification <- check_clonetracer_details()

# Run analyses
results <- list()

# Full dataset analysis - Standard LASSO
cat("\n--- Running Full Dataset Analysis (Standard) ---\n")
results$full_raw <- run_lasso_analysis(include_immune = TRUE, use_pca = FALSE, pca_integrated = pca_integrated)
results$full_pca <- run_lasso_analysis(include_immune = TRUE, use_pca = TRUE, pca_integrated = pca_integrated)

# Full dataset analysis - Weighted LASSO (NEW - as requested by mentor)
cat("\n--- Running Full Dataset Analysis (Patient-Weighted with sqrt dampening) ---\n")
results$full_raw_weighted <- run_lasso_analysis(include_immune = TRUE, use_pca = FALSE, pca_integrated = pca_integrated, use_patient_weights = TRUE)
results$full_pca_weighted <- run_lasso_analysis(include_immune = TRUE, use_pca = TRUE, pca_integrated = pca_integrated, use_patient_weights = TRUE)

# Non-immune dataset analysis - Standard LASSO
cat("\n--- Running Non-Immune Dataset Analysis (Standard) ---\n")
results$nonimmune_raw <- run_lasso_analysis(include_immune = FALSE, use_pca = FALSE, pca_integrated = pca_integrated)
results$nonimmune_pca <- run_lasso_analysis(include_immune = FALSE, use_pca = TRUE, pca_integrated = pca_integrated)

# Non-immune dataset analysis - Weighted LASSO (NEW - as requested by mentor)
cat("\n--- Running Non-Immune Dataset Analysis (Patient-Weighted with sqrt dampening) ---\n")
results$nonimmune_raw_weighted <- run_lasso_analysis(include_immune = FALSE, use_pca = FALSE, pca_integrated = pca_integrated, use_patient_weights = TRUE)
results$nonimmune_pca_weighted <- run_lasso_analysis(include_immune = FALSE, use_pca = TRUE, pca_integrated = pca_integrated, use_patient_weights = TRUE)

# Non-immune dataset analysis - Class-Weighted LASSO (NEW - for imbalance handling)
cat("\n--- Running Non-Immune Dataset Analysis (Class-Weighted for imbalance) ---\n")
results$nonimmune_raw_class_weighted <- run_lasso_analysis(include_immune = FALSE, use_pca = FALSE, pca_integrated = pca_integrated, use_class_weights = TRUE)
results$nonimmune_pca_class_weighted <- run_lasso_analysis(include_immune = FALSE, use_pca = TRUE, pca_integrated = pca_integrated, use_class_weights = TRUE)

# Non-immune dataset analysis - Combined Weighting + Threshold Optimization (ULTIMATE)
cat("\n--- Running Non-Immune Dataset Analysis (Combined weighting + threshold optimization) ---\n")
results$nonimmune_raw_combined <- run_lasso_analysis(include_immune = FALSE, use_pca = FALSE, pca_integrated = pca_integrated, 
                                                    use_patient_weights = TRUE, use_class_weights = TRUE, optimize_threshold = TRUE)
results$nonimmune_pca_combined <- run_lasso_analysis(include_immune = FALSE, use_pca = TRUE, pca_integrated = pca_integrated, 
                                                    use_patient_weights = TRUE, use_class_weights = TRUE, optimize_threshold = TRUE)

# MENTOR'S REQUEST: Full dataset with differential weighting (higher weights for non-immune healthy only)
cat("\n--- Running Full Dataset Analysis with Differential Weighting + Patient-Based CV (as requested by mentor) ---\n")
results$full_raw_differential <- run_lasso_analysis(include_immune = TRUE, use_pca = FALSE, pca_integrated = pca_integrated, 
                                                   use_differential_weights = TRUE, use_patient_cv = TRUE)
results$full_pca_differential <- run_lasso_analysis(include_immune = TRUE, use_pca = TRUE, pca_integrated = pca_integrated, 
                                                   use_differential_weights = TRUE, use_patient_cv = TRUE)

# MENTOR'S URGENT REQUEST: Full dataset with differential weighting, evaluated on non-immune subset only
cat("\n--- Running Full Dataset Analysis with Differential Weighting + Non-Immune Subset Evaluation + Patient-Based CV (URGENT - as requested by mentor) ---\n")
cat("This trains on full dataset with differential weighting but evaluates metrics only on non-immune cells\n")
results$full_raw_differential_nonimmune_eval <- run_lasso_analysis(include_immune = TRUE, use_pca = FALSE, pca_integrated = pca_integrated, 
                                                                  use_differential_weights = TRUE, evaluate_nonimmune_subset = TRUE, use_patient_cv = TRUE)
results$full_pca_differential_nonimmune_eval <- run_lasso_analysis(include_immune = TRUE, use_pca = TRUE, pca_integrated = pca_integrated, 
                                                                  use_differential_weights = TRUE, evaluate_nonimmune_subset = TRUE, use_patient_cv = TRUE)

# MENTOR'S REQUEST: Full dataset with combined patient + differential weighting + Patient-Based CV
cat("\n--- Running Full Dataset Analysis with Combined Patient + Differential Weighting + Patient-Based CV ---\n")
results$full_raw_patient_differential <- run_lasso_analysis(include_immune = TRUE, use_pca = FALSE, pca_integrated = pca_integrated, 
                                                           use_patient_weights = TRUE, use_differential_weights = TRUE, use_patient_cv = TRUE)
results$full_pca_patient_differential <- run_lasso_analysis(include_immune = TRUE, use_pca = TRUE, pca_integrated = pca_integrated, 
                                                           use_patient_weights = TRUE, use_differential_weights = TRUE, use_patient_cv = TRUE)

# Summarize results
cat("\n=== LASSO Results Summary ===\n")
summary_df <- data.frame()

for (result_name in names(results)) {
  result <- results[[result_name]]
  
  # Create method name based on all weighting options
  method_name <- "LASSO"
  if (result$use_patient_weights && result$use_class_weights) {
    method_name <- "LASSO_Combined_Weighted"
  } else if (result$use_patient_weights && result$use_differential_weights) {
    method_name <- "LASSO_Patient_Differential_Weighted"
  } else if (result$use_patient_weights) {
    method_name <- "LASSO_Patient_Weighted"
  } else if (result$use_class_weights) {
    method_name <- "LASSO_Class_Weighted"
  } else if (result$use_differential_weights) {
    method_name <- "LASSO_Differential_Weighted"
  }
  
  if (result$optimize_threshold) {
    method_name <- paste0(method_name, "_OptThreshold")
  }
  
  summary_row <- data.frame(
    method = method_name,
    dataset = result$dataset_type,
    evaluation_subset = "full",
    use_pca = result$use_pca,
    use_patient_weights = result$use_patient_weights,
    use_class_weights = ifelse(is.null(result$use_class_weights), FALSE, result$use_class_weights),
    use_differential_weights = ifelse(is.null(result$use_differential_weights), FALSE, result$use_differential_weights),
    optimize_threshold = ifelse(is.null(result$optimize_threshold), FALSE, result$optimize_threshold),
    n_train_cells = result$n_train_cells,
    n_test_cells = result$n_test_cells,
    n_features_selected = length(result$model$selected_features),
    accuracy = result$predictions$accuracy,
    sensitivity = result$predictions$sensitivity,
    specificity = result$predictions$specificity,
    auc = result$predictions$auc,
    stringsAsFactors = FALSE
  )
  
  summary_df <- rbind(summary_df, summary_row)
  
  # Add non-immune subset results if available (MENTOR'S REQUEST)
  if (!is.null(result$nonimmune_subset_predictions)) {
    subset_row <- data.frame(
      method = paste0(method_name, "_NonImmuneEval"),
      dataset = result$dataset_type,
      evaluation_subset = "non_immune",
      use_pca = result$use_pca,
      use_patient_weights = result$use_patient_weights,
      use_class_weights = ifelse(is.null(result$use_class_weights), FALSE, result$use_class_weights),
      use_differential_weights = ifelse(is.null(result$use_differential_weights), FALSE, result$use_differential_weights),
      optimize_threshold = ifelse(is.null(result$optimize_threshold), FALSE, result$optimize_threshold),
      n_train_cells = result$n_train_cells,
      n_test_cells = sum(result$nonimmune_subset_predictions$subset_mask),
      n_features_selected = length(result$model$selected_features),
      accuracy = result$nonimmune_subset_predictions$accuracy,
      sensitivity = result$nonimmune_subset_predictions$sensitivity,
      specificity = result$nonimmune_subset_predictions$specificity,
      auc = result$nonimmune_subset_predictions$auc,
      stringsAsFactors = FALSE
    )
    
    summary_df <- rbind(summary_df, subset_row)
  }
}

print(summary_df)

# Save results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_file <- paste0("lasso_regression_results_", timestamp, ".RDS")
saveRDS(results, results_file)
cat("Results saved to:", results_file, "\n")

# Save summary
summary_file <- paste0("lasso_regression_summary_", timestamp, ".csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Summary saved to:", summary_file, "\n")

# Print key findings
cat("\n=== Key Findings ===\n")
best_result <- summary_df[which.max(summary_df$accuracy), ]
cat("Best performing configuration:\n")
cat("- Dataset:", best_result$dataset, "\n")
cat("- Method:", ifelse(best_result$use_pca, "PCA", "Raw genes"), "\n")
cat("- Accuracy:", sprintf("%.4f", best_result$accuracy), "\n")
cat("- Sensitivity:", sprintf("%.4f", best_result$sensitivity), "\n")
cat("- Specificity:", sprintf("%.4f", best_result$specificity), "\n")
cat("- AUC:", sprintf("%.4f", best_result$auc), "\n")
cat("- Features selected:", best_result$n_features_selected, "\n")

# Show top selected genes for best model
best_model_name <- rownames(best_result)
best_model <- results[[best_model_name]]
if (length(best_model$model$feature_names) > 0) {
  cat("Top selected genes:", paste(head(best_model$model$feature_names, 10), collapse = ", "), "\n")
}

# Extract and save selected genes for all models (as requested by mentor)
cat("\n=== Selected Genes Analysis ===\n")
selected_genes_summary <- list()

for (result_name in names(results)) {
  result <- results[[result_name]]
  if (length(result$model$feature_names) > 0) {
    # Create method name for genes summary
    genes_method_name <- "LASSO"
    if (result$use_patient_weights && result$use_class_weights) {
      genes_method_name <- "LASSO_Combined_Weighted"
    } else if (result$use_patient_weights && result$use_differential_weights) {
      genes_method_name <- "LASSO_Patient_Differential_Weighted"
    } else if (result$use_patient_weights) {
      genes_method_name <- "LASSO_Patient_Weighted"
    } else if (result$use_class_weights) {
      genes_method_name <- "LASSO_Class_Weighted"
    } else if (result$use_differential_weights) {
      genes_method_name <- "LASSO_Differential_Weighted"
    }
    
    selected_genes_summary[[result_name]] <- list(
      method = genes_method_name,
      dataset = result$dataset_type,
      use_pca = result$use_pca,
      use_patient_weights = result$use_patient_weights,
      use_class_weights = ifelse(is.null(result$use_class_weights), FALSE, result$use_class_weights),
      use_differential_weights = ifelse(is.null(result$use_differential_weights), FALSE, result$use_differential_weights),
      optimize_threshold = ifelse(is.null(result$optimize_threshold), FALSE, result$optimize_threshold),
      n_genes = length(result$model$feature_names),
      genes = result$model$feature_names,
      accuracy = result$predictions$accuracy,
      sensitivity = result$predictions$sensitivity,
      specificity = result$predictions$specificity
    )
    
    cat("\n", result_name, ":\n")
    cat("- Method:", genes_method_name, "\n")
    cat("- Dataset:", result$dataset_type, "\n")
    cat("- Features:", ifelse(result$use_pca, "PCA", "Raw genes"), "\n")
    cat("- Genes selected:", length(result$model$feature_names), "\n")
    cat("- Accuracy:", sprintf("%.4f", result$predictions$accuracy), "\n")
    cat("- Sensitivity:", sprintf("%.4f", result$predictions$sensitivity), "\n")
    cat("- Specificity:", sprintf("%.4f", result$predictions$specificity), "\n")
    if (!result$use_pca && length(result$model$feature_names) <= 20) {
      cat("- Selected genes:", paste(result$model$feature_names, collapse = ", "), "\n")
    } else if (!result$use_pca) {
      cat("- Top 10 genes:", paste(head(result$model$feature_names, 10), collapse = ", "), "\n")
    }
  }
}

# Save selected genes summary
genes_file <- paste0("lasso_selected_genes_", timestamp, ".RDS")
saveRDS(selected_genes_summary, genes_file)
cat("\nSelected genes summary saved to:", genes_file, "\n")

end_time <- Sys.time()
cat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# =============================================================================
# MENTOR'S REQUESTED COMPARISON: Full vs Non-immune Models
# =============================================================================

cat("\n=== MENTOR'S REQUESTED COMPARISON ===\n")
cat("Comparing Full Dataset (differential weighting) vs Non-immune Only models\n")
cat("Both approaches predict on NON-IMMUNE cells only for fair comparison\n\n")

# Find the relevant models for comparison
full_diff_raw <- NULL
full_diff_pca <- NULL
nonimmune_combined_raw <- NULL
nonimmune_combined_pca <- NULL

for (name in names(results)) {
  result <- results[[name]]
  if (result$dataset_type == "full" && result$use_differential_weights && !result$use_pca) {
    full_diff_raw <- result
    full_diff_raw_name <- name
  } else if (result$dataset_type == "full" && result$use_differential_weights && result$use_pca) {
    full_diff_pca <- result
    full_diff_pca_name <- name
  } else if (result$dataset_type == "non_immune" && result$use_class_weights && result$use_patient_weights && !result$use_pca) {
    nonimmune_combined_raw <- result
    nonimmune_combined_raw_name <- name
  } else if (result$dataset_type == "non_immune" && result$use_class_weights && result$use_patient_weights && result$use_pca) {
    nonimmune_combined_pca <- result
    nonimmune_combined_pca_name <- name
  }
}

# Compare Raw Gene Models
if (!is.null(full_diff_raw) && !is.null(nonimmune_combined_raw)) {
  cat("RAW GENES COMPARISON:\n")
  cat("Full Dataset (Differential Weighting):\n")
  cat("  - Accuracy:", sprintf("%.4f", full_diff_raw$predictions$accuracy), "\n")
  cat("  - Sensitivity:", sprintf("%.4f", full_diff_raw$predictions$sensitivity), "\n")
  cat("  - Specificity:", sprintf("%.4f", full_diff_raw$predictions$specificity), "\n")
  cat("  - AUC:", sprintf("%.4f", full_diff_raw$predictions$auc), "\n")
  cat("  - Features:", length(full_diff_raw$model$selected_features), "\n")
  
  cat("Non-immune Only (Combined Weighting):\n")
  cat("  - Accuracy:", sprintf("%.4f", nonimmune_combined_raw$predictions$accuracy), "\n")
  cat("  - Sensitivity:", sprintf("%.4f", nonimmune_combined_raw$predictions$sensitivity), "\n")
  cat("  - Specificity:", sprintf("%.4f", nonimmune_combined_raw$predictions$specificity), "\n")
  cat("  - AUC:", sprintf("%.4f", nonimmune_combined_raw$predictions$auc), "\n")
  cat("  - Features:", length(nonimmune_combined_raw$model$selected_features), "\n")
  
  # Determine winner
  if (full_diff_raw$predictions$accuracy > nonimmune_combined_raw$predictions$accuracy) {
    cat("🏆 WINNER (Raw): Full Dataset with Differential Weighting\n")
  } else {
    cat("🏆 WINNER (Raw): Non-immune Only Model\n")
  }
}

cat("\n")

# Compare PCA Models
if (!is.null(full_diff_pca) && !is.null(nonimmune_combined_pca)) {
  cat("PCA COMPARISON:\n")
  cat("Full Dataset (Differential Weighting):\n")
  cat("  - Accuracy:", sprintf("%.4f", full_diff_pca$predictions$accuracy), "\n")
  cat("  - Sensitivity:", sprintf("%.4f", full_diff_pca$predictions$sensitivity), "\n")
  cat("  - Specificity:", sprintf("%.4f", full_diff_pca$predictions$specificity), "\n")
  cat("  - AUC:", sprintf("%.4f", full_diff_pca$predictions$auc), "\n")
  cat("  - Features:", length(full_diff_pca$model$selected_features), "\n")
  
  cat("Non-immune Only (Combined Weighting):\n")
  cat("  - Accuracy:", sprintf("%.4f", nonimmune_combined_pca$predictions$accuracy), "\n")
  cat("  - Sensitivity:", sprintf("%.4f", nonimmune_combined_pca$predictions$sensitivity), "\n")
  cat("  - Specificity:", sprintf("%.4f", nonimmune_combined_pca$predictions$specificity), "\n")
  cat("  - AUC:", sprintf("%.4f", nonimmune_combined_pca$predictions$auc), "\n")
  cat("  - Features:", length(nonimmune_combined_pca$model$selected_features), "\n")
  
  # Determine winner
  if (full_diff_pca$predictions$accuracy > nonimmune_combined_pca$predictions$accuracy) {
    cat("🏆 WINNER (PCA): Full Dataset with Differential Weighting\n")
  } else {
    cat("🏆 WINNER (PCA): Non-immune Only Model\n")
  }
}

cat("\nThis addresses mentor's request: 'Which approach is better - full dataset with differential weighting or non-immune only?'\n")

cat("\n=== LASSO Regression Analysis Complete ===\n")
