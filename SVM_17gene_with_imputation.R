#!/usr/bin/env Rscript

# =============================================================================
# 17-GENE SVM ANALYSIS - WITH MISSING GENE IMPUTATION
# =============================================================================
# STRATEGY:
# - Train on all 17 genes (reference dataset)
# - Test on healthy reference (all 17 genes available)
# - Test on PETTI (ARMH1 missing → impute as zero)
# 
# RATIONALE:
# This tests the actual 17-gene signature as published while handling
# the realistic scenario of missing gene data. Imputing missing genes as
# zero represents "no expression" and tests the robustness of the signature.
# =============================================================================

cat("=== 17-Gene SVM Analysis with Missing Gene Imputation ===\n")
start_time <- Sys.time()

# =============================================================================
# 1. LIBRARY LOADING
# =============================================================================

cat("\n[1/9] Loading required libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(e1071)
  library(dplyr)
  library(Matrix)
  library(data.table)
  library(caret)
  library(pROC)
})
cat("  ✓ All libraries loaded\n")

# =============================================================================
# 2. DEFINE 17-GENE SIGNATURE WITH ALTERNATIVE NAMES
# =============================================================================

cat("\n[2/9] Defining 17-gene signature...\n")

SEVENTEEN_GENES <- c(
  "ARMH1", "CLEC11A", "NREP", "AZU1", "PRAME", "IFITM2", "CFD", "DUSP6", 
  "SCN3A", "CD163", "SOX4", "IFI30", "HOXA9", "CD44", "SRGN", "UBE2C", "CXCL8"
)

# Alternative gene names (for cross-dataset compatibility)
# Source: https://www.genecards.org/
GENE_ALTERNATIVES <- list(
  "ARMH1" = c("ARMH1", "C1orf228", "P40", "NCRNA00082"),
  "CXCL8" = c("CXCL8", "IL8", "IL-8")
)

cat("Full 17-gene signature:\n")
for (i in seq_along(SEVENTEEN_GENES)) {
  cat(sprintf("  %2d. %s", i, SEVENTEEN_GENES[i]))
  if (SEVENTEEN_GENES[i] %in% names(GENE_ALTERNATIVES)) {
    alts <- GENE_ALTERNATIVES[[SEVENTEEN_GENES[i]]][-1]  # Exclude primary name
    if (length(alts) > 0) {
      cat(" (aka:", paste(alts, collapse = ", "), ")")
    }
  }
  cat("\n")
}

# =============================================================================
# 3. HELPER FUNCTION: MAP GENE NAMES TO AVAILABLE NAMES
# =============================================================================

map_gene_names <- function(target_genes, available_genes, alternatives = NULL) {
  # Maps target gene names to their actual names in the dataset
  # Returns a named vector: names = target genes, values = actual gene names in data
  
  gene_mapping <- setNames(rep(NA_character_, length(target_genes)), target_genes)
  
  for (gene in target_genes) {
    # Try exact match first
    if (gene %in% available_genes) {
      gene_mapping[gene] <- gene
      next
    }
    
    # Try alternative names if provided
    if (!is.null(alternatives) && gene %in% names(alternatives)) {
      alt_names <- alternatives[[gene]]
      for (alt in alt_names) {
        if (alt %in% available_genes) {
          gene_mapping[gene] <- alt
          cat("  Mapped", gene, "→", alt, "\n")
          break
        }
      }
    }
  }
  
  return(gene_mapping)
}

# =============================================================================
# 4. LOAD TRAINING DATA
# =============================================================================

cat("\n[3/9] Loading TRAINING data (reference dataset only)...\n")

reference_seurat <- readRDS("so_signature_h_220125.RDS")
malignant_info <- readRDS("malignant_info_combined(1).RDS")

cat("  Setting assay to RNA_raw...\n")
DefaultAssay(reference_seurat) <- "RNA_raw"
reference_seurat <- NormalizeData(reference_seurat, normalization.method = "LogNormalize", verbose = FALSE)

cat("  Assigning labels...\n")
common_cells <- intersect(colnames(reference_seurat), rownames(malignant_info))
reference_seurat@meta.data$malignant_status <- NA
reference_seurat@meta.data[common_cells, "malignant_status"] <- 
  malignant_info[common_cells, "classification_combined"]

reference_seurat@meta.data <- reference_seurat@meta.data %>%
  mutate(classification_combined_hr = ifelse(WHO_24 == "Healthy control", "Healthy",
                                              ifelse(malignant_status == "Malignant", "Malignant", NA)))

labeled_cells <- !is.na(reference_seurat@meta.data$classification_combined_hr)
reference_seurat <- reference_seurat[, labeled_cells]

cat("  ✓ Training dataset prepared:\n")
cat("    - Total cells:", ncol(reference_seurat), "\n")
table_ref <- table(reference_seurat@meta.data$classification_combined_hr)
cat("    - Healthy:", table_ref["Healthy"], "\n")
cat("    - Malignant:", table_ref["Malignant"], "\n")

# =============================================================================
# 5. CHECK GENE AVAILABILITY IN REFERENCE
# =============================================================================

cat("\n[4/9] Checking 17-gene availability in training data...\n")

ref_genes <- rownames(reference_seurat)
ref_gene_mapping <- map_gene_names(SEVENTEEN_GENES, ref_genes, GENE_ALTERNATIVES)

n_found <- sum(!is.na(ref_gene_mapping))
cat("  ✓ Found", n_found, "out of 17 genes in reference\n")

if (n_found < 17) {
  missing_in_ref <- names(ref_gene_mapping)[is.na(ref_gene_mapping)]
  cat("    Missing in reference:", paste(missing_in_ref, collapse = ", "), "\n")
  stop("ERROR: Cannot train 17-gene signature - genes missing in training data!")
}

cat("    All 17 genes available! ✓\n")

# =============================================================================
# 6. EXTRACT 17-GENE EXPRESSION FOR TRAINING
# =============================================================================

cat("\n[5/9] Extracting 17-gene expression for training...\n")

# Extract using mapped gene names
actual_gene_names <- ref_gene_mapping[!is.na(ref_gene_mapping)]
train_expression <- t(as.matrix(LayerData(reference_seurat, 
                                          assay = "RNA_raw", 
                                          layer = "data")[actual_gene_names, ]))

# Ensure column order matches SEVENTEEN_GENES order
colnames(train_expression) <- names(actual_gene_names)

train_labels <- ifelse(reference_seurat@meta.data$classification_combined_hr == "Malignant", 1, 0)

cat("  ✓ Training matrix:", nrow(train_expression), "cells × 17 genes\n")
cat("    - Malignant:", sum(train_labels), "\n")
cat("    - Healthy:", sum(1 - train_labels), "\n")

# =============================================================================
# 6. TRAIN SVM MODEL
# =============================================================================

cat("\n[6/9] Training SVM model on 17 genes...\n")

train_data <- data.frame(train_expression)
colnames(train_data) <- paste0("Gene_", seq_len(17))  # Fixed: always 17 genes
train_data$label <- factor(train_labels, levels = c(0, 1), labels = c("Healthy", "Malignant"))

class_weights <- c(
  "Healthy" = 1 / sum(train_labels == 0),
  "Malignant" = 1 / sum(train_labels == 1)
)
class_weights <- class_weights / min(class_weights)

cat("  Class weights - Healthy:", sprintf("%.4f", class_weights["Healthy"]), 
    "Malignant:", sprintf("%.4f", class_weights["Malignant"]), "\n")

svm_model <- svm(
  label ~ .,
  data = train_data,
  kernel = "linear",
  class.weights = class_weights,
  probability = TRUE,
  scale = TRUE
)

cat("  ✓ SVM trained - Support vectors:", svm_model$tot.nSV, "\n")

# =============================================================================
# 7. EXTERNAL VALIDATION 1: HEALTHY REFERENCE
# =============================================================================

cat("\n[7/9] External Validation 1: Healthy Reference...\n")

cat("  Loading healthy reference expression...\n")
healthy_expression_dt <- fread(cmd = "gunzip -c healthy_reference_expression_CORRECTED_20250918_163850.csv.gz", 
                                showProgress = FALSE)

if (is.character(healthy_expression_dt[[1]])) {
  gene_names <- healthy_expression_dt[[1]]
  healthy_expression_dt <- healthy_expression_dt[, -1, with = FALSE]
} else {
  stop("ERROR: First column should be gene names!")
}

healthy_expression_matrix <- as.matrix(healthy_expression_dt)
rownames(healthy_expression_matrix) <- gene_names
rm(healthy_expression_dt)
gc(verbose = FALSE)

cat("  ✓ Loaded:", nrow(healthy_expression_matrix), "genes ×", ncol(healthy_expression_matrix), "cells\n")

# Check gene availability with mapping
healthy_genes <- rownames(healthy_expression_matrix)
healthy_gene_mapping <- map_gene_names(SEVENTEEN_GENES, healthy_genes, GENE_ALTERNATIVES)

n_healthy_found <- sum(!is.na(healthy_gene_mapping))
healthy_missing <- names(healthy_gene_mapping)[is.na(healthy_gene_mapping)]

cat("  Gene availability:", n_healthy_found, "/ 17\n")

if (length(healthy_missing) > 0) {
  cat("    WARNING: Missing genes:", paste(healthy_missing, collapse = ", "), "\n")
  cat("    These will be imputed as zero\n")
} else {
  cat("    ✓ All 17 genes found!\n")
}

# Extract with imputation if needed
healthy_test_expr <- matrix(0, nrow = ncol(healthy_expression_matrix), ncol = 17)
colnames(healthy_test_expr) <- SEVENTEEN_GENES

for (i in seq_along(SEVENTEEN_GENES)) {
  target_gene <- SEVENTEEN_GENES[i]
  actual_gene <- healthy_gene_mapping[target_gene]
  
  if (!is.na(actual_gene)) {
    # Gene found - extract expression
    healthy_test_expr[, i] <- healthy_expression_matrix[actual_gene, ]
  }
  # else: stays as 0 (imputed)
}

healthy_test_data <- data.frame(healthy_test_expr)
colnames(healthy_test_data) <- paste0("Gene_", seq_len(17))

cat("  Making predictions...\n")
healthy_predictions <- predict(svm_model, healthy_test_data, probability = TRUE)
healthy_probabilities <- attr(healthy_predictions, "probabilities")[, "Malignant"]

healthy_pred_binary <- ifelse(healthy_predictions == "Malignant", 1, 0)
healthy_true_binary <- rep(0, length(healthy_predictions))

healthy_conf <- confusionMatrix(
  factor(healthy_pred_binary, levels = c(0, 1), labels = c("Healthy", "Malignant")),
  factor(healthy_true_binary, levels = c(0, 1), labels = c("Healthy", "Malignant")),
  positive = "Malignant"
)

cat("\n  ✓ HEALTHY REFERENCE RESULTS:\n")
cat("    - Specificity:", sprintf("%.3f", healthy_conf$byClass["Specificity"]), "\n")
cat("    - False Positive Rate:", sprintf("%.3f", 1 - healthy_conf$byClass["Specificity"]), "\n")
cat("    - Correctly classified:", sum(healthy_pred_binary == 0), "/", length(healthy_pred_binary), "\n")

# =============================================================================
# 8. EXTERNAL VALIDATION 2: PETTI DATASET
# =============================================================================

cat("\n[8/9] External Validation 2: PETTI Dataset...\n")

load_10x_data <- function(sample_id, base_path = "Dataset/expression_matrices/") {
  matrix_path <- file.path(base_path, paste0(sample_id, ".matrix.mtx.gz"))
  barcodes_path <- file.path(base_path, paste0(sample_id, ".barcodes.tsv.gz"))
  features_path <- file.path(base_path, paste0(sample_id, ".genes.tsv.gz"))
  
  counts <- readMM(matrix_path)
  barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)$V1
  gene_info <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  gene_names <- make.unique(gene_info$V2)
  rownames(counts) <- gene_names
  colnames(counts) <- paste0(sample_id, "_", barcodes)
  
  return(counts)
}

sample_ids <- c("508084", "548327", "721214", "782328", "809653")
cat("  Loading", length(sample_ids), "PETTI samples...\n")

petti_counts <- load_10x_data(sample_ids[1])
for (i in 2:length(sample_ids)) {
  cat("    [", i, "/", length(sample_ids), "] ", sample_ids[i], "...\n", sep = "")
  sample_counts <- load_10x_data(sample_ids[i])
  common_genes <- intersect(rownames(petti_counts), rownames(sample_counts))
  petti_counts <- petti_counts[common_genes, ]
  sample_counts <- sample_counts[common_genes, ]
  petti_counts <- cbind(petti_counts, sample_counts)
  rm(sample_counts)
  gc(verbose = FALSE)
}

cat("  ✓ Combined:", nrow(petti_counts), "genes ×", ncol(petti_counts), "cells\n")

cat("  Creating Seurat object...\n")
petti_seurat <- CreateSeuratObject(counts = petti_counts, project = "PETTI", min.cells = 0, min.features = 0)
petti_seurat <- NormalizeData(petti_seurat, normalization.method = "LogNormalize", verbose = FALSE)
rm(petti_counts)
gc(verbose = FALSE)

# Load mutation labels
cat("  Loading CB-Sniffer mutation labels...\n")
mutation_files <- list.files(pattern = "final_mutation_labels_corrected.*\\.RDS")
if (length(mutation_files) == 0) {
  stop("ERROR: No mutation labels found!")
}
mutation_labels <- readRDS(mutation_files[which.max(file.mtime(mutation_files))])
cat("    ✓ Loaded", nrow(mutation_labels), "label records\n")

common_cells <- intersect(colnames(petti_seurat), mutation_labels$cell_id)
cat("  Matched", length(common_cells), "cells with labels\n")

# Check gene availability in PETTI with mapping
petti_genes <- rownames(petti_seurat)
petti_gene_mapping <- map_gene_names(SEVENTEEN_GENES, petti_genes, GENE_ALTERNATIVES)

n_petti_found <- sum(!is.na(petti_gene_mapping))
petti_missing <- names(petti_gene_mapping)[is.na(petti_gene_mapping)]

cat("  Gene availability in PETTI:", n_petti_found, "/ 17\n")
if (length(petti_missing) > 0) {
  cat("    ⚠ Missing genes:", paste(petti_missing, collapse = ", "), "\n")
  cat("    → Will be imputed as ZERO (no expression)\n")
} else {
  cat("    ✓ All 17 genes found (with alternative names)!\n")
}

# Extract with imputation if needed
petti_expr <- matrix(0, nrow = length(common_cells), ncol = 17)
colnames(petti_expr) <- SEVENTEEN_GENES
rownames(petti_expr) <- common_cells

for (i in seq_along(SEVENTEEN_GENES)) {
  target_gene <- SEVENTEEN_GENES[i]
  actual_gene <- petti_gene_mapping[target_gene]
  
  if (!is.na(actual_gene)) {
    # Gene found - extract expression
    petti_expr[, i] <- as.matrix(LayerData(petti_seurat, assay = "RNA", layer = "data")[actual_gene, common_cells])
  }
  # else: stays as 0 (imputed)
}

ground_truth <- mutation_labels[match(common_cells, mutation_labels$cell_id), ]

# Filter to labeled cells
labeled_mask <- ground_truth$mutation_status %in% c("malignant", "healthy")
petti_expr <- petti_expr[labeled_mask, , drop = FALSE]
ground_truth <- ground_truth[labeled_mask, ]

cat("  ✓ PETTI labeled cells:", nrow(petti_expr), "\n")
cat("    - Malignant:", sum(ground_truth$mutation_status == "malignant"), "\n")
cat("    - Healthy:", sum(ground_truth$mutation_status == "healthy"), "\n")

# Predict
petti_test_data <- data.frame(petti_expr)
colnames(petti_test_data) <- paste0("Gene_", seq_len(17))

cat("  Making predictions...\n")
petti_predictions <- predict(svm_model, petti_test_data, probability = TRUE)
petti_probabilities <- attr(petti_predictions, "probabilities")[, "Malignant"]

petti_pred_binary <- ifelse(petti_predictions == "Malignant", 1, 0)
petti_true_binary <- ifelse(ground_truth$mutation_status == "malignant", 1, 0)

petti_conf <- confusionMatrix(
  factor(petti_pred_binary, levels = c(0, 1), labels = c("Healthy", "Malignant")),
  factor(petti_true_binary, levels = c(0, 1), labels = c("Healthy", "Malignant")),
  positive = "Malignant"
)

petti_roc <- roc(petti_true_binary, petti_probabilities, quiet = TRUE)
petti_auc <- auc(petti_roc)

cat("\n  ✓ PETTI RESULTS:\n")
cat("    - Accuracy:", sprintf("%.3f", petti_conf$overall["Accuracy"]), "\n")
cat("    - Sensitivity:", sprintf("%.3f", petti_conf$byClass["Sensitivity"]), "\n")
cat("    - Specificity:", sprintf("%.3f", petti_conf$byClass["Specificity"]), "\n")
cat("    - AUC:", sprintf("%.3f", petti_auc), "\n")

# Per-patient
ground_truth$patient <- sapply(strsplit(ground_truth$cell_id, "_"), `[`, 1)
per_patient <- ground_truth %>%
  mutate(pred = petti_pred_binary, truth = petti_true_binary) %>%
  group_by(patient) %>%
  summarise(
    n_cells = n(),
    n_malignant = sum(truth),
    n_healthy = n() - sum(truth),
    accuracy = mean(pred == truth),
    sensitivity = ifelse(sum(truth) > 0, sum(pred == 1 & truth == 1) / sum(truth), NA),
    specificity = ifelse(sum(1 - truth) > 0, sum(pred == 0 & truth == 0) / sum(1 - truth), NA),
    .groups = 'drop'
  )

cat("\n  Per-patient breakdown:\n")
print(per_patient)

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

cat("\n[9/9] Saving results...\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

results <- list(
  model = svm_model,
  genes_signature = SEVENTEEN_GENES,
  gene_alternatives = GENE_ALTERNATIVES,
  gene_mappings = list(
    reference = ref_gene_mapping,
    healthy = healthy_gene_mapping,
    petti = petti_gene_mapping
  ),
  imputation_strategy = list(
    method = "gene_name_mapping_with_zero_imputation",
    description = "Alternative gene names mapped (e.g., ARMH1→C1orf228); missing genes imputed as zero",
    healthy_missing = healthy_missing,
    petti_missing = petti_missing
  ),
  training = list(
    dataset = "so_signature_h_220125.RDS",
    n_cells = nrow(train_expression),
    n_genes = 17,
    n_malignant = sum(train_labels),
    n_healthy = sum(1 - train_labels)
  ),
  healthy_validation = list(
    dataset = "healthy_reference (external)",
    n_cells = ncol(healthy_expression_matrix),
    n_genes_available = n_healthy_found,
    n_genes_imputed = length(healthy_missing),
    specificity = healthy_conf$byClass["Specificity"],
    false_positive_rate = 1 - healthy_conf$byClass["Specificity"]
  ),
  petti_validation = list(
    dataset = "PETTI (external)",
    n_cells = nrow(petti_expr),
    n_genes_available = n_petti_found,
    n_genes_imputed = length(petti_missing),
    imputed_genes = petti_missing,
    n_malignant = sum(petti_true_binary),
    n_healthy = sum(1 - petti_true_binary),
    accuracy = petti_conf$overall["Accuracy"],
    sensitivity = petti_conf$byClass["Sensitivity"],
    specificity = petti_conf$byClass["Specificity"],
    auc = petti_auc,
    per_patient = per_patient
  )
)

results_file <- paste0("svm_17gene_with_imputation_results_", timestamp, ".RDS")
saveRDS(results, results_file)
cat("  ✓ Results:", results_file, "\n")

summary_df <- data.frame(
  Dataset = c("Training", "Healthy_External", "PETTI_External"),
  N_Cells = c(nrow(train_expression), ncol(healthy_expression_matrix), nrow(petti_expr)),
  Genes_Available = c(17, n_healthy_found, n_petti_found),
  Genes_Imputed = c(0, length(healthy_missing), length(petti_missing)),
  Specificity = c(NA, healthy_conf$byClass["Specificity"], petti_conf$byClass["Specificity"]),
  Sensitivity = c(NA, NA, petti_conf$byClass["Sensitivity"]),
  Accuracy = c(NA, NA, petti_conf$overall["Accuracy"]),
  AUC = c(NA, NA, petti_auc)
)

summary_file <- paste0("svm_17gene_with_imputation_summary_", timestamp, ".csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("  ✓ Summary:", summary_file, "\n")

# =============================================================================
# 10. FINAL SUMMARY
# =============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("17-GENE SVM ANALYSIS COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Runtime:", sprintf("%.2f", runtime), "minutes\n\n")
cat("METHOD: 17-Gene SVM Signature (Ohlstrom et al., 2023)\n")
cat("STRATEGY: Gene name mapping + zero imputation for any missing genes\n\n")
cat("TRAINING:\n")
cat("- Dataset: so_signature_h_220125.RDS\n")
cat("- Cells:", nrow(train_expression), "(", sum(1 - train_labels), "healthy,", sum(train_labels), "malignant)\n")
cat("- Genes: All 17 from signature ✓\n\n")
cat("EXTERNAL VALIDATION 1 (Healthy Reference):\n")
cat("- Cells:", ncol(healthy_expression_matrix), "\n")
cat("- Genes found:", n_healthy_found, "/ 17")
if (length(healthy_missing) > 0) {
  cat(" (missing:", paste(healthy_missing, collapse = ", "), ")")
} else {
  cat(" ✓")
}
cat("\n- Specificity:", sprintf("%.1f%%", healthy_conf$byClass["Specificity"] * 100), "\n")
cat("- False Positive Rate:", sprintf("%.1f%%", (1 - healthy_conf$byClass["Specificity"]) * 100), "\n\n")
cat("EXTERNAL VALIDATION 2 (PETTI):\n")
cat("- Cells:", nrow(petti_expr), "(", sum(1 - petti_true_binary), "healthy,", sum(petti_true_binary), "malignant)\n")
cat("- Genes found:", n_petti_found, "/ 17")
if (length(petti_missing) > 0) {
  cat(" (missing:", paste(petti_missing, collapse = ", "), ")")
} else {
  cat(" ✓ (with alternative names)")
}
cat("\n- Accuracy:", sprintf("%.1f%%", petti_conf$overall["Accuracy"] * 100), "\n")
cat("- Sensitivity:", sprintf("%.1f%%", petti_conf$byClass["Sensitivity"] * 100), "\n")
cat("- Specificity:", sprintf("%.1f%%", petti_conf$byClass["Specificity"] * 100), "\n")
cat("- AUC:", sprintf("%.3f", petti_auc), "\n\n")
if (length(petti_missing) > 0) {
  cat("NOTE: ", length(petti_missing), " gene(s) missing in PETTI, imputed as zero.\n")
  cat("      This tests the signature's robustness to incomplete gene coverage.\n")
} else {
  cat("NOTE: All 17 genes successfully mapped using alternative names.\n")
  cat("      Gene name mapping: ARMH1 (reference) → C1orf228 (PETTI)\n")
}
cat(paste(rep("=", 70), collapse = ""), "\n")

