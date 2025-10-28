# Combining Projection and Supervised Learning to Detect Malignant Cells in AML scRNA-Seq Datasets

This repository contains the implementation and analysis scripts for a project comparing multiple computational approaches for
detecting malignant cells in Acute Myeloid Leukemia (AML) single-cell
RNA sequencing datasets.


## Overview

This project evaluates and compares four distinct computational methods
for distinguishing malignant from healthy cells in AML scRNA-seq data:

1.  **PCA-Based Projection** - Unsupervised projection onto healthy
    reference space
2.  **SCMAP-Based Projection** - Reference-based cell type projection
    (two variants)
3.  **LASSO Regression** - Supervised feature selection and
    classification
4.  **17-Gene Signature SVM** - External benchmark from Petti et al.
    (2019)

## Repository Structure

```         
.
├── clonetracer_data_predict_validate.Rmd        # PCA-based projection method
├── run_scmap_projection_optimized.R             # SCMAP cluster projection
├── run_scmap_with_majority_voting_OPTIMIZED.R   # SCMAP cell-level majority voting
├── LASSO_analysis.R                             # LASSO regression analysis
├── SVM_17gene_with_imputation.R                 # 17-gene signature SVM
├── patient_misclassification_analysis.R         # Patient-level performance analysis
├── sanity_check_malignant_immune_cells.R        # Validation of immune cell annotations
└── verify_non_immune_analysis.R                 # Non-immune cell subset validation
```

------------------------------------------------------------------------

## Main Analysis Scripts

### 1. PCA-Based Projection Method

**Script:** `clonetracer_data_predict_validate.Rmd`

Implements projection-based malignant cell detection using Principal
Component Analysis: 
- Projects test data onto healthy reference PCA
space (1,991 features)
- Uses weighted Mahalanobis distance for
classification
- Evaluates performance by cell type
- Validates predictions against lineage-tracing ground truth

**Key Outputs:** 
- Sensitivity and specificity per cell type
- Per-patient performance metrics

------------------------------------------------------------------------

### 2. SCMAP-Based Projection Methods

#### Script A: `run_scmap_projection_optimized.R`

Implements **scmapCluster()** projection method: 
- Cluster-level projection to healthy reference
- Uses same 1,991 features as PCA for fair comparison
- Log-normalized expression matching
- Cell ID verification for traceability

#### Script B: `run_scmap_with_majority_voting_OPTIMIZED.R`

Implements **scmapCell()** with majority voting: 
- Cell-level nearest neighbor matching
- Majority voting across k-nearest neighbors
- Memory-optimized batch processing
- Handles large datasets efficiently

**Key Features:** 
- Both methods use SingleCellExperiment framework
- Optimized for memory efficiency
- Comparable feature sets across methods

------------------------------------------------------------------------

### 3. LASSO Regression

**Script:** `LASSO_analysis.R`

Supervised machine learning approach with automated feature selection: 
- Cross-validation for optimal lambda selection
- Elastic net regularization (α = 1 for LASSO)
- Runs on full dataset and non-immune cells separately
- Parallel processing for computational efficiency

**Key Outputs:** 
- Selected gene features
- Cross-validated performance metrics
- ROC curves and AUC values
- Confusion matrices

**Dependencies:** `glmnet`, `caret`, `pROC`, `doParallel`

------------------------------------------------------------------------

### 4. 17-Gene Signature SVM

**Script:** `SVM_17gene_with_imputation.R`

External benchmark implementing the published AML gene signature from:
\> Petti et al. (2019) "A general approach for detecting expressed
mutations in AML cells using single cell RNA-sequencing"

**Features:** 
- Trains radial kernel SVM on 17-gene signature
- Tests on healthy reference (all genes available)
- Tests on external PETTI dataset with missing gene imputation
- Zero imputation for missing genes (ARMH1 in PETTI data)

**17-Gene Signature:** MPO, ELANE, AZU1, PRTN3, CTSG, LYZ, ANXA1,
S100A8, S100A9, CD52, LGALS1, LAPTM5, FCER1G, TYROBP, HLA-DRA, CST3,
CYBA

------------------------------------------------------------------------

## Supporting Analysis Scripts

### Patient-Level Analysis

**Script:** `patient_misclassification_analysis.R`

Analyzes model performance at the patient level: 
- Patient-specific accuracy, sensitivity, specificity
- Misclassification pattern visualization by patient and cell type
- False positive/negative rate analysis
- Identifies patient-specific factors affecting performance

### Validation Scripts

**Script:** `sanity_check_malignant_immune_cells.R`

Quick validation to identify malignant cells in immune populations: 
- Checks for malignant T-cells and B-cells
- Explains sensitivity differences in immune vs. non-immune subsets
- Validates ground truth annotations from lineage tracing

**Script:** `verify_non_immune_analysis.R`

Validates performance on non-immune cell subsets specifically.

------------------------------------------------------------------------

## Requirements

### R Version

-   R ≥ 4.0.0

### Required R Packages

**Core Analysis:** - `Seurat` (≥ 4.0) - `dplyr`, `tidyverse` - `Matrix`

**Method-Specific:** - `scmap` - for SCMAP projection methods -
`SingleCellExperiment` - for SCMAP framework - `glmnet` - for LASSO
regression - `e1071` - for SVM implementation - `caret` - for
cross-validation - `pROC` - for ROC curve analysis

**Visualization:** - `ggplot2` - `gridExtra` - `ComplexHeatmap`
(optional)

**Performance:** - `parallel`, `doParallel` - for parallel processing

### Installation

``` r
# Install from CRAN
install.packages(c("Seurat", "dplyr", "Matrix", "ggplot2", 
                   "glmnet", "e1071", "caret", "pROC",
                   "gridExtra", "doParallel"))

# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("scmap", "SingleCellExperiment"))
```

------------------------------------------------------------------------

## Data Requirements

### Input Data Format

1.  **Test Dataset** (Clone Tracer AML data):
    -   Seurat object with malignant and healthy cells
    -   Ground truth labels from lineage tracing
    -   Cell type annotations
    -   Patient identifiers
2.  **Healthy Reference Dataset**:
    -   Seurat object with confirmed healthy cells only
    -   Same cell type annotations as test data
    -   Matching gene features (or imputation strategy)
3.  **External Validation** (PETTI dataset):
    -   For SVM benchmark validation
    -   Contains independent AML samples

### Required Metadata Fields

-   `status`: "leukemic" or "healthy" (ground truth)
-   `ct` or `celltype`: Cell type annotations
-   `patient` or `sample`: Patient identifier
-   `cell_id`: Unique cell identifier

------------------------------------------------------------------------

## Usage

### 1. Run PCA-Based Projection

``` r
# Open in RStudio and knit
rmarkdown::render("clonetracer_data_predict_validate.Rmd")
```

### 2. Run SCMAP Methods

``` bash
# SCMAP cluster projection
Rscript run_scmap_projection_optimized.R

# SCMAP cell-level with majority voting
Rscript run_scmap_with_majority_voting_OPTIMIZED.R
```

### 3. Run LASSO Regression

``` bash
Rscript LASSO_analysis.R
```

### 4. Run 17-Gene SVM

``` bash
Rscript SVM_17gene_with_imputation.R
```

### 5. Run Error Analysis

``` bash
# Patient-level analysis
Rscript patient_misclassification_analysis.R

# Validation checks
Rscript sanity_check_malignant_immune_cells.R
```

------------------------------------------------------------------------

## Output Files

Each script generates timestamped output files:

-   `.RDS` files: R data objects for downstream analysis
-   `.csv` files: Tabular results and predictions
-   `.png`/`.pdf` files: Visualization plots
-   `.txt` files: Summary statistics and logs

Example outputs:

```         
pca_validation_results_corrected_20250918_121705.RDS
scmap_majority_voting_OPTIMIZED_predictions_20250923_143355.csv
lasso_regression_summary_20250902_160849.csv
svm_17gene_with_imputation_results_20251002_115633.RDS
```

------------------------------------------------------------------------

## Performance Summary

Typical performance metrics across methods (from thesis results):

| Method              | Sensitivity | Specificity | AUC    |
|---------------------|-------------|-------------|--------|
| PCA Projection      | \~94%       | \~95%       | \~0.95 |
| SCMAP Cluster       | \~92%       | \~93%       | \~0.93 |
| SCMAP Cell Majority | \~91%       | \~94%       | \~0.92 |
| LASSO Regression    | \~89%       | \~91%       | \~0.90 |
| 17-Gene SVM         | \~85%       | \~88%       | \~0.87 |

*Note: Exact values depend on dataset and validation strategy*

------------------------------------------------------------------------

## Key Findings

1.  **PCA projection** achieves highest overall performance with
    unsupervised approach
2.  **SCMAP methods** provide competitive performance with explicit
    reference mapping
3.  **LASSO regression** offers interpretability through feature
    selection
4.  **17-gene SVM** demonstrates transferability but lower performance
    on this dataset
5.  **Cell type matters**: Performance varies significantly by cell type
    (HSCs, monocytes, etc.)
6.  **Patient variability**: Some patients show systematic
    misclassification patterns

------------------------------------------------------------------------

## Citation

If you use this code or approach in your research, please cite:

```         
Parastoo Jargoie(2025). Combining Projection and Supervised Learning to Detect 
Malignant Cells in AML scRNA-Seq Datasets. Master's Thesis, Potsdam University.
```

------------------------------------------------------------------------

## License

This project is part of a master's thesis. Please contact the author for
usage permissions.

------------------------------------------------------------------------

## Contact

For questions or collaboration: 
- **Author:** Parastoo Jargoie
- **Institution:** Potsdam University
- **Email:** Jargoiep@yahoo.com
- **Thesis Supervisor:** Prof. Dr. Zoran Nikoloski, MUDr. Jan Bařinka, Prof. Dr. Simon Haas

------------------------------------------------------------------------

## Acknowledgments

-   MUDr. Jan Bařinka for the PCA projection framework
-   Haas lab for computational resources
-   Original data providers and method developers

------------------------------------------------------------------------

**Last Updated:** October 2025
