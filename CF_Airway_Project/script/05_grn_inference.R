# -------------------------------
# Load Libraries
# -------------------------------
library(GENIE3)
library(Seurat)

# -------------------------------
# Load normalized Seurat object
# -------------------------------
combined <- readRDS("results/combined_normalized.rds")

# Extract expression matrix
expr_mat <- as.matrix(GetAssayData(combined, slot = "data"))

# Optional: subset to highly variable genes
hvgs <- VariableFeatures(combined)
expr_mat_hvg <- expr_mat[hvgs, ]

# -------------------------------
# Run GENIE3
# -------------------------------
set.seed(123)
cat("Running GENIE3, this may take some time...\n")
weight_matrix <- GENIE3(expr_mat_hvg)

# -------------------------------
# Create folders if not exist
# -------------------------------
dir.create("results/tables/GRN_importance_scores", recursive = TRUE, showWarnings = FALSE)

# Save importance scores
saveRDS(weight_matrix, file = "results/tables/GRN_importance_scores/GRN_weights.rds")
cat("GENIE3 importance scores saved!\n")