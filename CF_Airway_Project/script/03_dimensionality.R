# -------------------------------
# Load Libraries
# -------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# -------------------------------
# Load normalized Seurat object
# -------------------------------
combined <- readRDS("results/combined_normalized.rds")

# -------------------------------
# Scale data
# -------------------------------
combined <- ScaleData(combined, features = rownames(combined))

# -------------------------------
# PCA
# -------------------------------
combined <- RunPCA(combined, features = VariableFeatures(combined))
ElbowPlot(combined)
ggsave("results/figures/UMAPs/PCA_ElbowPlot.png")

# -------------------------------
# Clustering
# -------------------------------
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

# -------------------------------
# UMAP visualization
# -------------------------------
combined <- RunUMAP(combined, dims = 1:20)

# Create folders
dir.create("results/figures/UMAPs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/FeaturePlots", recursive = TRUE, showWarnings = FALSE)

# Save UMAPs
umap_sample <- DimPlot(combined, reduction = "umap", group.by = "SampleID") + ggtitle("UMAP by Sample")
umap_patient <- DimPlot(combined, reduction = "umap", group.by = "Patient") + ggtitle("UMAP by Patient")
ggsave("results/figures/UMAPs/UMAP_by_Sample.png", umap_sample, width = 7, height = 5)
ggsave("results/figures/UMAPs/UMAP_by_Patient.png", umap_patient, width = 7, height = 5)

# FeaturePlots: top 6 variable genes
top6 <- head(VariableFeatures(combined), 6)
for (gene in top6) {
  fp <- FeaturePlot(combined, features = gene)
  ggsave(paste0("results/figures/FeaturePlots/FeaturePlot_", gene, ".png"), fp, width = 5, height = 4)
}

# Save object
saveRDS(combined, file = "results/combined_dimensionality.rds")