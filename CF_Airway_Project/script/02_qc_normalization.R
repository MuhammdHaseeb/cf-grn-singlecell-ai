# -------------------------------
# Load Libraries
# -------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# -------------------------------
# Load merged Seurat object
# -------------------------------
combined <- readRDS("results/combined_raw.rds")

# -------------------------------
# QC metrics
# -------------------------------
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# -------------------------------
# Create folders if not exist
# -------------------------------
dir.create("results/figures/ViolinPlots", recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Violin plots: nFeature_RNA, nCount_RNA, percent.mt
# -------------------------------
vln <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("results/figures/ViolinPlots/QC_violin_plot.png", vln, width = 10, height = 5)

# -------------------------------
# Filter cells
# -------------------------------
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

# -------------------------------
# Normalize data
# -------------------------------
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# Save processed object
saveRDS(combined, file = "results/combined_normalized.rds")