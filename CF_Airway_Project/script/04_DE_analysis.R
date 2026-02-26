# -------------------------------
# Load Libraries
# -------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# -------------------------------
# Load Seurat object
# -------------------------------
combined <- readRDS("results/combined_dimensionality.rds")

# -------------------------------
# Set CF vs Control metadata
# -------------------------------
combined$Condition <- ifelse(grepl("^CF", combined$SampleID), "CF", "Control")

# -------------------------------
# Differential expression (CF vs Control)
# -------------------------------
de_markers <- FindMarkers(combined, ident.1 = "CF", ident.2 = "Control", group.by = "Condition")

# Save DE table
dir.create("results/tables/DE_results", recursive = TRUE, showWarnings = FALSE)
write.csv(de_markers, file = "results/tables/DE_results/DE_CF_vs_Control.csv")

# -------------------------------
# Create folders for plots
# -------------------------------
dir.create("results/figures/ViolinPlots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/FeaturePlots", recursive = TRUE, showWarnings = FALSE)

# Optional: plot top 10 DE genes
top10 <- head(rownames(de_markers), 10)

# Violin plots
vln <- VlnPlot(combined, features = top10, group.by = "Condition", pt.size = 0.1)
ggsave("results/figures/ViolinPlots/Top10_DE_genes_violin.png", vln, width = 12, height = 6)

# Feature plots
for (gene in top10) {
  fp <- FeaturePlot(combined, features = gene)
  ggsave(paste0("results/figures/FeaturePlots/FeaturePlot_DE_", gene, ".png"), fp, width = 5, height = 4)
}