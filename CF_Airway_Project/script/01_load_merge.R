# -------------------------------
# Load Libraries
# -------------------------------
library(Seurat)
library(dplyr)
library(Matrix)

# -------------------------------
# Set working directory
# -------------------------------
setwd("C:/Users/muham/OneDrive/Desktop/CF_Airway_Project")

# -------------------------------
# Identify all sample folders
# -------------------------------
sample_paths <- list.files("raw_data", full.names = TRUE)
sample_names <- basename(sample_paths)

# -------------------------------
# Create Seurat objects for each sample
# -------------------------------
seurat_list <- list()

for (i in seq_along(sample_paths)) {
  path <- sample_paths[i]
  name <- sample_names[i]
  
  cat("Processing sample:", name, "\n")
  
  # Read 10X data
  counts <- Read10X(data.dir = path)
  
  # Create Seurat object
  obj <- CreateSeuratObject(counts = counts, project = name, min.cells = 3, min.features = 200)
  
  # Add sample info
  obj$SampleID <- name
  obj$Patient <- sapply(strsplit(name, "_"), `[`, 1)
  obj$Tissue <- sapply(strsplit(name, "_"), `[`, 2)
  
  seurat_list[[name]] <- obj
}

# -------------------------------
# Merge all Seurat objects
# -------------------------------
combined <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)

cat("Merged Seurat object contains", ncol(combined), "cells and", nrow(combined), "genes.\n")

# Save RDS
saveRDS(combined, file = "results/combined_raw.rds")