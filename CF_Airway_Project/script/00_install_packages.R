# -------------------------------
# Install required packages
# -------------------------------
packages <- c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "DoubletFinder", "harmony",
  "cowplot", "pheatmap", "GENIE3", "Matrix"
)

install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, install_if_missing)

# GENIE3 is on Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("GENIE3", quietly = TRUE)){
  BiocManager::install("GENIE3")
}

cat("All packages are installed and loaded.\n")