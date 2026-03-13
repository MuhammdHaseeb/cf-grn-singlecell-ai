

# **Reconstructed gene regulatory networks in cystic fibrosis airway epithelium using AI-driven analysis of publicly available single-cell RNA-seq dataset**
OR
Cell-Type-Specific Gene Regulatory Network Disruption in Cystic Fibrosis Airway Epithelium: A Machine Learning-Assisted Single-Cell Transcriptomic Analysis

---

 ## Project Overview

This project develops an AI-driven framework to reconstruct gene regulatory networks (GRNs) from single-cell RNA sequencing (scRNA-seq) data of airway epithelial cells in cystic fibrosis (CF).

Using machine learning–based network inference, this study aims to identify key transcriptional regulators driving CF airway pathophysiology.

---
- **Dataset:** GSE271984
- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Tissue Focus:** Bronchoalveolar lavage (BAL), Bronchial, and Tracheal airway cells
- **Disease Context:** Cystic Fibrosis (CF)

---

## Objectives

- Process and integrate single-cell RNA-seq datasets from CF and healthy controls

- Perform dimensionality reduction and clustering of airway epithelial cells

- Reconstruct gene regulatory networks using AI-based methods

- Identify high-centrality transcription factors driving CF-specific phenotypes

- Provide a reproducible computational pipeline for GRN inference



---

##  Dataset Overview

| Sample Type | Number of Samples | Tissue |
|-------------|-----------------|--------|
| CF Patients | 3 | BAL, Bronchial, Tracheal |
| Healthy Controls | 4 | BAL, Bronchial, Tracheal |

# Methodology Pipeline

---

## 1️⃣ Data Acquisition

- Downloaded raw expression matrices from **GEO (GSE271984)**
- Selected airway samples:
  - BAL (Bronchoalveolar Lavage)
  - Bronchial
  - Tracheal
- Created structured metadata including:
  - Disease condition (CF vs Control)
  - Tissue origin

---

## 2️⃣ Preprocessing (Seurat – R)

Performed single-cell RNA-seq preprocessing using **Seurat**:

- Quality control filtering
- Removed low-quality cells
- Filtered cells with high mitochondrial gene expression
- Log-normalization of gene expression
- Identification of top 2000 highly variable genes
- Data scaling
- Dimensionality reduction

---

## 3️⃣ Exploratory Analysis

- Principal Component Analysis (PCA)
- UMAP visualization
- Graph-based cell clustering
- Annotation of epithelial subtypes:
  - Basal cells
  - Secretory cells
  - Ciliated cells

---

## 4️⃣ AI-Based Gene Regulatory Network (GRN) Inference

- Applied **GENIE3** to normalized expression matrix
- Used curated transcription factor list as candidate regulators
- Generated weighted gene regulatory network
- Ranked regulator importance scores

---

## 5️⃣ Network Analysis

- Constructed regulatory network graph
- Computed centrality metrics:
  - Degree centrality
  - Betweenness centrality
- Identified potential regulatory hubs in CF airway epithelial cells
- Network visualization using:
  - R (igraph / ggraph)
  - Cytoscape

---

## 📁 Folder Structure

- CF_Airway_Project/
- │
- ├── raw_data/
- │   │
- │   ├── CF1_BAL/
- │   │   ├── barcodes.tsv.gz
- │   │   ├── features.tsv.gz
- │   │   └── matrix.mtx.gz
- │   │
- │   ├── CF1_Bronchial/
- │   ├── CF1_Tracheal/
- │   │
- │   ├── CF2_BAL/
- │   ├── CF2_Bronchial/
- │   ├── CF2_Tracheal/
- │   │
- │   ├── CO1_BAL/
- │   ├── CO1_Bronchial/
- │   ├── CO1_Tracheal/
- │   │
- │   ├── CO2_BAL/
- │   ├── CO2_Bronchial/
- │   ├── CO2_Tracheal/
- │   │
- │   └── CO3_BAL/
- │
- ├── scripts/
- │   │
- │   ├── 00_install_packages.R
- │   ├── 01_load_and_merge.R        # Load raw data, create Seurat objects, merge, initial QC
- │   ├── 02_qc_normalization.R      # QC filtering, normalization, mitochondrial percentage
- │   ├── 03_dimensionality.R        # PCA, scaling, clustering, UMAP
- │   ├── 04_DE_analysis.R           # Differential expression (CF vs Control)
- │   ├── 05_grn_inference.R         # GENIE3 network inference
- │   └── 06_network_analysis.R      # Network visualization & hub TF identification
- │
- ├── metadata/
- │   └── sample_info.csv
- │
- ├── results/
- │   │
- │   ├── figures/
- │   │   ├── UMAPs/
- │   │   ├── ViolinPlots/
- │   │   ├── FeaturePlots/
- │   │   └── GRN_networks/
- │   │
- │   └── tables/
- │       ├── DE_results/
- │       └── GRN_importance_scores 




## 📌 Disclaimer

This is an **ongoing research project**. Figures, analyses, and interpretations may evolve as additional AI-based GRN work is completed.

