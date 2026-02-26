# cf-grn-singlecell-ai

**AI-Driven Reconstruction of Gene Regulatory Networks in Cystic Fibrosis Airway Epithelium Using Single-Cell RNA Sequencing**

---

 ## Project Overview

This project develops an AI-driven framework to reconstruct gene regulatory networks (GRNs) from single-cell RNA sequencing (scRNA-seq) data of airway epithelial cells in cystic fibrosis (CF).

Using machine learning–based network inference, this study aims to identify key transcriptional regulators driving CF airway pathophysiology.

---
Dataset: GSE271984
Source: NCBI Gene Expression Omnibus (GEO)
Tissue Focus: Bronchoalveolar lavage (BAL), Bronchial, and Tracheal airway cells
Disease Context: Cystic Fibrosis (CF)

---

## Objectives

-Process and integrate single-cell RNA-seq datasets from CF and healthy controls

-Perform dimensionality reduction and clustering of airway epithelial cells

-Reconstruct gene regulatory networks using AI-based methods

-Identify high-centrality transcription factors driving CF-specific phenotypes

Provide a reproducible computational pipeline for GRN inference



---

##  Dataset Overview

| Sample Type | Number of Samples | Tissue |
|-------------|-----------------|--------|
| CF Patients | 6 | BAL, Bronchial, Tracheal |
| Healthy Controls | 7 | BAL, Bronchial, Tracheal |

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

## 📌 Disclaimer

This is an **ongoing research project**. Figures, analyses, and interpretations may evolve as additional AI-based GRN work is completed.

