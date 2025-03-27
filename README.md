# Cross-Species Integration of Tendon scRNA-seq Data

This repository contains code and workflows for integrating single-cell RNA-seq (scRNA-seq) datasets from tendon tissue across three species: **mouse**, **horse**, and **human**. The aim is to uncover conserved cellular populations and investigate interspecies similarities in tendon biology.

## 📁 Repository Structure

- `data/` – Input `.rds` files for mouse, horse, and human scRNA-seq data
- `results/` – Output files including integrated Seurat objects and UMAP plots
- `scripts/` – Main R script for integration and analysis
- `README.md` – This documentation file

## 🧪 Integration Methods Included

- **SCTransform normalization**
- **CCA Integration**
- **RPCA Integration**
- **Harmony**
- **FastMNN**
- **scVI (via scvi-tools & reticulate)**

Each integration method is followed by UMAP visualization and clustering.

## ⚙️ Requirements

### R Packages

Make sure these are installed:

- `Seurat` (v5.0+)
- `SeuratWrappers`
- `tidyverse`
- `reshape2`
- `Matrix`
- `patchwork`
- `reticulate`
- `SeuratDisk` (for Seurat object compatibility and updates)

### Python Environment

For **scVI integration**, a Python environment (via Conda) must be set up with:

- `scvi-tools`
- `scanpy`
- `anndata`

You can activate it in R via `reticulate::use_condaenv()`.

## 🚀 How to Use

1. **Prepare Data**

   Place your `.rds` files in the `data/` directory:
   - `mouse.rds`
   - `horse.rds`
   - `human.rds`

2. **Run Integration**

   From the `scripts/` folder, execute:

   ```r
   Rscript dataIntegration_with_SCTransform.R
