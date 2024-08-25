# https://satijalab.org/seurat/articles/seurat5_integration
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(reticulate)

options(future.globals.maxSize = 1e9)

# Because the precious command did not work, the database was imported 
# in the local computer and copied here.

setwd('/Users/juanjovel/jj/data_analysis/juanJovel/pipelines/singleCell/dataIntegration_wSeurat/horse_and_mouse')

obj <- readRDS("mouse-horse_merged_seurat_object.rds")

# select only cells with at least 1000 features
obj <- subset(obj, nFeature_RNA > 1000)

obj <- JoinLayers(obj, assay = NULL, layers = NULL, new = NULL)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# Visualize a standard analysis without integration
# While a UMAP analysis is just a visualization of this, clustering 
# this dataset would return predominantly batch-specific clusters. 
# Especially if previous cell-type annotations were not available, 
# this would make downstream analysis extremely challenging.
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Visualize by batch and cell type annotation
# cell type annotation were previously added by Azimuth
png("UMAPplot_by_method.png")

DimPlot(obj, reduction = "umap.unintegrated",
        group.by = "orig.ident")
dev.off()

### Try another normalization step here


# CCA integration
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  dims = 1:10, # Reduce the maximum dimensions to 10
  verbose = FALSE
)

# RPCA integration
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  dims = 1:10, # Reduce the maximum dimensions to 10
  k.anchor = 30, # Increasing k.anchor
  verbose = FALSE
)

# Harmony integration
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# FastMNN integration 
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

# scVI integration 
# https://scvi-tools.org/
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "~/miniforge3/envs/scvi", verbose = FALSE
)

##########################################################
#         Generate figures and markdown report           #
##########################################################

---
title: "Merging data from different tendon RNAseq experiments"
author: "Juan Jovel"
date: "2024-05-07"
output: 
  html_document:
  fig_width: 8
fig_height: 6
---
  
# Introduction
  
Here, several Seurat scripts have been used to compile this script. The contant 
is mainly from 



# CCA_____________________
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:10)
obj <- FindClusters(obj, resolution = 2, cluster.name = "cca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:10, reduction.name = "umap.cca")

p1 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("Method", "predicted.celltype.l2", "cca_clusters"),
  combine = FALSE, label.size = 2
)

png("UMAPplot_CCAintegration.png")
p1
dev.off()

# RPCA_____________________
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:10)
obj <- FindClusters(obj, resolution = 2, cluster.name = "rpca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:10, reduction.name = "umap.rpca")
p2 <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = c("orig.ident", "rpca_clusters"),
  combine = FALSE, label.size = 2
)

png("UMAPplot_CCAintegration.png")
p2
dev.off()


# Harmony__________________

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p2 <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = c("orig_identity", "harmony_clusters"),
  combine = FALSE, label.size = 2
)

png("UMAPplot_harmonyintegration.png")
p2
dev.off()


# FastMNN__________________

obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "mnn_clusters")

obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
p2 <- DimPlot(
  obj,
  reduction = "umap.mnn",
  group.by = c("orig_identity", "mnn_clusters"),
  combine = FALSE, label.size = 2
)

png("UMAPplot_SCVIintegration.png")
p2
dev.off()



# SCVI__________________

obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "scvi_clusters")

obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
p2 <- DimPlot(
  obj,
  reduction = "umap.scvi",
  group.by = c("orig_identity", "scvi_clusters"),
  combine = FALSE, label.size = 2
)

png("UMAPplot_SCVIintegration.png")
p2
dev.off()



#wrap_plots(c(p1, p2), ncol = 2, byrow = F)


# compare the expression of biological markers based on 
# different clustering solutions, or visualize one method’s 
# clustering solution on different UMAP visualizations.
p1 <- VlnPlot(
  obj,
  features = "rna_CD8A", group.by = "unintegrated_clusters"
) + NoLegend() + ggtitle("CD8A - Unintegrated Clusters")
p2 <- VlnPlot(
  obj, "rna_CD8A",
  group.by = "cca_clusters"
) + NoLegend() + ggtitle("CD8A - CCA Clusters")
p3 <- VlnPlot(
  obj, "rna_CD8A",
  group.by = "scvi_clusters"
) + NoLegend() + ggtitle("CD8A - scVI Clusters")

png("violinPlot_CD8A_unIntegrated.png")
p1
dev.off()

png("violinPlot_CD8A_CCAintegration.png")
p2
dev.off()

png("violinPlot_CD8A_scVIintegration.png")
p3
dev.off()

# UMAP clustering
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p4 <- DimPlot(obj, reduction = "umap.unintegrated", group.by = c("cca_clusters"))
p5 <- DimPlot(obj, reduction = "umap.rpca", group.by = c("cca_clusters"))
p6 <- DimPlot(obj, reduction = "umap.scvi", group.by = c("cca_clusters"))

png("UMAPplot_unintegrated.png")
p4
dev.off()

png("UMAPplot_RPCAintegration.png")
p5
dev.off()

png("UMAPplot_scVIintegration.png")
p6
dev.off()

# Once integrative analysis is complete, you can rejoin the layers
obj <- JoinLayers(obj)
obj

# Lastly, users can also perform integration using sctransform-normalized data
# https://htmlpreview.github.io/?https://github.com/satijalab/sctransform/blob/supp_html/supplement/seurat.html
options(future.globals.maxSize = 3e+09)
obj <- SCTransform(obj)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)

# obj <- FindNeighbors(obj, dims = 1:30, reduction = "integrated.dr")
# that reduction method was not found by the interpreter
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2)

# obj <- RunUMAP(obj, dims = 1:30, reduction = "integrated.dr")
# that reduction method was not found by the interpreter
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca")


DimPlot(obj, reduction = "pca", 
        group.by = "predicted.celltype.l2")
