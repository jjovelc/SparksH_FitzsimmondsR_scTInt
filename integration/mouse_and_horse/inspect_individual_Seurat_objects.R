library(tidyverse)
library(Seurat)
library(patchwork)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
library(BPCells)

setwd('/Users/juanjovel/jj/data_analysis/juanJovel/pipelines/singleCell/dataIntegration_wSeurat/horse_and_mouse')
dir <- getwd()


obj <- readRDS("obj_zamboulis_horse_young_old.rds")
nrow(obj@meta.data) # 2700 cells



# Pre-processing of data some QC procedures, data normalization, scalling and detection 
# of highly-variable features

# Relevant QC features include

# 1. The number of unique genes detected in each cell.
# 2. Total number of molecules detected within a cell
# 3. The percentage of reads that map to the mitochondrial genome

# So we can store the QC report in the metadata table
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Show QC metrics for the first 10 cells
head(obj@meta.data, 10)

# QC metrics can be visualized in a violin plot
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        col=c('dodgerblue', 'firebrick'))


# Join or merge the two layers within the mouse dataset
# Inspect the layers within the RNA assay before merging


# Merge the layers within the RNA assay
obj[["RNA"]] <- merge(x = obj[["RNA"]], y = obj[["RNA"]], merge.data = TRUE)




# FeatureScatter is typically used to visualize feature-feature relationships
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                                          )
plot2

# The violin plots would suggest to filter out data that contains less than 200 
# features (detected genes) or more than 2000 features. 
# HOW WE AUTOMATE THIS?
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# Join or merge the two layers within the mouse dataset


# Let's visualize again after filtering of data
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Normalizing the data
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes 
# the feature expression measurements for each cell by the total expression, multiplies this 
# by a scale factor (10,000 by default), and log-transforms the result.

# In Seurat v5, Normalized values are stored in pbmc[["RNA"]]$data.
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Normalization and scaling should not affect data distribution
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
# Other methods of normalization are available. See them with ?NormalizeData

# Identification of highly variable features (feature selection)

# The mean-variance relationship is analysed with the FindVariableFeatures() function. 
# By default, 2,000 features are returned per dataset. These will be used in 
# downstream analysis, like PCA.
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0,
                     ynudge = 0)
plot1 + plot2

# Scaling the data to have mean 0 and variance 1, to reduce effect of 
# highly expressing genes
# The results of this are stored in pbmc[["RNA"]]$scale.data

# By default, only variable features are scaled.
obj <- ScaleData(obj)

# The features argument can be used to scale additional features 
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)

# ScaleData() function can also be used to remove unwanted sources,
# for example that fraction of variance associated with the amount
# of mitochondrial contamination.


# Linear dimensional reduction
# Next we perform PCA on the scaled data. By default, only the previously 
# determined variable features are used as input

# For the first principal components, Seurat outputs a list of genes with the 
# most positive and negative loadings
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

# Examine and visualize PCA results a few different ways
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(obj, dims = 1:2, reduction = "pca")

DimPlot(obj, reduction = "pca") + NoLegend()

# Though clearly a supervised analysis, we find this to be a valuable tool for 
# exploring correlated feature sets.
DimHeatmap(obj, dims = 1:2, cells = 500, balanced = TRUE)

# Let's now do it for the first 15 dimensions
DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(obj)

# The FindClusters() function implements this procedure, and contains a resolution 
# parameter that sets the ‘granularity’ of the downstream clustering, with increased 
# values leading to a greater number of clusters.

# We find that setting this parameter between 0.4-1.2 typically returns good results 
# for single-cell datasets of around 3K cells. Optimal resolution often increases 
# for larger datasets. The clusters can be found using the Idents() function.
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(obj), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# The goal of these algorithms is to learn underlying structure in the dataset, 
# in order to place similar cells together in low-dimensional space. 
obj <- RunUMAP(obj, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(obj, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)
# FindAllMarkers() automates this process for all clusters, but you can 
# also test groups of clusters vs. each other, or against all cells.

obj <- JoinLayers(obj, assay = NULL, layers = NULL, new = NULL)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(obj, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(obj, ident.1 = 5, ident.2 = c(0, 3))

# find markers for every cluster compared to all remaining cells, report 
# only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Differential expression
# Wilcoxon test run by default
cluster0.markers <- FindMarkers(obj, ident.1 = 0, logfc.threshold = 0.25, 
                                only.pos = TRUE)


cluster0.markers <- FindMarkers(obj, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)

VlnPlot(obj, features = c("COMP", "LOX", "THBS4"))

# you can plot raw counts as well
#VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
#                               "FCGR3A", "LYZ", "PPBP", "CD8A"))

#FeaturePlot(pbmc, features = c("MS4A1", "GNLY"))

# DoHeatmap() generates an expression heatmap for given cells and features. 
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(obj, features = top10$gene) + NoLegend()

# Assigning cell type identity to clusters
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", 
#                     "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet", "others")

#names(new.cluster.ids) <- levels(obj)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)
#DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 4.5) + 
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

print(plot)

saveRDS(obj, file = "horse.rds")
