library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/hollySparks/dataIntegration_wSeurat/horse_and_mouse')
dir <- file.path(getwd(), "quant")

# Read the matrix.mtx file
counts_matrix <- readMM(file = file.path(dir, "matrix.mtx"))

# Read the genes.tsv file
genes <- read.table(file = file.path(dir, "genes.tsv"), header = FALSE, sep = "\t")
gene_names <- genes$V1

# Read the barcodes.tsv file
barcodes <- read.table(file = file.path(dir, "barcodes.tsv"), header = FALSE, sep = "\t")
barcode_names <- barcodes$V1

# Assign row and column names to the counts matrix
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- barcode_names

# Create a new Seurat object
mouse_seurat <- CreateSeuratObject(counts = counts_matrix)

horse <- readRDS("horse.rds")
horse <- JoinLayers(horse)

obj <- merge(mouse_seurat, y = c(horse), add.cell.ids = c("mouse", "horse"), 
             project = "obj_merged")

# export object to be able to use it in ARC
# saveRDS(obj.merged, "mouse-horse_merged_seurat_object.rds")

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
png("UMAPplot_by_experiment_unintegrated.png")

DimPlot(obj, reduction = "umap.unintegrated",
        group.by = "orig.ident")
dev.off()


# CCA integration
 obj <- IntegrateLayers(
   object = obj, method = CCAIntegration,
   orig.reduction = "pca", new.reduction = "integrated.cca",
   dims = 1:10, # Reduce the maximum dimensions to 10
   verbose = FALSE
 )

# RPCA integration
 options(future.globals.maxSize = 5 * 1024^3)
 obj <- IntegrateLayers(
   object = obj, method = RPCAIntegration,
   orig.reduction = "pca", new.reduction = "integrated.rpca",
   dims = 1:10, # Reduce the maximum dimensions to 10
   k.anchor = 30, # Increasing k.anchor
   verbose = FALSE
 )

# # Harmony integration
obj <- IntegrateLayers(
   object = obj, method = HarmonyIntegration,
   orig.reduction = "pca", new.reduction = "harmony",
   verbose = FALSE
 )

# FastMNN integration 
library(SeuratWrappers)
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

### Create plots

plotIntegratedData <- function(obj, keyword_string){
  if (keyword_string == "harmony"){
    red_method <- keyword_string
  } else {
    red_method <- paste0("integrated.", keyword_string)
  }
  
  clust_name <- paste0(keyword_string, "_clusters")
  red_name   <- paste0("umap.", keyword_string)
  
  obj <- FindNeighbors(obj, reduction = red_method, dims = 1:30)
  obj <- FindClusters(obj, resolution = 2, cluster.name = clust_name)
  
  obj <- RunUMAP(obj, reduction = red_method, dims = 1:30, reduction.name = red_name)
  
  p1 <- DimPlot(obj, reduction = red_name, group.by = c(clust_name),combine = FALSE, label.size = 2)
  p2 <- DimPlot(obj, reduction = red_name, group.by = c("orig.ident"),combine = FALSE, label.size = 2)
  
  combined_plot <- p1[[1]] + p2[[1]]
  
  umapplot_filename <- paste0(keyword_string, "_UMAPplot.png")
  png(umapplot_filename, width = 1600, height = 800)
  print(combined_plot)
  dev.off()
  cat("UMAP plot was saved to file: ", umapplot_filename)
}

plotIntegratedData(obj, "cca")
plotIntegratedData(obj, "rpca")
plotIntegratedData(obj, "harmony")
plotIntegratedData(obj, "mnn")
plotIntegratedData(obj, "scvi")

saveRDS(obj, file = "int_mouse_horse_seurat_obj.rds")


scvi_obj <- readRDS("scvi_reduction.rds")
scvi_obj

