# https://satijalab.org/seurat/articles/seurat5_integration
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(reticulate)
library(biomaRt)
library(tidyverse)
library(Matrix)

options(future.globals.maxSize = 1e9)

# Because the precious command did not work, the database was imported 
# in the local computer and copied here.

setwd('/Users/juanjovel/jj/data_analysis/hollySparks/dataIntegration_wSeurat/human_and_mouse')
       
# Read/import RDS objects
mouse  <- readRDS("mouse.rds")
human_d <- readRDS("obj_kendal_human_diseased.rds")
human_h <- readRDS("obj_kendal_human_healthy.rds")


##################
# Extract count layers
count_layers <- grep("counts", names(human_d[["RNA"]]@layers), value = TRUE)

# Sum the count layers to create a consolidated matrix
consolidated_counts <- Reduce("+", lapply(count_layers, function(layer) human_d[["RNA"]]@layers[[layer]]))

# Extract gene names and barcodes
human_genes <- rownames(consolidated_counts)
human_barcodes <- colnames(consolidated_counts)

##################




# Extract the mouse gene symbols
mouse_genes <- rownames(mouse[["RNA"]]$counts)

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

# write.csv(mouse_human_genes, 'mouse_human_genes.rpt')

# Alternatively
#mouse_human_genes <- read.csv('mouse_human_genes.rpt')

# Create a mapping dictionary
mouse_to_human <- mouse_human_genes %>%
  filter(Common.Organism.Name == "mouse, laboratory") %>%
  select(DB.Class.Key, Symbol) %>%
  rename(Mouse_Symbol = Symbol)

human_genes <- mouse_human_genes %>%
  filter(Common.Organism.Name == "human") %>%
  select(DB.Class.Key, Symbol) %>%
  rename(Human_Symbol = Symbol)

gene_mapping <- merge(mouse_to_human, human_genes, by = "DB.Class.Key")

convert_mouse_to_human <- function(gene_list, mapping_data) {
  mapping_dict <- setNames(mapping_data$Human_Symbol, mapping_data$Mouse_Symbol)
  converted_genes <- mapping_dict[gene_list]
  
  # Create a data frame to keep the original order
  result_df <- data.frame(
    Mouse_Gene = gene_list,
    Human_Gene = converted_genes,
    stringsAsFactors = FALSE
  )
  
  # Replace NA with Mouse_Gene
  result_df$Final_Gene <- ifelse(is.na(result_df$Human_Gene), result_df$Mouse_Gene, result_df$Human_Gene)
  
  result_df
}

# Convert mouse genes to human genes
conversion_result <- convert_mouse_to_human(mouse_genes, gene_mapping)
final_gene_list <- conversion_result$Final_Gene

# Check for and make unique names in final_gene_list
final_gene_list <- make.unique(final_gene_list)
# Rename the rows of the counts matrix
rownames(mouse@assays[["RNA"]]$counts) <- final_gene_list


obj <- merge(mouse, y = c(human_d, human_h), add.cell.ids = c("mouse", "human_d", "human_h"), 
             project = "obj_merged")

# export object to be able to use it in ARC
# saveRDS(obj.merged, "mouse-horse_merged_seurat_object.rds")

# select only cells with at least 1000 features
obj <- subset(obj, nFeature_RNA > 1000)

obj <- JoinLayers(obj, assay = NULL, layers = NULL, new = NULL)
# Join the already split 'data' layer in the RNA assay
obj[["RNA"]]$data <- JoinLayers(obj[["RNA"]]$data)
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

setwd('/Users/juanjovel/jj/data_analysis/hollySparks/dataIntegration_wSeurat/horse_and_mouse/ensembl_renaming')

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

