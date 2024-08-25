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

options(future.globals.maxSize = 1e9)

# Because the precious command did not work, the database was imported 
# in the local computer and copied here.

setwd('/Users/juanjovel/jj/data_analysis/hollySparks/dataIntegration_wSeurat/horse_and_mouse')
       
# Read/import RDS objects
mouse  <- readRDS("mouse.rds")
horse   <- readRDS("horse.rds")

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

# Display the result
cat("Number of converted mouse genes to human genes:", sum(conversion_result$Human_Gene != "UNMATCHED"), "\n")
cat("Number of mouse genes with no human match:", sum(is.na(conversion_result$Human_Gene)), "\n")

# Ensure final_gene_list has no duplicates and is in the correct order
final_gene_list <- make.unique(final_gene_list)

# Rename the rows of the counts matrix
rownames(mouse@assays[["RNA"]]$counts) <- final_gene_list


obj <- merge(mouse, y = c(horse), add.cell.ids = c("mouse", "horse"), 
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

