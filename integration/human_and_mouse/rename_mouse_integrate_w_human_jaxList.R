library(Seurat)
library(tidyverse)
library(reshape2)
library(Matrix)
library(patchwork)
library(SeuratWrappers)

setwd("/Users/juanjovel/OneDrive/jj/UofC/data_analysis/hollySparks/dataIntegration_wSeurat/human_and_mouse")

mouse <- readRDS("mouse.rds")

mouse <- JoinLayers(mouse)

# Regenerate quantification files
# Extract the RNA assay
rna_assay <- mouse[["RNA"]]

# Extract the counts matrix
counts_matrix <- rna_assay$counts

# Extract barcodes (cell names)
barcodes <- colnames(counts_matrix)

# Define the directory path
dir_path <- "quant"

# Check if the directory exists
if (!dir.exists(dir_path)) {
  # If the directory does not exist, create it
  dir.create(dir_path)
  cat("Directory 'quant' created.\n")
} else {
  cat("Directory 'quant' already exists.\n")
}

# Save the counts_matrix to the 'quant' directory
writeMM(counts_matrix, file = file.path(dir_path, "matrix.mtx"))

# Write barcodes to barcodes.tsv
write.table(data.frame(barcodes), file = file.path(dir_path, "barcodes.tsv"), sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Extract the mouse gene symbols
mouse_genes <- rownames(mouse[["RNA"]]$counts)

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

# It is a good idea to save the file locally, instead of downloading it every time
# write.csv(mouse_human_genes, 'mouse_human_genes.rpt')

# If the mouse_human_genes file is saved locally, you can read it in like this
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
final_gene_list  <- conversion_result$Final_Gene

# Ensure final_gene_list has no duplicates and is in the correct order
final_gene_list <- make.unique(final_gene_list)

# Write the gene list to genes.tsv
write.table(data.frame(final_gene_list), file = file.path(dir_path, "genes.tsv"), sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Set up path for quantification files
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

# Import data for horse and human
human.h <- readRDS("obj_kendal_human_healthy.rds")
human.d <- readRDS("obj_kendal_human_diseased.rds")

# Join layers to get a single counts layers
human.h <- JoinLayers(human.h)
human.d <- JoinLayers(human.d)

obj <- merge(mouse_seurat, y = c(human.h, human.d), 
             add.cell.ids = c("mouse", "human.healthy", "human.diseased"), 
             project = "obj_merged")

# export merged object here if you need to
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
png("UMAPplot_mouse-human_unintegrated.png")

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
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

# scVI integration 
# Remember to isntall scVI conda environment, 
# according to instructions in the scVI website
# https://scvi-tools.org/
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "~/miniforge3/envs/scvi", verbose = FALSE
)

### Inspect number of anchors
# Find anchors for CCA
# Preprocess each object before integration
objects <- list(mouse_seurat, human.h, human.d)

for (i in 1:length(objects)) {
  objects[[i]] <- NormalizeData(objects[[i]])
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", nfeatures = 2000)
  objects[[i]] <- ScaleData(objects[[i]])
}

# Find integration anchors using CCA
anchors.cca <- FindIntegrationAnchors(object.list = objects, anchor.features = 2000, reduction = "cca")
num_anchors_cca <- nrow(anchors.cca@anchors)
print(paste("Number of CCA anchors:", num_anchors_cca))

# Find integration anchors using RPCA
anchors.rpca <- FindIntegrationAnchors(object.list = objects, anchor.features = 2000, reduction = "rpca")
num_anchors_rpca <- nrow(anchors.rpca@anchors)
print(paste("Number of RPCA anchors:", num_anchors_rpca))

# Integrate data using the anchors found by CCA
integrated_data_cca <- IntegrateData(anchorset = anchors.cca, dims = 1:30)

### Create plots for each integration method
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
