library(Seurat)
library(tidyverse)
library(reshape2)
library(Matrix)
library(patchwork)
library(SeuratWrappers)

# Set your working directory
setwd("/work/vetmed_data/jj/projects/hollySparks/dataIntegration/horse_mouse_human_Feb-2025/seurat_analysis")

# Define a prefix for file naming
prefix <- "horse-human-mouse"

# ---------------------------
# Process Mouse Data
# ---------------------------

# Load mouse data
mouse <- readRDS("mouse.rds")

# Join layers if necessary
mouse <- JoinLayers(mouse)

# Extract the RNA assay and counts matrix
counts_matrix_mouse <- mouse[["RNA"]]$counts
barcodes_mouse <- colnames(counts_matrix_mouse)

# Define the directory path for mouse quantification files
dir_path_mouse <- paste0(prefix, "_mouse_quant")

# Create the directory if it doesn't exist
if (!dir.exists(dir_path_mouse)) {
  dir.create(dir_path_mouse)
  cat("Directory ", dir_path_mouse, " created.\n")
} else {
  cat("Directory ", dir_path_mouse, " already exists.\n")
}

# Save the counts matrix and barcodes for mouse
writeMM(counts_matrix_mouse, file = file.path(dir_path_mouse, "matrix.mtx"))
write.table(data.frame(barcodes_mouse), file = file.path(dir_path_mouse, "barcodes.tsv"), sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Extract the mouse gene symbols
mouse_genes <- rownames(counts_matrix_mouse)

# Read your custom mouse-human ortholog mapping table
# Replace 'human_mouse_symbols_mapping.tsv' with your actual filename
orthologs_mouse <- read.table("human_mouse_symbols_mapping.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create a mapping data frame for mouse
mapping_data_mouse <- orthologs_mouse %>% select(Mouse_Symbol, Human_Symbol)

# Define the conversion function for mouse genes
convert_mouse_to_human <- function(gene_list, mapping_data) {
  # Create a mapping dictionary from Mouse_Symbol to Human_Symbol
  mapping_dict <- setNames(mapping_data$Human_Symbol, mapping_data$Mouse_Symbol)
  
  # Map mouse genes to human genes
  converted_genes <- mapping_dict[gene_list]
  
  # Create a data frame to keep the original order
  result_df <- data.frame(
    Original_Gene = gene_list,
    Human_Gene = converted_genes,
    stringsAsFactors = FALSE
  )
  
  # Replace NA (genes without a human ortholog) with the original gene name
  result_df$Final_Gene <- ifelse(is.na(result_df$Human_Gene), result_df$Original_Gene, result_df$Human_Gene)
  result_df
}

# Convert mouse genes to human genes
conversion_result_mouse <- convert_mouse_to_human(mouse_genes, mapping_data_mouse)
final_gene_list_mouse <- conversion_result_mouse$Final_Gene

# Ensure the final gene list has unique names
final_gene_list_mouse <- make.unique(final_gene_list_mouse)

# Write the gene list to 'genes.tsv' for mouse
write.table(data.frame(final_gene_list_mouse), file = file.path(dir_path_mouse, "genes.tsv"), sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Read the counts matrix and associated files for mouse
counts_matrix_mouse <- readMM(file = file.path(dir_path_mouse, "matrix.mtx"))
genes_mouse <- read.table(file = file.path(dir_path_mouse, "genes.tsv"), header = FALSE, sep = "\t")
gene_names_mouse <- genes_mouse$V1
barcodes_mouse <- read.table(file = file.path(dir_path_mouse, "barcodes.tsv"), header = FALSE, sep = "\t")
barcode_names_mouse <- barcodes_mouse$V1

# Assign row and column names to the counts matrix
rownames(counts_matrix_mouse) <- gene_names_mouse
colnames(counts_matrix_mouse) <- barcode_names_mouse

# Create a Seurat object for mouse with updated gene names
mouse_seurat <- CreateSeuratObject(
  counts = counts_matrix_mouse,
  meta.data = mouse@meta.data  
)


# ---------------------------
# Process Horse Data
# ---------------------------

# Load horse data
horse <- readRDS("horse.rds")

# Join layers if necessary
horse <- JoinLayers(horse)

# Extract the RNA assay and counts matrix
counts_matrix_horse <- horse[["RNA"]]$counts
barcodes_horse <- colnames(counts_matrix_horse)

# Define the directory path for horse quantification files
dir_path_horse <- paste0(prefix, "_horse_quant")

# Create the directory if it doesn't exist
if (!dir.exists(dir_path_horse)) {
  dir.create(dir_path_horse)
  cat("Directory ", dir_path_horse, " created.\n")
} else {
  cat("Directory ", dir_path_horse, " already exists.\n")
}

# Save the counts matrix and barcodes for horse
writeMM(counts_matrix_horse, file = file.path(dir_path_horse, "matrix.mtx"))
write.table(data.frame(barcodes_horse), file = file.path(dir_path_horse, "barcodes.tsv"), sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Extract the horse gene symbols
horse_genes <- rownames(counts_matrix_horse)

# Read your custom horse-human ortholog mapping table
# Replace 'human_horse_symbols_mapping.tsv' with your actual filename
orthologs_horse <- read.table("human_horse_symbols_mapping.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create a mapping data frame for horse
mapping_data_horse <- orthologs_horse %>% select(Horse_Symbol, Human_Symbol)

# Define the conversion function for horse genes
convert_horse_to_human <- function(gene_list, mapping_data) {
  # Create a mapping dictionary from Horse_Symbol to Human_Symbol
  mapping_dict <- setNames(mapping_data$Human_Symbol, mapping_data$Horse_Symbol)
  
  # Map horse genes to human genes
  converted_genes <- mapping_dict[gene_list]
  
  # Create a data frame to keep the original order
  result_df <- data.frame(
    Original_Gene = gene_list,
    Human_Gene = converted_genes,
    stringsAsFactors = FALSE
  )
  
  # Replace NA (genes without a human ortholog) with the original gene name
  result_df$Final_Gene <- ifelse(is.na(result_df$Human_Gene), result_df$Original_Gene, result_df$Human_Gene)
  result_df
}

# Convert horse genes to human genes
conversion_result_horse <- convert_horse_to_human(horse_genes, mapping_data_horse)
final_gene_list_horse <- conversion_result_horse$Final_Gene

# Ensure the final gene list has unique names
final_gene_list_horse <- make.unique(final_gene_list_horse)

# Write the gene list to 'genes.tsv' for horse
write.table(data.frame(final_gene_list_horse), file = file.path(dir_path_horse, "genes.tsv"), sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Read the counts matrix and associated files for horse
counts_matrix_horse <- readMM(file = file.path(dir_path_horse, "matrix.mtx"))
genes_horse <- read.table(file = file.path(dir_path_horse, "genes.tsv"), header = FALSE, sep = "\t")
gene_names_horse <- genes_horse$V1
barcodes_horse <- read.table(file = file.path(dir_path_horse, "barcodes.tsv"), header = FALSE, sep = "\t")
barcode_names_horse <- barcodes_horse$V1

# Assign row and column names to the counts matrix
rownames(counts_matrix_horse) <- gene_names_horse
colnames(counts_matrix_horse) <- barcode_names_horse

# Create a Seurat object for horse with updated gene names
horse_seurat <- CreateSeuratObject(
    counts = counts_matrix_horse,
    meta.data = horse@meta.data)

# ---------------------------
# Load Human Data
# ---------------------------

# Load human datasets
human <- readRDS("human.rds")
human_seurat <- CreateSeuratObject(
  counts =  human@assays$RNA$counts,
  meta.data = human@meta.data
)

human_seurat <- AddMetaData(human_seurat, metadata = "Human_healthy_ST",
                            col.name = "orig.ident")

library(biomaRt)

# Use biomaRt to get the mapping
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get a conversion table from Ensembl IDs to HGNC symbols
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart
)

# Ensure no duplicated Ensembl IDs in the mapping
gene_mapping <- gene_mapping[!duplicated(gene_mapping$ensembl_gene_id), ]

write.table(gene_mapping, "gene_mapping_human_ENSG-HGNC.tsv", sep = "\t",
            quote = F, row.names = F)

# Extract current gene names (Ensembl IDs)
current_gene_names <- rownames(human_seurat)

# Create a named vector for mapping (Ensembl ID → HGNC Symbol)
gene_map <- setNames(gene_mapping$hgnc_symbol, gene_mapping$ensembl_gene_id)

# Ensure empty HGNC symbols are replaced with NA
gene_map[gene_map == ""] <- NA

# Create a set to track already assigned gene symbols
assigned_symbols <- c()

# Initialize new gene names vector
new_gene_names <- current_gene_names  # Default to Ensembl IDs

# Assign HGNC symbols while ensuring uniqueness
for (i in seq_along(current_gene_names)) {
  ensembl_id <- current_gene_names[i]
  hgnc_symbol <- gene_map[ensembl_id]
  
  # If HGNC symbol exists and hasn't been used before, assign it
  if (!is.na(hgnc_symbol) && !(hgnc_symbol %in% assigned_symbols)) {
    new_gene_names[i] <- hgnc_symbol
    assigned_symbols <- c(assigned_symbols, hgnc_symbol)  # Mark as used
  }
}

# Assign the new unique names to the Seurat object
rownames(human_seurat) <- new_gene_names

# ---------------------------
# Merge All Datasets
# ---------------------------

# Before merging make sure no duplicated gene names are included
sum(duplicated(rownames(horse_seurat)))
sum(duplicated(rownames(human_seurat)))
sum(duplicated(rownames(mouse_seurat)))

# Merge mouse, horse, and human datasets
obj <- merge(mouse_seurat, y = list(horse_seurat, human_seurat), 
             add.cell.ids = c("mouse", "horse", "human"), 
             project = "obj_merged")

# Subset to include only cells with at least 1000 features
obj <- subset(obj, nFeature_RNA > 1000)

# Join layers in the merged object before splitting
obj <- JoinLayers(obj)

# Split the RNA assay by 'orig.ident'
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

# ---------------------------
# Data Preprocessing
# ---------------------------

# Normalize, identify variable features, and scale the data
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# ---------------------------
# Initial Clustering Without Integration
# ---------------------------

# Perform initial clustering
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Save UMAP plot of unintegrated data
png(paste0(prefix, "_UMAPplot_unintegrated.png"), width = 1600, height = 800)
DimPlot(obj, reduction = "umap.unintegrated", group.by = "orig.ident")
dev.off()

# ---------------------------
# Data Integration
# ---------------------------

# ---------------------------
# Visualization Function
# ---------------------------

# Function to create UMAP plots for each integration method
plotIntegratedData <- function(obj, keyword_string, prefix){
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

  p1 <- DimPlot(obj, reduction = red_name, group.by = clust_name, label = TRUE, label.size = 2.5)
  p2 <- DimPlot(obj, reduction = red_name, group.by = "orig.ident", label = TRUE, label.size = 2.5)

  combined_plot <- p1 + p2 + plot_layout(ncol = 2)

  umapplot_filename <- paste0(prefix, "_", keyword_string, "_UMAPplot.png")
  png(umapplot_filename, width = 3200, height = 1600)
  print(combined_plot)
  dev.off()
  cat("UMAP plot was saved to file: ", umapplot_filename, "\n")
}

#
## Perform CCA integration
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  dims = 1:10,
  verbose = FALSE
)

plotIntegratedData(obj, "cca", prefix)

## Perform RPCA integration
options(future.globals.maxSize = 12 * 1024^3)
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  dims = 1:10,
  k.anchor = 30,
  verbose = FALSE
)

plotIntegratedData(obj, "rpca", prefix)

## Perform Harmony integration
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

plotIntegratedData(obj, "harmony", prefix)

## Perform FastMNN integration
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

plotIntegratedData(obj, "mnn", prefix)


library(reticulate)
use_condaenv("/work/vetmed_data/mamba/envs/scvi", required = TRUE)

# Perform scVI integration
# Ensure scVI environment is set up as per scVI-tools instructions
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/work/vetmed_data/mamba/envs/scvi/", verbose = FALSE
)

plotIntegratedData(obj, "scvi", prefix)


# Save integrated RDS object
saveRDS(obj, paste0(prefix, "object_postIntegration.rds"))

