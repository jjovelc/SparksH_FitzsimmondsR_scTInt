library(biomaRt)


setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/hollySparks/dataIntegration_wSeurat/horse_and_mouse')

horse_obj <- readRDS('horse.rds')
human_obj <- readRDS('human.rds')

# Set up Ensembl databases
horse_mart <- useMart("ensembl", dataset = "ecaballus_gene_ensembl")
horse_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

horse_genes <- rownames(horse_obj)
horse_genes <- head(horse_genes, n=100)

# Get horse to human orthologs
orthologs <- getLDS(attributes = c("external_gene_name"), 
                    filters = "external_gene_name", 
                    values = horse_genes, 
                    mart = horse_mart, 
                    attributesL = c("external_gene_name", "ensembl_gene_id"), 
                    martL = human_mart)

# Rename columns for clarity
colnames(orthologs) <- c("horse_gene", "human_gene", "human_ensembl_id")

# Create a named vector for easy mapping
gene_map <- setNames(orthologs$human_gene, orthologs$horse_gene)

# Get the current gene names
current_genes <- rownames(horse_obj)

# Map to human genes, keeping original names if no match is found
new_genes <- gene_map[current_genes]
new_genes[is.na(new_genes)] <- current_genes[is.na(new_genes)]

# Update the gene names in the Seurat object
rownames(horse_obj) <- new_genes