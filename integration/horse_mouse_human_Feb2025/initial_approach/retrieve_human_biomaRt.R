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