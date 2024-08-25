library(Seurat)
library(tidyverse)
library(reshape2)
library(Azimuth)

setwd("/Users/juanjovel/jj/data_analysis/hollySparks/dataIntegration_wSeurat/horse_and_mouse")

obj <- readRDS("mouse.rds")

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Extract the relevant data
meta_data <- obj@meta.data

# Select the columns for features and counts
data_to_plot <- meta_data[, c("nFeature_RNA", "nCount_RNA")]

# Add an identifier column
data_to_plot$Cell <- rownames(data_to_plot)

# Reshape the data
melted_data <- melt(data_to_plot, id.vars = "Cell")

# Generate the violin plot using ggplot2
ggplot(melted_data, aes(x = variable, y = value, fill = variable)) +
  geom_violin() +
  geom_jitter(width = 0.5, size = 0.025, alpha = 0.6) +
  theme_minimal() +
  labs(title = "Violin Plot of Features and Counts", x = "Metric", y = "Value") +
  scale_fill_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# Capitalize gene names in the counts matrix
counts_matrix <- GetAssayData(obj, layer = "counts")
rownames(counts_matrix) <- toupper(rownames(counts_matrix))
obj <- SetAssayData(obj, layer = "counts", new.data = counts_matrix)

# Capitalize gene names in the data matrix
data_matrix <- GetAssayData(obj, layer = "data")
rownames(data_matrix) <- toupper(rownames(data_matrix))
obj <- SetAssayData(obj, layer = "data", new.data = data_matrix)

# Capitalize gene names in the scale.data matrix
scale_data_matrix <- GetAssayData(obj, layer = "scale.data")
rownames(scale_data_matrix) <- toupper(rownames(scale_data_matrix))
obj <- SetAssayData(obj, layer = "scale.data", new.data = scale_data_matrix)

# Verify the changes
head(rownames(GetAssayData(obj, layer = "counts")), n=200)

# Check if CD79A is among the gene names
"F12" %in% rownames(GetAssayData(obj, layer = "counts"))

# Generate UMAP plot for the CD79A gene (now capitalized)
FeaturePlot(obj, features = "BC034090", reduction = "umap") + 
  ggtitle("UMAP Plot for CD79A Expression") +
  theme_minimal()
