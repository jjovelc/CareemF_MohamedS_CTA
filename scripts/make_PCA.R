# Load required libraries
library(ggplot2)
library(tidyverse)
library(ggrepel)

# Read the data and metadata
data <- read.table("deregulated_abundance_4PCA.tsv", header = TRUE, row.names = 1)
metadata <- read.table("deregulated_metadata_4PCA.tsv", header = TRUE)

# Transpose the data matrix since now samples are columns
data_t <- t(data)

# Scale the data before PCA (this is crucial)
data_scaled <- scale(data_t)

# Perform PCA
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# Calculate variance explained correctly
var_explained <- round((pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100, 3)

# Print variance explained by first few PCs
print("Variance explained by first 4 PCs:")
print(head(var_explained, 4))

# Create a data frame for plotting
plot_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2]
)

# Add metadata information
plot_data$group <- metadata$group
plot_data$label <- metadata$label

# Set manual colors
group_colors <- c("kidney" = "#cc1e06", 
                  "ovary" = "#ecbf1c", 
                  "oviduct" = "#3160f6")

# Create the PCA plot with variance percentages in axis labels
p <- ggplot(plot_data, aes(x = PC1, y = PC2, fill = group, label = label)) +
  geom_point(size = 5, shape = 21, color = "black", stroke = 1) +  # Added yellow stroke
  geom_text_repel(size = 4, show.legend = FALSE) +
  scale_fill_manual(values = group_colors) +
  xlab(paste0("PC1 (", format(var_explained[1], nsmall = 3), "%)")) +
  ylab(paste0("PC2 (", format(var_explained[2], nsmall = 3), "%)")) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 24),  # Doubled axis title size
    axis.text = element_text(size = 20),   # Doubled axis text size
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  ggtitle("PCA of Differentially Expressed Genes Across Samples")

# Save the plot
ggsave("PCA_deregulated_abundance.pdf", p, width = 10, height = 8)

# Print summary of PCA to see variance explained by each component
print(summary(pca_result))
