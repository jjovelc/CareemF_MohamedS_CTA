# Load required libraries
library(ggplot2)
library(reshape2)

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sufnaMohamed/assembly_ind_and_comb/kallisto_quants/blast_filtering/DESeq2_analysis/heatmaps')

# Create the data frame
data <- read.table("deregulated_abundance.tsv", header = TRUE,
                   row.names = 1, sep = "\t")

# Melt the data for ggplot2
melted_data <- melt(as.matrix(data))
colnames(melted_data) <- c("Group", "Sample", "Value")

# Create the heatmap
p <- ggplot(melted_data, aes(x = Sample, y = Group, fill = Value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Value), size = 5) +
  scale_fill_gradient(low = "white", high = "purple") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # doubled size from default
    axis.text.y = element_text(size = 12),  # doubled size from default
    axis.title.x = element_text(size = 14),  # doubled size for axis titles
    axis.title.y = element_text(size = 14),  # doubled size for axis titles
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  labs(
    x = "Libraries",
    y = "Assemblies",
    fill = "Value"
  ) +
  coord_equal()

# Save the plot
ggsave("deregulated_heatmap.pdf", p, width = 12, height = 8)

