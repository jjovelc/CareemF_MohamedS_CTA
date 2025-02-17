setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sufnaMohamed/assembly_ind_and_comb/kallisto_quants/blast_filtering/DESeq2_analysis/heatmaps_PCAs')

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# Read the data
data <- read.table("deregulated_abundance.tsv", header = TRUE, row.names = 1)

# Convert to long format for ggplot
data_long <- melt(data, variable.name = "Library", value.name = "DEG_Count")

# Ensure rownames are stored in a column for proper grouping
data_long$Assembly <- rep(rownames(data), times = ncol(data))

# Classify assemblies into three groups: Control (C), DMV (D), Mass (M), or mark as "Other"
data_long$Group <- case_when(
  grepl("^C", data_long$Assembly) ~ "Control",
  grepl("^D", data_long$Assembly) ~ "DMV",
  grepl("^M", data_long$Assembly) ~ "Mass",
  TRUE ~ NA_character_  # Mark irrelevant samples as NA
)

# Remove irrelevant rows (where Group is NA)
data_filtered <- filter(data_long, !is.na(Group))

# Convert Group to a factor with ordered levels
data_filtered$Group <- factor(data_filtered$Group, levels = c("Control", "DMV", "Mass"))

# Verify correct classification
print(table(data_filtered$Group))  # Should show only three categories

# Perform Kruskal-Wallis test (non-parametric ANOVA)
kruskal_test <- kruskal.test(DEG_Count ~ Group, data = data_filtered)
print(kruskal_test)

# Violin plot
violin_plot <- ggplot(data_filtered, aes(x = Group, y = DEG_Count, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(y = "Number of Differentially Expressed Genes",
       x = "Assembly Group") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16 * 1.3),  # Increase axis text by 1.5
    axis.title = element_text(size = 18 * 1.3), # Increase axis title by 1.5
    legend.text = element_text(size = 12),      # Slightly larger legend text
    legend.title = element_text(size = 14)      # Slightly larger legend title
  ) +
  scale_fill_manual(values = c("Control" = "blue", "DMV" = "red", "Mass" = "green"))

# Boxplot
box_plot <- ggplot(data_filtered, aes(x = Group, y = DEG_Count, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  labs(y = "Number of Differentially Expressed Genes",
       x = "Assembly Group") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16 * 1.3),  # Increase axis text by 1.5
    axis.title = element_text(size = 18 * 1.3), # Increase axis title by 1.5
    legend.text = element_text(size = 12),      # Slightly larger legend text
    legend.title = element_text(size = 14)      # Slightly larger legend title
  ) +
  scale_fill_manual(values = c("Control" = "blue", "DMV" = "red", "Mass" = "green"))

# Save plots
ggsave("violin_plot_DEG_3groups.png", violin_plot, width = 6, height = 5)
ggsave("box_plot_DEG_3groups.png", box_plot, width = 6, height = 5)

# Display plots
print(violin_plot)
print(box_plot)

