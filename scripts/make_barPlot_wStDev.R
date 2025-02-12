# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Read the data
data <- read.table("downregulated_DE_transc.tsv", header=TRUE, sep="\t")

# Reshape data from wide to long format for grouped bars
plot_data <- data %>%
  filter(Metric %in% c("Average", "Ensembl")) %>%
  gather(key="Sample", value="Value", -Metric) %>%
  mutate(Metric = factor(Metric, levels = c("Average", "Ensembl")))

# Get StDev data separately
stdev_data <- data %>%
  filter(Metric == "StDev") %>%
  gather(key="Sample", value="StDev", -Metric) %>%
  select(Sample, StDev)

# Ensure error bars align with Average bars only
average_data <- plot_data %>%
  filter(Metric == "Average") %>%
  left_join(stdev_data, by="Sample")

# Define dodge width
dodge_width <- 0.8
bar_width <- 0.85  # Increase bar width

# Create the plot
p <- ggplot(plot_data, aes(x=Sample, y=Value, fill=Metric)) +
  # Add grouped bars with black stroke
  geom_bar(stat="identity", 
           position=position_dodge(width=dodge_width),
           width=bar_width,
           color="black") +
  # Add error bars only for Average and ensure alignment
  geom_errorbar(data=average_data,
                aes(ymin=Value-StDev, ymax=Value+StDev, 
                    x=as.numeric(as.factor(Sample)) - dodge_width/4),
                width=0.2, linewidth=0.9) +
  # Customize colors
  #scale_fill_manual(values=c("Average"="firebrick1", "Ensembl"="gray80")) +
  scale_fill_manual(values=c("Average"="green3", "Ensembl"="gray80")) +
  # Customize theme
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  # Labels
  labs(x = "Samples",
       y = "Count",
       title = "Average vs Ensembl Counts with Standard Deviation",
       subtitle = "Error bars shown for Average values only") +
  # Adjust y-axis limits to accommodate error bars
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Save the plot
ggsave("down_paired_barplot.pdf", p, width = 12, height = 8)
ggsave("down_paired_barplot.png", p, width = 12, height = 8, dpi = 300)
