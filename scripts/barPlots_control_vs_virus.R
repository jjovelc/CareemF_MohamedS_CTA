# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create mean data (from user input)
data_mean <- data.frame(
  Group = c("Control", "DMV", "Mass"),
  KD4 = c(253.33, 268.33, 274),
  KD11 = c(222, 226.33, 222),
  KM4 = c(244.67, 251.33, 245.33),
  KM11 = c(113.67, 119.67, 121.67),
  OvaD4 = c(356, 375.67, 374),
  OvaD11 = c(907.33, 950, 929.67),
  OvaM4 = c(152.33, 146, 143),
  OvaM11 = c(214, 224.33, 219),
  OviD4 = c(1361, 1432.33, 1410),
  OviD11 = c(894, 930.33, 892.67),
  OviM4 = c(1008.67, 1029.67, 1024.33),
  OviM11 = c(693.67, 718.33, 703.33)
)

# Create standard deviation data (from user input)
data_std <- data.frame(
  Group = c("Control", "DMV", "Mass"),
  KD4 = c(5.51, 15.95, 12.12),
  KD11 = c(9.54, 7.09, 6.24),
  KM4 = c(8.96, 19.86, 26.56),
  KM11 = c(2.89, 9.07, 5.51),
  OvaD4 = c(19.16, 10.02, 24.25),
  OvaD11 = c(24.17, 37.27, 22.12),
  OvaM4 = c(14.74, 5.29, 13),
  OvaM11 = c(6, 5.86, 11.53),
  OviD4 = c(37.64, 75.37, 66.09),
  OviD11 = c(15.72, 27.02, 46.11),
  OviM4 = c(38.79, 60.45, 35.92),
  OviM11 = c(17.01, 40.99, 46.70)
)

# Convert data to long format for ggplot
data_long <- pivot_longer(data_mean, cols = -Group, names_to = "Library", values_to = "Mean")
std_long <- pivot_longer(data_std, cols = -Group, names_to = "Library", values_to = "SD")

# Merge mean and standard deviation data
data_plot <- left_join(data_long, std_long, by = c("Group", "Library"))

# Define colors (matching violin plot)
group_colors <- c("Control" = rgb(0, 0, 1, alpha = 0.7),  # Blue with transparency
                  "DMV" = rgb(1, 0, 0, alpha = 0.7),      # Red with transparency
                  "Mass" = rgb(0, 1, 0, alpha = 0.7))     # Green with transparency

# Create bar plot with error bars
p <- ggplot(data_plot, aes(x = Library, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.7), width = 0.3) +
  scale_fill_manual(values = group_colors) +
  labs(x = "Libraries",
       y = "Number of Differentially Expressed Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12 * 1.5),  # Rotate x-axis labels & enlarge
    axis.text.y = element_text(size = 12 * 1.5),  # Enlarge y-axis labels
    axis.title = element_text(size = 14 * 1.5),  # Enlarge axis titles
    legend.text = element_text(size = 12),  # Enlarge legend text
    legend.title = element_text(size = 14)  # Enlarge legend title
  )

# Save the plot
ggsave("barplot_DEG_with_SD.png", p, width = 10, height = 6)

# Display the plot
print(p)

