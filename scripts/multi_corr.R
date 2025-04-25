# Load libraries
library(ggplot2)
library(reshape2)
library(ggpubr)

# Ens 44938 209 169 166 84 294 748 107 162 1145 711 845 535

# Create the data frame
df <- read.table(header = TRUE, text = "
sample N_transcripts KD4 KD11 KM4 KM11 OvaD4 OvaD11 OvaM4 OvaM11 OviD4 OviD11 OviM4 OviM11
DK 55931 276 220 247 110 368 930 140 222 1353 910 964 671
Dova 57699 250 225 234 128 387 993 150 231 1441 920 1042 742
Dovi 59547 279 234 273 121 372 927 148 220 1503 961 1083 742
MK 54768 261 227 230 119 346 909 135 210 1366 858 994 662
Mova 56357 276 215 230 118 388 927 158 215 1378 875 1015 694
Movi 58187 285 224 276 128 388 953 136 232 1486 945 1064 754
CK 55296 253 228 239 112 334 885 141 208 1318 876 981 681
Cova 56811 248 211 240 112 369 933 169 214 1377 901 992 713
Covi 58392 259 227 255 117 365 904 147 220 1388 905 1053 687
K 52693 246 213 218 108 347 864 138 211 1285 849 922 641
Ova 54032 255 209 222 109 335 902 142 215 1320 846 957 651
Ovi 55348 274 221 224 108 369 923 148 226 1418 884 1033 703
")

# Reshape data for plotting
df_long <- melt(df, id.vars = c("sample", "N_transcripts"),
                variable.name = "Condition", value.name = "Expression")

# Create the scatter plot with Pearson correlation
ggscatter <- ggplot(df_long, aes(x = N_transcripts, y = Expression)) +
  geom_point(size = 2, color = "firebrick1") +
  facet_wrap(~Condition, scales = "free_y") +
  stat_cor(method = "pearson", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), size = 3) +
  theme_bw() +
  labs(title = "Scatter plots: N_transcripts vs Conditions",
       x = "Size of transcriptome",
       y = "Number DE transcripts") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))

# Print the plot
print(ggscatter)

setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/sufnaMohamed/assembly_ind_and_comb/kallisto_quants/blast_filtering/DESeq2_analysis/heatmaps_PCAs')
ggsave("scatter_panels.pdf", ggscatter, width = 8, height = 9, units = "in", dpi = 300)
