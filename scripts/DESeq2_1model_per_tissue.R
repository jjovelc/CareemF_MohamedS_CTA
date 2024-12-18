# Differential expression analysis
# Project: Long shank neonate lim
#
# By Juan Jovel (use it at your own risk)
# (juan.jovel@ucalgary.ca)
# 
# Last revision: Dec 16, 2024

library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(pheatmap)
library(biomaRt)
library(EnhancedVolcano)
library(ggrepel)

rm(list = ls())

#################### FUNCTIONS START ####################
# Function to parse custom annotations safely
parse_custom_annotations <- function(annotations_file) {
  # Read the annotation file
  annotations <- read.table(
    annotations_file,
    sep = "\t",
    header = TRUE,
    comment.char = "",
    quote = ""
  )
  
  print(colnames(annotations))
  print(head(annotations))
  
  # Extract required fields
  annot_parsed <- data.frame(
    ensembl_transcript_id = annotations$gene_id,  # Use gene_id as transcript_id
    ensembl_gene_id = annotations$gene_id,        # Use gene_id
    entrezgene_id = NA,                           # No entrezgene_id available
    external_gene_name = annotations$gene_id,     # Repeat gene_id
    wikigene_description = ifelse(
      grepl("Full=", annotations$sprot_Top_BLASTX_hit), 
      sub(".*Full=([^;]+).*", "\\1", annotations$sprot_Top_BLASTX_hit), 
      NA
    ), # Extract 'Full=...' field
    GO_group = ifelse(
      grepl("^GO", annotations$gene_ontology_BLASTX),
      sub("^([^`]+)`.*", "\\1", annotations$gene_ontology_BLASTX),
      NA
    ), # First GO term
    GO_definition = ifelse(
      grepl("\\^", annotations$gene_ontology_BLASTX),
      sub("^.*\\^([^`]+)\\^.*", "\\1", annotations$gene_ontology_BLASTX),
      NA
    ), # GO description
    Ontology = ifelse(
      grepl("\\^", annotations$gene_ontology_BLASTX),
      sub("^.*\\^([^`]+)\\^.*", "\\1", annotations$gene_ontology_BLASTX),
      NA
    ) # GO ontology
  )
  
  # Replace missing descriptions with NA
  annot_parsed$wikigene_description[annot_parsed$wikigene_description == annot_parsed$ensembl_transcript_id] <- NA
  
  return(annot_parsed)
}



# Function to parse custom annotations
parse_custom_annotations <- function(annot_file) {
  annot <- read.table(annot_file, sep = "\t", header = TRUE, comment.char = "", quote = "")
  
  # Extract required fields
  annot_parsed <- data.frame(
    ensembl_transcript_id = annot$gene_id,  # Use gene_id as transcript_id
    ensembl_gene_id = annot$gene_id,        # Use gene_id
    entrezgene_id = NA,                     # No entrezgene_id available
    external_gene_name = annot$gene_id,     # Repeat gene_id
    wikigene_description = sub(".*Full=([^;]+).*", "\\1", annot$sprot_Top_BLASTX_hit), # Extract 'Full=...' field
    GO_group = sub("^([^`]+)`.*", "\\1", annot$gene_ontology_BLASTX),   # First GO term
    GO_definition = sub("^.*\\^([^`]+)\\^.*", "\\1", annot$gene_ontology_BLASTX), # GO description
    Ontology = sub("^.*\\^([^`]+)\\^.*", "\\1", annot$gene_ontology_BLASTX) # GO ontology
  )
  
  # Replace missing descriptions with NA
  annot_parsed$wikigene_description[annot_parsed$wikigene_description == annot$gene_id] <- NA
  
  return(annot_parsed)
}


combine_annotations <- function(biomart_annot, custom_annot, transcripts) {
  # Concatenate BioMart and Custom annotations
  concatenated <- rbind(
    biomart_annot,
    custom_annot
  )
  
  # Deduplicate concatenated annotations based on `ensembl_transcript_id`
  concatenated <- concatenated[!duplicated(concatenated$ensembl_transcript_id), ]
  
  # Create a dataframe with all transcripts_clean as base
  combined <- data.frame(ensembl_transcript_id = transcripts, stringsAsFactors = FALSE)
  
  # Merge concatenated annotations with the base
  combined <- merge(combined, concatenated, by = "ensembl_transcript_id", all.x = TRUE)
  
  # Fill missing values with `NA` where necessary
  combined[is.na(combined)] <- NA
  
  return(combined)
}



# Function to extract records
pullRecords <- function(attributes, mart, filter_values){
  records <- getBM(attributes = attributes, filters = "ensembl_transcript_id", 
                   values = filter_values, mart = mart)
  
  return(records)
}


# Extract first hit only
getFirstMatch <- function(records, transcripts) {
  first_annotation_df <- data.frame()
  for (transcript in transcripts) {
    hits <- which(records$ensembl_transcript_id == transcript)
    if (length(hits) > 0) {
      first <- hits[1]
      first_annotation <- records[first,]
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    } else {
      first_annotation <- c(transcript, "_", "-", "-", "-", "-", "-", "-")
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    }
  }
  return(first_annotation_df)
}

renameColumns <- function(df){
  colnames(df)[6] <- "GO_group"
  colnames(df)[7] <- "GO_definition"
  colnames(df)[8] <- "Ontology"
  return(df)
}


makeVolcanoPlot <- function(df, vp_file, lfcthresh = 1, sigthresh = 0.05, 
                            width = 800, height = 600, 
                            point_size = 0.8, highlight_size = 1.1, node_size = 1.3,
                            custom_colors = list(
                              nonsig = "black",
                              lowfc = "dodgerblue",
                              upleg = "red",
                              downreg = "forestgreen",
                              nodes = "yellow"
                            )) {
  # Input validation
  if (!all(c("log2FoldChange", "padj") %in% colnames(df))) {
    stop("DataFrame must contain 'log2FoldChange' and 'padj' columns")
  }
  
  # Remove NA values
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  # Identify NODE_XXXX rows
  df$is_node <- grepl("NODE_", rownames(df))
  
  # Calculate plot limits with padding
  max_lfc <- max(abs(df$log2FoldChange), na.rm = TRUE)
  max_pval <- max(-log10(df$padj), na.rm = TRUE)
  xlim <- c(-max_lfc * 1.1, max_lfc * 1.1)
  ylim <- c(0, max_pval * 1.1)
  
  # Counts for legend
  total_genes <- nrow(df)
  nonsig_count <- sum(df$padj >= sigthresh, na.rm = TRUE)
  lowfc_count <- sum(df$padj < sigthresh & abs(df$log2FoldChange) <= lfcthresh, na.rm = TRUE)
  up_count <- sum(df$padj < sigthresh & df$log2FoldChange > lfcthresh, na.rm = TRUE)
  down_count <- sum(df$padj < sigthresh & df$log2FoldChange < -lfcthresh, na.rm = TRUE)
  novel_count <- sum(df$padj < sigthresh & df$is_node & abs(df$log2FoldChange) > lfcthresh, na.rm = TRUE)
  
  # Open PNG device
  png(vp_file, width = width, height = height, res = 120)
  
  # Base plot
  plot(df$log2FoldChange, -log10(df$padj), 
       type = "n", 
       main = sprintf("Volcano Plot\n(Total genes: %d)", total_genes),
       xlab = expression("log"[2]*"(Fold Change)"),
       ylab = expression("-log"[10]*"(Adjusted P-value)"),
       xlim = xlim, ylim = ylim)
  grid(lty = 2, col = "grey90")
  abline(h = -log10(sigthresh), v = c(-lfcthresh, lfcthresh), lty = 2, col = "grey50")
  
  # Non-significant points
  nonsig <- df$padj >= sigthresh
  points(df$log2FoldChange[nonsig], -log10(df$padj[nonsig]),
         pch = 21, bg = custom_colors$nonsig, col = 'black', cex = point_size)
  
  # Significant but low FC
  lowfc <- df$padj < sigthresh & abs(df$log2FoldChange) <= lfcthresh
  points(df$log2FoldChange[lowfc], -log10(df$padj[lowfc]),
         pch = 21, bg = custom_colors$lowfc, col = 'black', cex = point_size)
  
  # Upregulated
  up <- df$padj < sigthresh & df$log2FoldChange > lfcthresh
  points(df$log2FoldChange[up], -log10(df$padj[up]),
         pch = 21, bg = custom_colors$upleg, col = 'black', cex = highlight_size)
  
  # Downregulated
  down <- df$padj < sigthresh & df$log2FoldChange < -lfcthresh
  points(df$log2FoldChange[down], -log10(df$padj[down]),
         pch = 21, bg = custom_colors$downreg, col = 'black', cex = highlight_size)
  
  # Significant NODE points
  significant_nodes <- df$padj < sigthresh & df$is_node & abs(df$log2FoldChange) > lfcthresh
  points(df$log2FoldChange[significant_nodes], -log10(df$padj[significant_nodes]),
         pch = 21, bg = custom_colors$nodes, col = 'black', cex = node_size)
  
  # Legend
  legend_text <- c(
    sprintf("Non-significant (%d)", nonsig_count),
    sprintf("Padj < %.2f; |FC| <= %.1f (%d)", sigthresh, 2^lfcthresh, lowfc_count),
    sprintf("Padj < %.2f; FC > %.1f (%d)", sigthresh, 2^lfcthresh, up_count),
    sprintf("Padj < %.2f; FC < -%.1f (%d)", sigthresh, 2^lfcthresh, down_count),
    sprintf("Novel transcripts (%d)", novel_count)
  )
  
  legend("topright", 
         legend = legend_text, 
         pch = 21, pt.bg = unlist(custom_colors), 
         col = 'black', pt.cex = 1, bty = "n", cex = 0.8)
  
  # Close PNG device
  dev.off()
}




extract_norm_data <- function (dds, result, prefix){
  allResfile <- paste(prefix, "withNormalizedCounts.tsv", sep = '_') 
  resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Transcript"
  write.table(resdata, file=allResfile, sep="\t", quote = F, row.names = F)
  cat("Results with normalized data were save in ", allResfile)
}



#################### FUNCTIONS END ######################


#################### IMPORT DATA STARTS ######################


setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sufnaMohamed/assembly_ind_and_comb/kallisto_results/blast/DESeq2_analysis')
#args <- commandArgs(trailingOnly = T)
#infile <- argv[1]
counts_file <- "kidney_DMV_counts.tsv"
annotations_file <- gsub('counts', 'annotations', counts_file)

custom_annotations <- parse_custom_annotations(annotations_file)
head(custom_annotations)


# Construct the full path
full_path <- file.path(getwd(), counts_file)

# Read the data
data <- read.table(full_path, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = TRUE)

metadata <- read.table("metadata.tsv", sep = '\t', header = T, row.names = 1)

#################### IMPORT DATA ENDS.  ######################

# Create a combined factor for Line and Tissue
metadata$group <- factor(paste0(metadata$tissue, "_", metadata$infection, "_", metadata$time_point))

data_matrix <- as.matrix(round(data), digits = 0) 

# Import data matrix into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data_matrix, colData = metadata, design =~ group)

# For Kidney tissue
# 4 dpi
metadata_kidney_4dpi <- metadata[metadata$tissue == "Kidney" & metadata$time_point == "4dpi" , ]
dds_kidney_4dpi <- DESeqDataSetFromMatrix(countData = data_matrix[, metadata$tissue == "Kidney" & metadata$time_point == "4dpi"], 
                                  colData = metadata_kidney_4dpi, design =~ group)
dds_kidney_4dpi <- DESeq(dds_kidney_4dpi, fitType='local')

# 11 dpi
metadata_kidney_11dpi <- metadata[metadata$tissue == "Kidney" & metadata$time_point == "11dpi", ]
dds_kidney_11dpi <- DESeqDataSetFromMatrix(countData = data_matrix[, metadata$tissue == "Kidney" & metadata$time_point == "11dpi"], 
                                          colData = metadata_kidney_11dpi, design =~ group)
dds_kidney_11dpi <- DESeq(dds_kidney_11dpi, fitType='local')

# For Ovary tissue
# 4 dpi
metadata_ovary_4dpi <- metadata[metadata$tissue == "Ovary" & metadata$time_point == "4dpi", ]
dds_ovary_4dpi <- DESeqDataSetFromMatrix(countData = data_matrix[, metadata$tissue == "Ovary" & metadata$time_point == "4dpi"], 
                                          colData = metadata_ovary_4dpi, design =~ group)
dds_ovary_4dpi <- DESeq(dds_ovary_4dpi, fitType='local')

# 11 dpi
metadata_ovary_11dpi <- metadata[metadata$tissue == "Ovary" & metadata$time_point == "11dpi", ]
dds_ovary_11dpi <- DESeqDataSetFromMatrix(countData = data_matrix[, metadata$tissue == "Ovary" & metadata$time_point == "11dpi"], 
                                         colData = metadata_ovary_11dpi, design =~ group)
dds_ovary_11dpi <- DESeq(dds_ovary_11dpi, fitType='local')

# For Oviduct tissue
# 4 dpi
metadata_oviduct_4dpi <- metadata[metadata$tissue == "Oviduct" & metadata$time_point == "4dpi", ]
dds_oviduct_4dpi <- DESeqDataSetFromMatrix(countData = data_matrix[, metadata$tissue == "Oviduct" & metadata$time_point == "4dpi"], 
                                         colData = metadata_oviduct_4dpi, design =~ group)
dds_oviduct_4dpi <- DESeq(dds_oviduct_4dpi, fitType='local')

# 11 dpi
metadata_oviduct_11dpi <- metadata[metadata$tissue == "Oviduct" & metadata$time_point == "11dpi", ]
dds_oviduct_11dpi <- DESeqDataSetFromMatrix(countData = data_matrix[, metadata$tissue == "Oviduct" & metadata$time_point == "11dpi"], 
                                           colData = metadata_oviduct_11dpi, design =~ group)
dds_oviduct_11dpi <- DESeq(dds_oviduct_11dpi, fitType='local')


# Now you can get results for each tissue comparison
resultsNames(dds_kidney_4dpi)
resultsNames(dds_kidney_11dpi)
resultsNames(dds_ovary_4dpi)
resultsNames(dds_ovary_11dpi)
resultsNames(dds_oviduct_4dpi)
resultsNames(dds_oviduct_11dpi)

res_kid_DMV_4dpi  <- lfcShrink(dds_kidney_4dpi, coef="group_Kidney_DMV_4dpi_vs_Kidney_Cont_4dpi", type="apeglm")
res_kid_Mass_4dpi <- lfcShrink(dds_kidney_4dpi, coef="group_Kidney_Mass_4dpi_vs_Kidney_Cont_4dpi", type="apeglm")
res_kid_DMV_11dpi  <- lfcShrink(dds_kidney_11dpi, coef="group_Kidney_DMV_11dpi_vs_Kidney_Cont_11dpi", type="apeglm")
res_kid_Mass_11dpi <- lfcShrink(dds_kidney_11dpi, coef="group_Kidney_Mass_11dpi_vs_Kidney_Cont_11dpi", type="apeglm")

res_ova_DMV_4dpi  <- lfcShrink(dds_ovary_4dpi, coef="group_Ovary_DMV_4dpi_vs_Ovary_Cont_4dpi", type="apeglm")
res_ova_Mass_4dpi <- lfcShrink(dds_ovary_4dpi, coef="group_Ovary_Mass_4dpi_vs_Ovary_Cont_4dpi", type="apeglm")
res_ova_DMV_11dpi  <- lfcShrink(dds_ovary_11dpi, coef="group_Ovary_DMV_11dpi_vs_Ovary_Cont_11dpi", type="apeglm")
res_ova_Mass_11dpi <- lfcShrink(dds_ovary_11dpi, coef="group_Ovary_Mass_11dpi_vs_Ovary_Cont_11dpi", type="apeglm")

res_ovi_DMV_4dpi  <- lfcShrink(dds_oviduct_4dpi, coef="group_Oviduct_DMV_4dpi_vs_Oviduct_Cont_4dpi", type="apeglm")
res_ovi_Mass_4dpi <- lfcShrink(dds_oviduct_4dpi, coef="group_Oviduct_Mass_4dpi_vs_Oviduct_Cont_4dpi", type="apeglm")
res_ovi_DMV_11dpi  <- lfcShrink(dds_oviduct_11dpi, coef="group_Oviduct_DMV_11dpi_vs_Oviduct_Cont_11dpi", type="apeglm")
res_ovi_Mass_11dpi <- lfcShrink(dds_oviduct_11dpi, coef="group_Oviduct_Mass_11dpi_vs_Oviduct_Cont_11dpi", type="apeglm")

extract_norm_data(dds_kidney_4dpi, res_kid_DMV_4dpi, 'kid_DMV_4dpi')
extract_norm_data(dds_kidney_11dpi, res_kid_DMV_11dpi, 'kid_DMV_11dpi')
extract_norm_data(dds_kidney_4dpi, res_kid_Mass_4dpi, 'kid_Mass_4dpi')
extract_norm_data(dds_kidney_11dpi, res_kid_Mass_11dpi, 'kid_Mass_11dpi')

extract_norm_data(dds_ovary_4dpi, res_ova_DMV_4dpi, 'ova_DMV_4dpi')
extract_norm_data(dds_ovary_11dpi, res_ova_DMV_11dpi, 'ova_DMV_11dpi')
extract_norm_data(dds_ovary_4dpi, res_ova_Mass_4dpi, 'ova_Mass_4dpi')
extract_norm_data(dds_ovary_11dpi, res_ova_Mass_11dpi, 'ova_Mass_11dpi')

extract_norm_data(dds_oviduct_4dpi, res_ovi_DMV_4dpi, 'ovi_DMV_4dpi')
extract_norm_data(dds_oviduct_11dpi, res_ovi_DMV_11dpi, 'ovi_DMV_11dpi')
extract_norm_data(dds_oviduct_4dpi, res_ovi_Mass_4dpi, 'ovi_Mass_4dpi')
extract_norm_data(dds_oviduct_11dpi, res_ovi_Mass_11dpi, 'ovi_Mass_11dpi')

# Put results in a list of dataframes
results_list <- list(
  res_kid_DMV_4dpi = as.data.frame(res_kid_DMV_4dpi),
  res_kid_DMV_11dpi = as.data.frame(res_kid_DMV_11dpi),
  res_kid_Mass_4dpi = as.data.frame(res_kid_Mass_4dpi),
  res_kid_Mass_11dpi = as.data.frame(res_kid_Mass_11dpi),
  res_ova_DMV_4dpi = as.data.frame(res_ova_DMV_4dpi),
  res_ova_DMV_11dpi = as.data.frame(res_ova_DMV_11dpi),
  res_ova_Mass_4dpi = as.data.frame(res_ova_Mass_4dpi),
  res_ova_Mass_11dpi = as.data.frame(res_ova_Mass_11dpi),
  res_ovi_DMV_4dpi = as.data.frame(res_ovi_DMV_4dpi),
  res_ovi_DMV_11dpi = as.data.frame(res_ovi_DMV_11dpi),
  res_ovi_Mass_4dpi = as.data.frame(res_ovi_Mass_4dpi),
  res_ovi_Mass_11dpi = as.data.frame(res_ovi_Mass_11dpi)
)

##### ANNOTATION OF DE FEATURES #####
library(biomaRt)

# Create remote connection
my_mart <- useEnsembl('ensembl', dataset = "ggallus_gene_ensembl")

# Make list of attributes to retrieve
attributes <- c("ensembl_transcript_id",
                "ensembl_gene_id",
                "entrezgene_id",
                "external_gene_name",
                "wikigene_description",
                "name_1006",
                "definition_1006",
                "namespace_1003")

transcripts       <- row.names(res_kid_DMV_4dpi)
transcripts_clean <- gsub("\\.\\d+$", "", transcripts) # remove transc version
transc_annotations <- pullRecords(attributes, my_mart, transcripts_clean)
biomart_annotations <- getFirstMatch(transc_annotations, transcripts_clean)
biomart_annotations <- renameColumns(biomart_annotations)
biomart_annotations <- biomart_annotations[!grepl("NODE", biomart_annotations$ensembl_transcript_id), ]

# Filter custom_annotations to include only rows with NODE IDs present in transcripts_clean
custom_filtered <- custom_annotations[
  custom_annotations$ensembl_transcript_id %in% transcripts_clean, ]

# Combine BioMart and custom annotations
final_annotations <- rbind(biomart_annotations, custom_filtered)

#################### ITERATE ALONG INFECTION-TIME POINT DATASETS ####################

# Iterate and annotate results
for (res_name in names(results_list)) {
  prefix <- gsub("res_", "", res_name)
  result <- results_list[[res_name]]
  
  # Clean transcript IDs in result
  result$ensembl_transcript_id_clean <- gsub("\\.\\d+$", "", rownames(result))
  
  sign_res <- subset(result,  padj < 0.05 & abs(log2FoldChange) >= 1)
  
  # Clean transcript IDs in final_annotations (ensure consistency)
  final_annotations$ensembl_transcript_id_clean <- gsub("\\.\\d+$", "", final_annotations$ensembl_transcript_id)
  
  # Perform the merge
  annotated_results <- merge(
    sign_res,
    final_annotations,
    by.x = "ensembl_transcript_id_clean",
    by.y = "ensembl_transcript_id_clean",
    all.x = TRUE
  )
  
  # Check for unmatched rows
  unmatched_rows <- annotated_results[is.na(annotated_results$ensembl_gene_id), ]
  
  if (nrow(unmatched_rows) > 0) {
    warning(paste(nrow(unmatched_rows), "rows in result do not match any entry in final_annotations."))
  }
  
  # Sort by log2FoldChange
  annotated_results <- annotated_results[order(-annotated_results$log2FoldChange), ]
  
  # Save the results
  file_name <- paste0(prefix, "_q0.05_apeglm_FC2_annotated.tsv")
  write.table(annotated_results, file_name, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Generate a volcano plot
  vp_file <- paste0(prefix, "_volcanoPlot.png")
  makeVolcanoPlot(result, vp_file)
  
  # Print number of significant transcripts
  significant_count <- nrow(subset(annotated_results, padj < 0.05 & abs(log2FoldChange) >= 1))
  print(paste("Number of significant transcripts deregulated:", significant_count))
}
