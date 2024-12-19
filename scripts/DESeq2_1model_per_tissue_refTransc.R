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

# Function to pull records from BioMart
pullRecords <- function(attributes, mart, filter_values) {
  records <- getBM(attributes = attributes, filters = "ensembl_transcript_id", 
                   values = filter_values, mart = mart)
  return(records)
}

# Function to extract the first annotation match
getFirstMatch <- function(records, transcripts) {
  first_annotation_df <- data.frame()
  for (transcript in transcripts) {
    hits <- which(records$ensembl_transcript_id == transcript)
    if (length(hits) > 0) {
      first <- hits[1]
      first_annotation <- records[first,]
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    }
  }
  return(first_annotation_df)
}

# Rename columns for clarity
renameColumns <- function(df) {
  colnames(df)[6] <- "GO_group"
  colnames(df)[7] <- "GO_definition"
  colnames(df)[8] <- "Ontology"
  return(df)
}

# Function to create a volcano plot
makeVolcanoPlot <- function(df, vp_file, lfcthresh = 1, sigthresh = 0.05, 
                            width = 800, height = 600, 
                            point_size = 0.8, highlight_size = 1.1,
                            custom_colors = list(
                              nonsig = "black",
                              lowfc = "dodgerblue",
                              upleg = "red",
                              downreg = "forestgreen"
                            )) {
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  png(vp_file, width = width, height = height, res = 120)
  plot(df$log2FoldChange, -log10(df$padj), 
       type = "n", 
       main = "Volcano Plot",
       xlab = expression("log"[2]*"(Fold Change)"),
       ylab = expression("-log"[10]*"(Adjusted P-value)"),
       xlim = c(-max(abs(df$log2FoldChange)) * 1.1, max(abs(df$log2FoldChange)) * 1.1),
       ylim = c(0, max(-log10(df$padj), na.rm = TRUE) * 1.1))
  grid(lty = 2, col = "grey90")
  abline(h = -log10(sigthresh), v = c(-lfcthresh, lfcthresh), lty = 2, col = "grey50")
  
  # Plot points by categories
  nonsig <- df$padj >= sigthresh
  points(df$log2FoldChange[nonsig], -log10(df$padj[nonsig]), 
         pch = 21, bg = custom_colors$nonsig, col = 'black', cex = point_size)
  
  lowfc <- df$padj < sigthresh & abs(df$log2FoldChange) <= lfcthresh
  points(df$log2FoldChange[lowfc], -log10(df$padj[lowfc]),
         pch = 21, bg = custom_colors$lowfc, col = 'black', cex = point_size)
  
  up <- df$padj < sigthresh & df$log2FoldChange > lfcthresh
  points(df$log2FoldChange[up], -log10(df$padj[up]),
         pch = 21, bg = custom_colors$upleg, col = 'black', cex = highlight_size)
  
  down <- df$padj < sigthresh & df$log2FoldChange < -lfcthresh
  points(df$log2FoldChange[down], -log10(df$padj[down]),
         pch = 21, bg = custom_colors$downreg, col = 'black', cex = highlight_size)
  
  dev.off()
}

extract_norm_data <- function (dds, result, prefix){
  allResfile <- paste(prefix, "withNormalizedCounts.tsv", sep = '_') 
  resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Transcript"
  write.table(resdata, file=allResfile, sep="\t", quote = F, row.names = F)
  cat("Results with normalized data were save in ", allResfile)
}

#################### FUNCTIONS END ####################


#################### IMPORT DATA STARTS ######################

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sufnaMohamed/assembly_ind_and_comb/kallisto_results/blast/DESeq2_analysis/with_ref_transc')

# Get list of all _counts.tsv files in current directory
counts_file <- "ensembl_counts.tsv"

outdir      <- gsub("_counts.tsv", "", counts_file)
  
  # Read the data
  data <- read.table(counts_file, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = TRUE)
  
  metadata <- read.table("../metadata.tsv", sep = '\t', header = T, row.names = 1)
  
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
  
  #################### ITERATE ALONG INFECTION-TIME POINT DATASETS ####################
  
  
  # Iterate and annotate results
  for (res_name in names(results_list)) {
    prefix <- gsub("res_", "", res_name)
    result <- results_list[[res_name]]
    
    # Clean transcript IDs in result
    result$ensembl_transcript_id <- gsub("\\.\\d+$", "", rownames(result))
    
    sign_res <- subset(result, padj < 0.05 & abs(log2FoldChange) >= 1)
    
    # Perform the merge
    annotated_results <- merge(
      sign_res,
      biomart_annotations,
      by.x = "ensembl_transcript_id",
      by.y = "ensembl_transcript_id",
      all.x = TRUE
    )
    
    # Check for unmatched rows
    unmatched_rows <- annotated_results[is.na(annotated_results$ensembl_gene_id), ]
    if (nrow(unmatched_rows) > 0) {
      warning(paste(nrow(unmatched_rows), "rows in result do not match any entry in final_annotations."))
    }
    
    # Sort by log2FoldChange
    annotated_results <- annotated_results[order(-annotated_results$log2FoldChange), ]
    
    # Save the annotated results to the outdir
    file_name <- paste0(prefix, "_q0.05_apeglm_FC2_annotated.tsv")
    write.table(annotated_results, file_name, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Generate a volcano plot in the outdir
    vp_file <- paste0(prefix, "_volcanoPlot.png")
    makeVolcanoPlot(result, vp_file)
    
    # Print the number of significant transcripts
    significant_count <- nrow(subset(annotated_results, padj < 0.05 & abs(log2FoldChange) >= 1))
    print(paste("Number of significant transcripts deregulated:", significant_count))
  }
