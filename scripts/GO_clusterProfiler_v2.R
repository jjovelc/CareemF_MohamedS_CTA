# Gene Ontology Analysis with clusterProfiler

library(clusterProfiler)
library(tidyverse)

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)  # Capture command-line arguments

# Ensure a directory is provided
if (length(args) == 0) {
  stop("Usage: Rscript GO_clusterProfiler_v2.R <directory>")
}

# Set the working directory to the specified directory
input_dir <- args[1]
regulation <- args[2]

setwd(input_dir)

# List all files matching the pattern
my_patt <- paste0("_", regulation, ".txt")
de_files <- list.files(pattern = my_patt)

# Iterate over each file
for (file in de_files) {
  cat("File being processed: ", file, "\n")
  outfile <- gsub(".txt", "_GOenrich.tsv", file)
  out_dotplot <- gsub(".txt", "_GOenrich.png", file)
  
  # Read DE transcripts and universe
  genes <- read.table(file, header = TRUE, stringsAsFactors = FALSE) %>%
    pull(transcript)
  universe <- read.table("../universe_with_pseudoIDs.tsv", header = TRUE, stringsAsFactors = FALSE) %>%
    pull(genes)
  gene_to_go <- read.table("../transc_2_GO_All_with_pseudoIDs.tsv", header = TRUE, sep = '\t', fill = TRUE, quote = "", stringsAsFactors = FALSE)

  # Limit GO terms per transcript to reduce excessive bias
  max_gos <- 255
  gene_to_go_filtered <- gene_to_go %>%
    group_by(LocusID) %>%
    slice_head(n = max_gos) %>%
    ungroup()

  # Add weights to GO terms
  weight_table <- table(gene_to_go_filtered$LocusID)
  gene_to_go_filtered <- gene_to_go_filtered %>%
    mutate(weight = 1 / weight_table[LocusID])

  # GO Enrichment Function
  GOEnrichment <- function(genes_of_interest, universe, gene_to_GO, ontology) {
    if (ontology == 'ALL') {
      TERM2GENE <- gene_to_GO[, c('GOterm', 'LocusID', 'weight')]
      TERM2NAME <- gene_to_GO[, c('GOterm', 'GOdesc')]
    } else {
      TERM2GENE <- gene_to_GO[gene_to_GO$Ontology == ontology, c('GOterm', 'LocusID', 'weight')]
      TERM2NAME <- gene_to_GO[gene_to_GO$Ontology == ontology, c('GOterm', 'GOdesc')]
    }

    results <- enricher(gene = genes_of_interest,
                        universe = universe,
                        TERM2GENE = TERM2GENE[, c('GOterm', 'LocusID')],
                        TERM2NAME = TERM2NAME,
                        pAdjustMethod = 'BH',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        minGSSize = 10,
                        maxGSSize = 300)

    # Adjust enrichment results based on weights
    if (!is.null(results)) {
      term_weights <- TERM2GENE %>%
        group_by(GOterm) %>%
        summarize(weight = mean(weight, na.rm = TRUE))

      results@result <- results@result %>%
        left_join(term_weights, by = c("ID" = "GOterm")) %>%
        mutate(weight_adjusted_pvalue = pvalue * weight,
               adjusted_qvalue = p.adjust(weight_adjusted_pvalue, method = 'BH')) %>%
        arrange(weight_adjusted_pvalue)
    }
    
    return(results)
  }

  
  # Perform GO enrichment
  results <- GOEnrichment(genes_of_interest = genes,
                          universe = universe,
                          gene_to_GO = gene_to_go_filtered,
                          ontology = 'ALL')
  
  if (!is.null(results) && nrow(results@result) > 0 && any(results@result$p.adjust < 0.05)) {
    # Save simplified results
    write.table(results@result, outfile, sep = '\t', quote = FALSE, row.names = FALSE)

    # Generate dot plot
    png(out_dotplot, width = 1200, height = 800)
    p <- dotplot(results, showCategory = min(10, nrow(results@result)), font.size = 12)
    print(p)
    dev.off()
  } else {
    cat("No significant GO terms to plot.\n")
  }
  
}

