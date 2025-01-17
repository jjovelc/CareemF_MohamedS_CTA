# Load necessary libraries
args <- commandArgs(trailingOnly = TRUE)  # Capture command-line arguments

# Ensure a directory is provided
if (length(args) == 0) {
  stop("Usage: Rscript convert_IDs.R <directory>")
}

# Set the working directory to the specified directory
input_dir <- args[1]
setwd(input_dir)

# List all files matching the pattern
de_files <- list.files(pattern = "_q0.05_apeglm_FC2_annotated.tsv")


# Read the mapping table
id_mapping <- read.table("../../id_mapping.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Iterate over each file
for (file in de_files) {
  
  # Read the DE results file
  de_results <- read.table(file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  transc <- as.character(de_results$ensembl_transcript_id_clean)
  
  # Replace only NODE_ IDs in the DE results
  de_results$ensembl_transcript_id_clean <- ifelse(
    grepl("^NODE_", transc),
    id_mapping$NewID[match(transc, id_mapping$OriginalID)],
    transc  # Retain original IDs for existing Ensembl transcripts
  )
  
  # Define the output file name
  outfile <- gsub("_annotated.tsv", "_annotated_renamed.tsv", file)
  
  # Write the updated DE results to a new file
  write.table(de_results, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Print a status message
  cat("Processed file:", file, "->", outfile, "\n")
}
