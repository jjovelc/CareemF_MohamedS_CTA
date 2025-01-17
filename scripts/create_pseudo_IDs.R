setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sufnaMohamed/assembly_ind_and_comb/kallisto_results/blast')

# Load universe and gene-to-GO mapping files
universe <- read.table("universe_sorted.txt", header = TRUE, stringsAsFactors = FALSE)
universe <- as.character(universe$genes)
tail(universe)
gene_to_GO <- read.table("transc_2_GO_All.tsv", header = TRUE, sep = '\t', fill = TRUE, quote = "", stringsAsFactors = FALSE)

# Identify novel transcripts explicitly by naming pattern
novel_transcripts <- unique(gene_to_GO$LocusID[grepl("^NODE_", gene_to_GO$LocusID)])

# Start numbering from the last Ensembl ID
last_ensembl_id <- 10072690
novel_ids <- paste0("ENSGALT000", seq(last_ensembl_id + 1, last_ensembl_id + length(novel_transcripts)))

# Create a mapping table
id_mapping <- data.frame(
  OriginalID = novel_transcripts,
  NewID = novel_ids
)

# Save the ID mapping table
write.table(id_mapping, "id_mapping.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Replace only novel IDs in gene_to_GO
gene_to_GO$LocusID <- ifelse(
  grepl("^NODE_", gene_to_GO$LocusID),
  id_mapping$NewID[match(gene_to_GO$LocusID, id_mapping$OriginalID)],
  gene_to_GO$LocusID  # Retain original IDs for existing Ensembl transcripts
)

# Replace only novel IDs in the universe
universe <- ifelse(
  grepl("^NODE_", universe),
  id_mapping$NewID[match(universe, id_mapping$OriginalID)],
  universe  # Retain original IDs for existing Ensembl transcripts
)


# Save the updated files
write.table(gene_to_GO, "transc_2_GO_All_with_pseudoIDs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(universe, "universe_with_pseudoIDs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
