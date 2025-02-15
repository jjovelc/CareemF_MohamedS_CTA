# This script will retrieve BioMart annotations


# Load the biomaRt package
library(biomaRt)

# Connect to the Ensembl BioMart
ensembl <- useMart("ensembl")

# Select the Gallus gallus (chicken) dataset
gallus <- useDataset("ggallus_gene_ensembl", mart = ensembl)

# Specify the attributes (fields) to retrieve.
# In this case, we'll retrieve the transcript ID, associated GO term, 
# and GO term annotations (name and namespace).
attributes <- c("ensembl_transcript_id",
                "go_id",
                "name_1006",
                "namespace_1003")

# Since we want *all* transcripts and their GOs, we do not set any filters.
results <- getBM(attributes = attributes,
                 mart = gallus)

# Save results
write.table(results, "ggallus_transcripts_and_GOs.txt", quote = F, row.names = F, sep = '\t')
