import pandas as pd
import sys

infile = sys.argv[1]
blast_report = sys.argv[2]
outfile = infile.replace("_contigs.fa", "_novel_transcripts.ids")

# Load BLAST output
blast_results = pd.read_csv(blast_report, sep="\t", header=None)
blast_results.columns = [
    "query_id", "subject_id", "identity", "alignment_length", "mismatches", 
    "gap_opens", "query_start", "query_end", "subject_start", "subject_end", 
    "evalue", "bit_score"
]

# Parse contig lengths from query_id
# Assuming the contig length is embedded in the query_id like "NODE_12611_length_3514_cov_..."
blast_results["contig_length"] = blast_results["query_id"].str.extract(r"_length_(\d+)").astype(int)

# Calculate alignment coverage as a percentage of the contig length
blast_results["alignment_coverage"] = (blast_results["alignment_length"] / blast_results["contig_length"]) * 100

# Filter results to retain only alignments covering at least 50% of the contig length
filtered_results = blast_results[blast_results["alignment_coverage"] >= 25]

# Get unique query IDs (matching contigs) from the filtered results
matching_contigs = set(filtered_results["query_id"])

# Load all contigs from the FASTA file
all_contigs = set()
with open(infile, "r") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            contig_id = line.split()[0][1:]  # Extract the contig ID (remove ">")
            all_contigs.add(contig_id)

# Identify non-matching contigs
non_matching_contigs = all_contigs - matching_contigs

# Save non-matching contigs to a file
with open(outfile, "w") as f:
    for contig in non_matching_contigs:
        f.write(f"{contig}\n")

