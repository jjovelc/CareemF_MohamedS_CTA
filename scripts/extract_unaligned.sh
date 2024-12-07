#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.sam> <output.fasta>"
    exit 1
fi

# Input SAM file and output FASTA file
INPUT_SAM=$1
OUTPUT_FASTA=$2

# Check if the input SAM file exists
if [ ! -f "$INPUT_SAM" ]; then
    echo "Error: Input file $INPUT_SAM does not exist."
    exit 1
fi

# Extract non-aligned reads
echo "Extracting non-aligned reads from $INPUT_SAM..."
awk 'BEGIN {OFS = "\n"} /^@/ {next} $3 == "*" {print ">" $1, $10}' "$INPUT_SAM" > "$OUTPUT_FASTA"

# Check if the output was successfully created
if [ -f "$OUTPUT_FASTA" ]; then
    echo "Non-aligned reads saved to $OUTPUT_FASTA."
else
    echo "Error: Failed to create output file."
    exit 1
fi
