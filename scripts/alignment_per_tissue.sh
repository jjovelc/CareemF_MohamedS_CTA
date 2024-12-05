#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -u  # Treat unset variables as errors

# Input control file
CONTROL_FILE=$1

# Derive file paths for DMV and Mass from the control file
DMV_FILE="${CONTROL_FILE/Con/DMV}"
MASS_FILE="${CONTROL_FILE/Con/Mass}"

# Extract tissue name from the control file
SUFFIX="${CONTROL_FILE#*_}"
TISSUE="${SUFFIX%%_*}"

# Validate input files
if [[ ! -f "$CONTROL_FILE" || ! -f "$DMV_FILE" || ! -f "$MASS_FILE" ]]; then
    echo "Error: One or more input files are missing!"
    echo "Control file: $CONTROL_FILE"
    echo "DMV file: $DMV_FILE"
    echo "Mass file: $MASS_FILE"
    exit 1
fi

# Check if bwa is installed
if ! command -v bwa &> /dev/null; then
    echo "Error: bwa is not installed or not in PATH!"
    exit 1
fi

echo "Processing tissue: $TISSUE"

# Index DMV to be used as reference
DMV_INDEX="${TISSUE}_DMV_INDEX"
bwa index -p $DMV_INDEX $DMV_FILE

# Run bwa alignments
# For Control
CONTROL_SAM="${CONTROL_FILE/fa/sam}"
bwa mem $DMV_INDEX $CONTROL_FILE > $CONTROL_SAM

# For Mass
MASS_SAM="${MASS_FILE/fa/sam}"
bwa mem $DMV_INDEX $MASS_FILE > $MASS_SAM

# Extract non-aligned reads
# For Control
CONTROL_NOT_DMV="${CONTROL_FILE/_fa/_not-in-DMV.fa}"
awk 'BEGIN {OFS = "\n"} /^@/ {next} $3 == "*" {print ">" $1, $10}' "$CONTROL_SAM" > "$CONTROL_NOT_DMV"

# For Mass
MASS_NOT_DMV="${MASS_FILE/_fa/_not-in-DMV.fa}"
awk 'BEGIN {OFS = "\n"} /^@/ {next} $3 == "*" {print ">" $1, $10}' "$MASS_SAM" > "$MASS_NOT_DMV"

# Index Mass non-aligned reads
MASS_NOT_DMV_INDEX="${TISSUE}_Mass_INDEX"
bwa index -p $MASS_NOT_DMV_INDEX $MASS_NOT_DMV

# Align Control non-aligned reads to Mass non-aligned reads
ONLY_CONTROL="${CONTROL_FILE/_fa/_not-in-DMV-Mass.fa}"
bwa mem $MASS_NOT_DMV_INDEX $CONTROL_NOT_DMV > $ONLY_CONTROL

# Combine all unique transcripts
FINAL_OUTPUT="${TISSUE}_novel_transcriptome.fa"
cat "$DMV_FILE" "$MASS_NOT_DMV" "$ONLY_CONTROL" > "$FINAL_OUTPUT"

echo "Final novel transcriptome for tissue $TISSUE saved to: $FINAL_OUTPUT"
