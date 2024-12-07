#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -u  # Treat unset variables as errors

eval "$(conda shell.bash hook)"
conda activate cd-hit

# Input control file
control_file=$1

# Derive file paths for DMV and Mass from the control file
dmv_file="${control_file/Con/DMV}"
mass_file="${control_file/Con/Mass}"

# Extract tissue name from the control file
suffix="${control_file#*_}"
tissue="${suffix%%_*}"

# Validate input files
if [[ ! -f "$control_file" || ! -f "$dmv_file" || ! -f "$mass_file" ]]; then
    echo "Error: One or more input files are missing!"
    echo "Control file: $control_file"
    echo "DMV file: $dmv_file"
    echo "Mass file: $mass_file"
    exit 1
fi

echo "Processing tissue: $tissue"

# Define output paths
control_clustered="${control_file/_annotated_novel_transcripts_na-ren.fa/}_clustered.fa"
dmv_clustered="${dmv_file/_annotated_novel_transcripts_na-ren.fa/}_clustered.fa"
mass_clustered="${mass_file/_annotated_novel_transcripts_na-ren.fa/}_clustered.fa"
final_output="${tissue}_final_novel_transcripts_clustered.fa"

# Cluster control file
echo "Clustering control file: $control_file"
cd-hit-est -i "$control_file" -o "$control_clustered" -c 0.95 -n 10 -M 16000 -d 0 || { echo "Error clustering control file"; exit 1; }

# Cluster DMV file
echo "Clustering DMV file: $dmv_file"
cd-hit-est -i "$dmv_file" -o "$dmv_clustered" -c 0.95 -n 10 -M 16000 -d 0 || { echo "Error clustering DMV file"; exit 1; }

# Cluster Mass file
echo "Clustering Mass file: $mass_file"
cd-hit-est -i "$mass_file" -o "$mass_clustered" -c 0.95 -n 10 -M 16000 -d 0 || { echo "Error clustering Mass file"; exit 1; }

# Combine clustered files
echo "Combining clustered files..."
cat "$control_clustered" "$dmv_clustered" "$mass_clustered" > combined_clustered.fa

# Final clustering
echo "Clustering combined file..."
cd-hit-est -i combined_clustered.fa -o "$final_output" -c 0.95 -n 10 -M 16000 -d 0 || { echo "Error clustering combined file"; exit 1; }

# Clean up intermediate files
rm -f "$control_clustered" "$dmv_clustered" "$mass_clustered" combined_clustered.fa

echo "Clustering and filtering complete for tissue: $tissue"
echo "Final output saved in: $final_output"
