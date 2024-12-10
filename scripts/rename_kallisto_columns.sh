#!/bin/bash

# Extract the second column from CHGI_and_sample_ids.txt and transpose it to create a header file
cut -f 2 CHGI_and_sample_ids.txt | bash /Users/juanjovel/jj/scripts/bash/transpose.bash > header

# Loop through directories matching the pattern kallistoRes_*
for DIR in kallistoRes_*
do
  # Remove the first line (current header) from all_samples_counts.tsv and prepend the new header
  sed '1d' "${DIR}/all_samples_counts.tsv" | cat header - > "${DIR}/tmp_counts.tsv"
  mv "${DIR}/tmp_counts.tsv" "${DIR}/all_samples_counts.tsv"

  # Remove the first line (current header) from all_samples_tpms.tsv and prepend the new header
  sed '1d' "${DIR}/all_samples_tpms.tsv" | cat header - > "${DIR}/tmp_tpms.tsv"
  mv "${DIR}/tmp_tpms.tsv" "${DIR}/all_samples_tpms.tsv"

  # Calculate the number of columns in the header
  NUM_COLS=$(($(awk 'END{print NF}' header) - 1))

  # Extract the directory name before the first '-' character
  DIR_NAME=$(echo "$DIR" | sed 's/-.*$//')

  # Trim the files to the number of columns in the header
  cut -f 1-"$NUM_COLS" "${DIR}/all_samples_counts.tsv" > "${DIR_NAME}_counts.tsv"
  cut -f 1-"$NUM_COLS" "${DIR}/all_samples_tpms.tsv" > "${DIR_NAME}_tpms.tsv"
done
