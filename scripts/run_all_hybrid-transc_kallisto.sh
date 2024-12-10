#!/usr/bin/bash

set -euo pipefail

# Index all hybrid transcriptomes

for FILE in *_hyb-transc.fa
do
	~/software/kallisto/kallisto index -i ${FILE/fa/idx} $FILE
done

DATA_DIR="/work/careem_lab/sufna"

# Reference transcriptomes (use spaces, not commas, between items)
REFS=(kidney_comb-assembly_hyb-transc.idx \
kidney_ind-assembly_hyb-transc.idx \
ovary_comb-assembly_hyb-transc.idx \
ovary_ind-assembly_hyb-transc.idx \
oviduct_comb-assembly_hyb-transc.idx \
oviduct_ind-assembly_hyb-transc.idx)

# Loop through each reference transcriptome index
for INDEX in "${REFS[@]}"
do
        echo "Processing reference: $INDEX"

        # Run alignments for each library
        for FILE in "$DATA_DIR"/*_R1.fq
        do
                PAIR_FILE="${FILE/_R1.fq/_R2.fq}"
                echo "Processing file: $FILE with paired file: $PAIR_FILE"

                # Submit each job with paired-end files
                sbatch submit_kallisto.slurm "$INDEX" "$FILE" "$PAIR_FILE"
        done
done
