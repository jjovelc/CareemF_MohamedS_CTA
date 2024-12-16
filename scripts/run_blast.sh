#!/usr/bin/bash

FASTA_FILE="$1"
DATABASE="ggallus_blast_db/gallus_transcriptome_db"
OUTFILE="${FASTA_FILE/_contigs.fa/_blast_matches.txt}"

blastn -query $FASTA_FILE -db $DATABASE -out $OUTFILE -evalue 1e-6 -outfmt 6
