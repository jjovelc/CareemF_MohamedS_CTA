#!/usr/bin/bash

for FILE in *_contigs.fa
do
	IDS="${FILE/_contigs.fa/_novel_transcripts.ids}"
	OUTFILE="${IDS/ids/fasta}"
	seqkit grep -f "$IDS" -o "$OUTFILE" "$FILE"
done




