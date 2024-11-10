#!/usr/bin/bash

for DIR in *files
do
	cut -f 1,2 ${DIR}/ggallus.fasta.transdecoder.complete_trinotate_annotations.txt | sed 's/\^.*$//' > ${DIR}/${DIR/files/node_and_protID.txt}
done
