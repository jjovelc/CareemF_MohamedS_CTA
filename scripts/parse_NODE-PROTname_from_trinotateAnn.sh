#!/usr/bin/bash

for DIR in *trinotate
do
	cut -f 1,2 ${DIR}/ggallus.fasta.transdecoder.complete_trinotate_annotations.txt | sed 's/\^.*$//' > ${DIR}/${DIR/trinotate/node_and_protID.txt}
done
