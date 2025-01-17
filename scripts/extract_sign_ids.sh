#!/usr/bin/bash

for FILE in *_q0.05_apeglm_FC2_annotated.tsv; do awk '$3 > 1' $FILE | cut -f 1 | sed 1d | sed '1s/^/transcript\n/' > ${FILE/_q0.05_apeglm_FC2_annotated.tsv/_up.txt}; done

for FILE in *_q0.05_apeglm_FC2_annotated.tsv; do awk '$3 < 1' $FILE | cut -f 1 | sed '1s/^/transcript\n/' > ${FILE/_q0.05_apeglm_FC2_annotated.tsv/_down.txt}; done


