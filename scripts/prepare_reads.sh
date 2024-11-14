#!/usr/bin/bash

TISSUES=(Kidney Ovary Oviduct)
DATA_DIRECTORY="/work/careem_lab/sufna/"

for TISSUE in ${TISSUES[@]}
do
        for DIR in "${DATA_DIRECTORY}/*_${TISSUE}"
        do
                OUTPUT_FILE="${DATA_DIRECTORY}/${TISSUE}_R1.fq"
                for FILE in ${DIR}/*_R1_25M.fq
                do
                        head -5500000 $FILE >> $OUTPUT_FILE
                done
                gzip ${TISSUE}_R1.fq

                OUTPUT_FILE="${DATA_DIRECTORY}/${TISSUE}_R2.fq"
                for FILE in ${DIR}/*_R2_25M.fq
                do
                        head -5500000 $FILE >> $OUTPUT_FILE
                done
                gzip ${TISSUE}_R2.fq
        done
done

