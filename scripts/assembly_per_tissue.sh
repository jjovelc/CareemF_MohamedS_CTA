#!/usr/bin/bash


### STEP #1: Assembly with rnaSPAdes

#***********************************************#
#********** CHECK OR CHANGE PARAMETER **********#
#***********************************************#

TISSUES=(Kidney Ovary Oviduct)
DATA_DIRECTORY="/work/careem_lab/sufna/"

for TISSUE in ${TISSUES[@]}
do
	R1="${DATA_DIRECTORY}/${TISSUE}_R1.fq"
	R2="${DATA_DIRECTORY}/${TISSUE}_R2.fq"
	OUTDIR="${DATA_DIRECTORY}/${TISSUE}_rnaSPAdes_assembly"

	sbatch submit_3jobs.slurm "$R1" "$R2" "$OUTDIR"
done


