#!/usr/bin/bash

for FILE in *_contigs.fa
do
	sbatch submit.slurm $FILE
done
