#!/usr/bin/bash

####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=256G

####### Set environment variables ###############

start_time=$(date +%s)
echo "Job started on: $(date)"
####### Run your script #########################
eval "$(conda shell.bash hook)"
conda activate spades-env

END1=$1
END2=$2
OUTDIR=$3

rnaspades.py -t 24  -1 "$END1"  -2 "$END2" -o "$OUTDIR"


end_time=$(date +%s)
echo "Job ended on: $(date)"

# Calculate the difference in times to get the duration
duration=$((end_time - start_time))

hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Job duration: $hours hours, $minutes minutes, $seconds seconds"

