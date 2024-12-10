#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G

start_time=$(date +%s)
echo "Job started on: $(date)"

####### Load environment and run kallisto ########

INDEX=$1
END1=$2
END2=$3
OUTDIR="/work/vetmed_data/jj/projects/juanJovel/pipelines/trinotate/alignments/quantification-2"
INDEX_DIR="/work/vetmed_data/jj/projects/juanJovel/pipelines/trinotate/alignments/quantification-2"

eval "$(conda shell.bash hook)"
conda activate kallisto

# Run kallisto quant with specified options
~/software/kallisto/kallisto quant -i "${INDEX_DIR}/${INDEX}" --bias -b 100 -o "${OUTDIR}/$(basename "${END1/_R1.fq/}_${INDEX/_hybrid_genome.idx/}_kallisto")" \
    --threads 4 $END1 $END2

end_time=$(date +%s)
echo "Job ended on: $(date)"

# Calculate and display job duration
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Job duration: $hours hours, $minutes minutes, $seconds seconds"

