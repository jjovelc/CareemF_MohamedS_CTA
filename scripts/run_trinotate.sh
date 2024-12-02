#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate trinotate

# Change this directory to reflect your pwd
### ******************* ###
DIR="/work/vetmed_data/jj/projects/juanJovel/pipelines/trinotate/$1"

cd "$DIR"
echo "Working directory: $DIR"
TRINOTATE_DATA_DIR="${DIR}/trinotate_db"
DB_PATH="${DIR}/ggallus_transc.db"
CPU=24

# Set TRINOTATE_DATA_DIR environment variable
export TRINOTATE_DATA_DIR="${TRINOTATE_DATA_DIR}"
export PYTHONPATH="/home/juan.jovel/mambaforge/envs/trinotate/bin/"

export PERL5LIB=$(mamba info --base)/envs/prokka_env/lib/site_perl/5.26.2

# Run EggnogMapper (if installed)
echo "Running emapper..."
/home/juan.jovel/mambaforge/envs/trinotate/bin/emapper.py -i "${DIR}/transcripts.fasta.transdecoder.complete.pep" \
	--dmnd_db trinotate_db/EGGNOG_DATA_DIR/eggnog_proteins.dmnd  --data_dir trinotate_db/EGGNOG_DATA_DIR/ \
-o "${DIR}/eggnog_mapper" --cpu "${CPU}"
echo "emapper run completed..."

# Run Infernal
echo "Running cmscan..."
cmscan --cpu "${CPU}" --tblout "${DIR}/infernal.out" \
    "${TRINOTATE_DATA_DIR}/Rfam.cm" "${DIR}/transcripts.fasta.transdecoder.complete.cds" > "${DIR}/cmscan.log"
echo "cmscan run completed..."

# Load results into the database
echo "Loading results into sqlite3 database..."
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --LOAD_swissprot_blastp "${DIR}/uniprot_sprot.ncbi.blastp.outfmt6"
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --LOAD_pfam "${DIR}/TrinotatePFAM.out"
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --LOAD_signalp6 "${DIR}/sigPout_summary.signalp5"
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --LOAD_EggnogMapper "${DIR}/eggnog_mapper.emapper.annotations"
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --LOAD_tmhmmv2 "${DIR}/tmhmm.v2.out"
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --LOAD_swissprot_blastx "${DIR}/uniprot_sprot.ncbi.blastx.outfmt6"
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --LOAD_infernal "${DIR}/infernal.out"
echo "Results loaded to the sqlite3 database..."

# Generate the annotation report
echo "Generating the trinotate report..."
~/software/Trinotate-Trinotate-v4.0.2/Trinotate --db "${DB_PATH}" --report > "${DIR}/ggallus_transc_Trinotate_report.tsv"
echo "Trinotate report was generated..."

# Extract annotations
echo "Parsing trinotate report..."
cut -f 1 "${DIR}/transcripts.fasta.transdecoder.complete.cds" | sed 's/>//' | \
grep -F -f - "${DIR}/ggallus_transc_Trinotate_report.tsv" | cut -f 1,3,13 | \
awk '$2 != "."' > "${DIR}/ggallus.fasta.transdecoder.complete_trinotate_annotations.txt"
echo "pipeline completed..."
echo "check file: ggallus.fasta.transdecoder.complete_trinotate_annotations.txt"
