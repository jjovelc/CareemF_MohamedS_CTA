# Transcriptome assembly pipeline

## Rational

In this project, libraries constructed for Saleh etal., 2025 () were used to conduct transcriptome assembly using the pipeline described below.

The goal of this project is to pursue the discovery of novel transcripts not included in the reference Gallus gallus in the Ensembl database, which is the most common repository used for the analysis of RNAseq libraries. The rational behind our hypothesis is that no transcriptome include all transcripts encoded by the genome. Expression of many transcripts is not constitutive but rather facultative. Since these libraries come from an experiment of virus infection in chicken it is possible that some transcripts expressed in response to virus infection are not included in the reference transcriptome. Other genetic and environmental factors may also be the cause of expression of transcripts not seen previously.

## Description of pipeline

### Preprocessing and assembly

1. Libraries were quality trimmed with fastq-mcf:

```bash
	TRIMMING_THRESHOLD=35

	# If files are not compressed change next line for: for FILE in *_R1.fq
	for FILE in *_R1.fq.gz
	do
    	fastq-mcf ~/useful_files/adapters.fa -o ${FILE/_R1.fq.gz/}_R1_trim${TRIMMING_THRESHOLD}.fq -o ${FILE/_R1.fq.gz/}_R2_trim${TRIMMING_THRESHOLD}.fq $FILE ${FILE/_R1/_R2} -k 0 -l 75  -w 3 -q $TRIMMING_THRESHOLD

	done
```

2. Crop reads required for achieving 100X coverage of transcriptome

```bash
	NUMBER_READS_NEEDED=116000000
	NUMBER_LINES=$(($NUMBER_READS_NEEDED * 4))
	NUMBER_IN_MILLIONS=${NUMBER_READS_NEEDED/000000/}

	for FILE in *_trim${TRIMMING_THRESHOLD}.fq
	do
        	head -${NUMBER_LINES} $FILE > ${FILE/.fq/}_${NUMBER_IN_MILLIONS}M.fq
	done
```

3. Concatenate all R1 and R2 files subsets (for each tissue)

```bash
	cat *R1*M.fq > all_R1.fq
	cat *R2*M.fq > all_R2.fq
```

4. Conduct assembly (for each tissue). For example:

```bash
	OUT="kidney_rnaSPAdes_assembly"
	rnaspades.py -t 48  -1 kidney_all_R1.fq -2 kidney_all_R2.fq -o $OUT
```

5. Predict transcripts with transdecoder. Run the following two command lines from inside the assebly results folder.

```bash
	TransDecoder.LongOrfs -t transcripts.fasta
	TransDecoder.Predict -t transcripts.fasta
```

The above commands will generate two key files, one for DNA and the other one for protein. For example, for the kidney assembly:

ovary_transcripts.fa.transdecoder.cds
ovary_transcripts.fa.transdecoder.pep

6. Extract only complete ORF, exlcuing 5' and 3' incomplete transcripts:

```bash
	perl extract_complete_ORFs.pl ovary_transcripts.fa.transdecoder.cds > ovary_transcripts.fa.transdecoder_complete.cds 
	perl extract_complete_ORFs.pl ovary_transcripts.fa.transdecoder.pep > ovary_transcripts.fa.transdecoder_complete.pep
``` 

7. Run the trinotate pipeline with script run_trinotate.sh in the directory were the two previous files with complete ORFs are.

### Postprocessing (after trinotate)

Once Trinotate has finished sucessfully:

1. Run script parse_NODE-PROTname_from_trinotateAnn.sh. The code may need a bit adjusment depending on the name of the directories to be processed. It will generate, inside each target directory, a file with suffix "_node_and_protID.txt". It will have two columns, the first one with the transcripts (contigs) IDs, and the second one with the proteins identified by BLAST. Example below. 

NODE_10000_length_3532_cov_22.214673_g5508_i0.p6	R212B_HUMAN
NODE_10000_length_3532_cov_22.214673_g5508_i0.p1	CARL3_MOUSE
NODE_10001_length_3532_cov_21.596346_g5336_i1.p1	CP131_DANRE
...

3. Retrieve annotations. Example:

```bash
	python ../retrieve_protein_annotations_uniprot.py proteins > ovary_uniprot_ann.txt
```


