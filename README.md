# Transcriptome assembly pipeline

## Rational

In this project, libraries constructed for Saleh etal., 2025 () were used to conduct transcriptome assembly using the pipeline described below.

The goal of this project is to pursue the discovery of novel transcripts not included in the reference Gallus gallus in the Ensembl database, which is the most common repository used for the analysis of RNAseq libraries. The rational behind our hypothesis is that no transcriptome include all transcripts encoded by the genome. Expression of many transcripts is not constitutive but rather facultative. Since these libraries come from an experiment of virus infection in chicken it is possible that some transcripts expressed in response to virus infection are not included in the reference transcriptome. Other genetic and environmental factors may also be the cause of expression of transcripts not seen previously.

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


