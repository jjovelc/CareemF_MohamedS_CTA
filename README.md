# Transcriptome assembly pipeline

## Rational

In this project, libraries constructed for Saleh etal., 2025 () were used to conduct transcriptome assembly using the pipeline described below.

The goal of this project is to pursue the discovery of novel transcripts not included in the reference Gallus gallus in the Ensembl database, which is the most common repository used for the analysis of RNAseq libraries. The rational behind our hypothesis is that no transcriptome include all transcripts encoded by the genome. Expression of many transcripts is not constitutive but rather facultative. Since these libraries come from an experiment of virus infection in chicken it is possible that some transcripts expressed in response to virus infection are not included in the reference transcriptome. Other genetic and environmental factors may also be the cause of expression of transcripts not seen previously.

Once Trinotate has finished sucessfully:

1. Run script parse_NODE-PROTname_from_trinotateAnn.sh. The code may need a bit adjusment depending on the name of the directories to be processed. It will generate, inside each target directory, a file with suffix "_node_and_protID.txt". It will have two columns, the first one with the transcripts (contigs) IDs, and the second one with the proteins identified by BLAST. Example below. 

NODE_10000_length_3532_cov_22.214673_g5508_i0.p6	R212B_HUMAN
NODE_10000_length_3532_cov_22.214673_g5508_i0.p1	CARL3_MOUSE
NODE_10001_length_3532_cov_21.596346_g5336_i1.p1	CP131_DANRE
NODE_10002_length_3532_cov_18.651727_g2741_i1.p1	CRBN_CHICK
NODE_10003_length_3532_cov_17.756494_g5509_i0.p1	DTX2_HUMAN
NODE_10004_length_3532_cov_10.634313_g1368_i1.p1	ZNFX1_HUMAN
NODE_10006_length_3531_cov_18.715877_g388_i3.p1	MED12_MOUSE
NODE_10007_length_3531_cov_13.457167_g1246_i3.p1	SUCO_HUMAN
NODE_1000_length_7867_cov_16.939781_g604_i0.p1	CYFP2_MOUSE
NODE_10010_length_3530_cov_25.010854_g1052_i4.p1	ASND1_CHICK 

2. Generate a list of protein to search annotations in Uniprot. Example:

```bash
	cut -f 2 ovary_node_and_protID.txt > proteins
```

3. Retrieve annotations. Example:

```bash
	python ../retrieve_protein_annotations_uniprot.py proteins > ovary_uniprot_ann.txt
```


