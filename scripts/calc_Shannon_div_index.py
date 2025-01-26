import numpy as np
from collections import Counter

def calculate_shannon_diversity(fastq_file, chunk_size=100000):
   from Bio import SeqIO
   from collections import defaultdict
   import math
   
   kmer_counts = defaultdict(int)
   total_kmers = 0
   k = 5  # k-mer size
   
   # Process in chunks
   with open(fastq_file) as f:
       for record in SeqIO.parse(f, "fastq"):
           seq = str(record.seq)
           # Generate k-mers
           for i in range(len(seq) - k + 1):
               kmer = seq[i:i+k]
               kmer_counts[kmer] += 1
               total_kmers += 1
               
           # Calculate Shannon index periodically
           if total_kmers >= chunk_size:
               shannon = 0
               for count in kmer_counts.values():
                   p = count / total_kmers
                   shannon -= p * math.log2(p)
               print(f"Current Shannon Index: {shannon}")
               
   return shannon

