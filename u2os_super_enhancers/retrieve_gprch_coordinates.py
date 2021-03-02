import os, pysam, math
import pandas as pd
from gen import gen
import pysam
import utilities

data_file="~/Dropbox/TEAD_paper/data/genomics/U2OS_aacr_super_enhancers/super_enhancer_list.csv"

data = pd.read_csv(data_file,sep=',',header=None)
data.columns=['chr', 'start', 'stop']

genome=gen("/media/ngs/data/genomes/GRCh38.p13/", "GCF_000001405.39_GRCh38.p13_genomic.fna")
genome = gen("/media/ngs/data/genomes/hg19/", "hg19.fa")


fasta = genome.get_genome_pysam()

# Loop over super-enhancers. Extract out region. Write to new fasta file.
seqs=[]
for idx, row in data.iterrows():
    seqs.append(fasta.fetch(row['chr'], row['start'], row['stop']))

out_fa="/home/matt/Dropbox/code_bank/sat3/u2os_super_enhancers/enhancers_seq.fa"
utilities.write_fasta_from_list(out_fa, seqs)