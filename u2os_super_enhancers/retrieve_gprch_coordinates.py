import os, pysam, math
import pandas as pd
from gen import gen
import pysam

data_file="~/Dropbox/TEAD_paper/data/genomics/U2OS_aacr_super_enhancers/super_enhancer_list.csv"

data = pd.read_csv(data_file,sep=',',header=None)
data.columns=['chr','start','stop']

genome=gen("/media/ngs/data/genomes/GRCh38.p13/", "GCF_000001405.39_GRCh38.p13_genomic.fna")


fasta = genome.get_genome_pysam()

