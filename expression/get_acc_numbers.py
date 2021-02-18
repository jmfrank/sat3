# Convert the transcript name to a GENE name. Need to read a gbff file to search for gene name.




import os, pysam
import pandas as pd
import Bio.GenBank.Scanner
#from BCBio.GFF import GFFExaminer
from Bio import SeqIO

# Extract all acc numbers.

# File
in_file = '/media/ngs/data/genomes/GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.gbff'

L = pd.read_csv('/media/ngs/py_projects/sat3/acc_vs_name.txt', delim_whitespace=True, header=None)

# open og bed file.
og_bed_file = '/media/ngs/data/genomes/chm13_v1.0/RNA_aligned.bed'

BData = pd.read_csv(og_bed_file,sep='\t',header=None)

BData[3].to_csv('/media/ngs/py_projects/sat3/all_acc_names.csv',index=False)
