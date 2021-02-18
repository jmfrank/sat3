# Convert the transcript name to a GENE name. Need to read a gbff file to search for gene name.




import os, pysam
import pandas as pd
import Bio.GenBank.Scanner
#from BCBio.GFF import GFFExaminer
from Bio import SeqIO


# File
in_file = '/media/ngs/data/genomes/GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.gbff'

L = pd.read_csv('/media/ngs/py_projects/sat3/acc_vs_name.txt', delim_whitespace=True, header=None)

# open og bed file.
og_bed_file = '/media/ngs/data/genomes/HG002_CCS_canu_paternal/RNA_aligned_to_HG002.bed'

BData = pd.read_csv(og_bed_file,sep='\t',header=None)

# Name list
names=list()

# Loop over entries, find gene name.
for index, row in BData.iterrows():
    # This accession.
    this_acc = row[3]
    # Find in list.
    this_name = L[L[0] == this_acc][1].iloc[0]
    print(index)
    names.append(this_name)

# Now replace column with names.
BData[3] = names

# Save as BED.
BData.to_csv('/media/ngs/data/genomes/HG002_CCS_canu_paternal/RNA_aligned_to_HG002_names.bed',sep='\t', index=False)
