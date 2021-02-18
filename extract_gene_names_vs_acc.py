# Convert the transcript name to a GENE name. Need to read a gbff file to search for gene name.




import os, pysam
import pandas as pd
import Bio.GenBank.Scanner
#from BCBio.GFF import GFFExaminer
from Bio import SeqIO


# File
in_file = '/media/ngs/data/genomes/GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.gbff'

L=list()
c=0
for seq_record in SeqIO.parse(in_file,'genbank'):

    # get accession number
    acc = seq_record.id
    # Get gene name.
    for feat in seq_record.features:
        if feat.type == 'gene':
            name = feat.qualifiers['gene']

    # Add to list.
    L.append( [acc,name[0]] )
    c=c+1

with open('acc_vs_name_list.txt', 'w') as f:
    for item in L:
        f.write("%s %s\n" %(item[0], item[1]))
