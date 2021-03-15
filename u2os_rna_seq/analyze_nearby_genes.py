# Look at tissue-specific expression of genes nearby Sat3 repeats.

import os, pysam, math
import pandas as pd
from gen import gen
import plotting
import numpy as np
# Define file system.
genomes_dir ='/media/ngs/data/genomes/chm13_v1.0/'
#genomes_dir ='/Users/franklin/Dropbox/TEAD_paper/data/genomics/chm13_v1.0/'
genome_fa = 'chm13.draft_v1.0.fasta'
genome = gen(genomes_dir, genome_fa)

# Add repeat bin file.
genome.add_file('binned_repeats','binned_CATTCC.bed')

# Add rna mapped to genome.
genome.add_file('rna','RNA_aligned.bam')

# add gene list.
genome.add_file('genes', 'all_matches_1M.dat')
# add tissue data.
genome.add_file('tissues','/media/ngs/data/ncbi_expression/summary_PRJEB4337.csv', type='full')
# read gene list.
genes=pd.read_csv(genome.genes, delim_whitespace=True)

# Let's look at distribution of distance from sat3 block to genes.
#plotting.histplot(genes['distance_to_repeat'].values,'dist to repeats', bin_count=400)

# Load sperm rna-seq.
rna_counts = pd.read_csv('/media/ngs/data/U2OS_riboZero_PRJNA251691/summary_table.csv', sep=',')

# Cut-off
max_dist=10000

genes_sel = genes.loc[ genes['distance_to_repeat']<max_dist]
# Reduce by unique.
genes_unique = genes_sel['names'].unique()
idx=[]
for gene in genes_unique:
    # Get index of gene.
    idx.append( rna_counts.index[ rna_counts['name']==gene].values[0])


# Perform averaging and std'ing.
rna_counts['mean_tpm']= rna_counts.iloc[:,1:5].mean(axis=1)
rna_counts['std_tpm'] = rna_counts.iloc[:,1:5].mean(axis=1)
rna_counts_red=rna_counts.loc[idx]

plotting.histplot(rna_counts_red['mean_tpm'],'tpm',bin_count=150)