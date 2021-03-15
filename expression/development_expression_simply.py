
import os, pysam, math, tqdm
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
genome.add_file('genes', 'all_matches_1000000.dat')
# add tissue data.
genome.add_file('tissues','/media/ngs/data/development_project_PRJNA270632/summary_PRJNA270632.csv', type='full')
# read gene list.
genes=pd.read_csv(genome.genes, delim_whitespace=True)

# Let's look at distribution of distance from sat3 block to genes.
#plotting.histplot(genes['distance_to_repeat'].values,'dist to repeats', bin_count=400)

# Cut-off
max_dist=1000000

genes_sel = genes.loc[ genes['distance_to_repeat']<max_dist]
# Reduce by unique.
genes_unique = genes_sel['names'].unique()
genes_unique = genes_sel.loc[genes_sel['names'].str.contains('LOC'),'names'].unique()
#plotting.histplot(genes_sel['distance_to_repeat'].values,'dist to repeats', bin_count=100)

# Load tissues.
tissue_data = pd.read_csv(genome.tissues, sep='\t')
tissue_names = tissue_data['tissue'].astype('string').unique()
tissue_max_counts={ x:0 for x in tissue_names}

for gene in tqdm.tqdm(genes_unique):
	DATA = tissue_data[ tissue_data['symbol']==gene]
	if DATA.empty:
		continue
	# Tissue with maximum expression.
	max_tissue=DATA.loc[DATA['mean'].idxmax(),'tissue']
	#print(max_tissue)
	tissue_max_counts[max_tissue]=tissue_max_counts[max_tissue]+1

plotting.bar(tissue_max_counts)