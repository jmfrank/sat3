# Look at tissue-specific expression of genes nearby Sat3 repeats.

import os, pysam, math
from tqdm import tqdm
import pandas as pd
from gen import gen
import plotting
import numpy as np

# Define file system.
genomes_dir = '/media/ngs/data/genomes/chm13_v1.0/'
# genomes_dir ='/Users/franklin/Dropbox/TEAD_paper/data/genomics/chm13_v1.0/'
genome_fa = 'chm13.draft_v1.0.fasta'
genome = gen(genomes_dir, genome_fa)

# Add repeat bin file.
genome.add_file('binned_repeats', 'binned_CATTCC.bed')

# Add rna mapped to genome.
genome.add_file('rna', 'RNA_aligned.bam')

# add gene list.
genome.add_file('genes', 'all_matches_1000000.dat')
# add tissue data.
genome.add_file('tissues', '/media/ngs/data/ncbi_expression/summary_PRJEB4337.csv', type='full')
# read gene list.
genes = pd.read_csv(genome.genes, delim_whitespace=True)


# Let's look at distribution of distance from sat3 block to genes.
# plotting.histplot(genes['distance_to_repeat'].values,'dist to repeats', bin_count=400)

# Funfction for loading average data.
def load_rna_data(file_name):
    rna_counts = pd.read_csv(file_name, sep=',')
    # Average over all data columns (first column is gene name).
    rna_counts['mean_tpm'] = rna_counts.iloc[:, 1:rna_counts.shape[1]].mean(axis=1)
    rna_counts['std_tpm'] = rna_counts.iloc[:, 1:rna_counts.shape[1]].mean(axis=1)
    return rna_counts


# Load 10A rna-seq.
mcf10A = load_rna_data('/media/ngs/data/10A_rna_seq_rRNA_depl_GSE150003/summary_table.csv')
NT2 = load_rna_data('/media/ngs/data/U2OS_riboZero_PRJNA251691/summary_table.csv')

# Cut-off
max_dist = 50000
genes_sel = genes.loc[genes['distance_to_repeat'] < max_dist]

# Reduce by unique.
genes_unique = genes_sel['names'].unique()
idx = []
for gene in tqdm(genes_unique):
    # Get index of gene.
    idx.append(mcf10A.index[mcf10A['name'] == gene].values[0])

mcf10A_red = mcf10A.loc[idx]
NT2_red = NT2.loc[idx]

ratio = NT2_red['mean_tpm'] / mcf10A_red['mean_tpm']
NT2_red['ratio'] = ratio

ratio = ratio.replace(np.inf, 10 ** 4)
ratio.loc[(ratio == 0) | (ratio.isna())] = 10 ** -4
plotting.histplot(np.log10(ratio), 'log10(NT2 / MCF10A)', bin_count=200, xlims=[-4.1, 4.1])

# For genes are only expressed in u2os, what are their properties?
NT2_red.loc[NT2_red['ratio'] == np.inf, 'name']
# Problem here is that each gene can have multiple matches to sat3. take closest? Largest?
dist = np.empty([NT2_red.shape[0], 1])

for c, p in enumerate(NT2_red.iterrows()):
    this_name = p[1]['name']
    matches = genes_sel.loc[genes['names'] == this_name]
    # Closest match?
    dist[c]=matches['distance_to_repeat'].min()

NT2_red['distance_to_repeat_min']=dist
dist_exclusive_expression=u2os_red.loc[NT2_red['ratio'] == np.inf, 'distance_to_repeat_min'].values
