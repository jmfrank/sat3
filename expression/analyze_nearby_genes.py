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
plotting.histplot(genes['distance_to_repeat'].values,'dist to repeats', bin_count=400)

# Cut-off
max_dist=1000000

genes_sel = genes.loc[ genes['distance_to_repeat']<max_dist]
# Reduce by unique.
genes_unique = genes_sel['names'].unique()
genes_unique = genes_sel.loc[genes_sel['names'].str.contains('LOC'),'names'].unique()
plotting.histplot(genes_sel['distance_to_repeat'].values,'dist to repeats', bin_count=100)

# Load tissues.
tissue_data = pd.read_csv(genome.tissues, sep='\t')

# for gene in genes_unique:
#     DATA = tissue_data[ tissue_data['symbol']==gene]
#     if DATA.shape[0] == 0:
#         print('Did not find gene in tissue data')
#         continue
#     print(DATA)
#     plotting.tissue_bar(DATA)
#     input()
# genes_unique

# Load sperm rna-seq.
rna_counts = pd.read_csv('/media/ngs/data/GSE144085_sperm_stem_cells/summary_table.csv', sep=',')
idx=[]
for gene in genes_unique:
    # Get index of gene.
    idx.append( rna_counts.index[ rna_counts['name']==gene].values[0])


# Perform averaging and std'ing.
rna_counts['mean_kit']= rna_counts.iloc[:,1:4].mean(axis=1)
rna_counts['std_kit'] = rna_counts.iloc[:,1:4].std(axis=1)
rna_counts['mean_PLPPR3']= rna_counts.iloc[:,5:8].mean(axis=1)
rna_counts['std_PLPPR3'] = rna_counts.iloc[:,5:8].std(axis=1)

rna_counts.loc[idx]

ratio=rna_counts.loc[idx,'mean_kit'] / rna_counts.loc[idx,'mean_PLPPR3']
ratio.loc[(ratio==0) | (ratio.isna())]=0.01
ratio=ratio.replace(np.inf,  100)
#plotting.histplot(ratio, 'ratio', bin_count=50)
plotting.histplot(np.log10(ratio), 'ratio', bin_count=50)

# Compare to all? Maybe counts aren't normalized / sample?
ratio_all=rna_counts['mean_kit'] / rna_counts['mean_PLPPR3']
ratio_all.loc[(ratio_all==0) | (ratio_all.isna())]=0.01
ratio_all=ratio_all.replace(np.inf,  100)
#plotting.histplot(ratio, 'ratio', bin_count=50)
plotting.histplot(np.log10(ratio_all), 'ratio', bin_count=500)
