# Look at tissue-specific expression of genes nearby Sat3 repeats.

import os, pysam, math
from tqdm import tqdm
import pandas as pd
from gen import gen
import plotting
import numpy as np
pd.set_option('display.max_rows', 1000)

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
def calc_mean(df):
    # Average over all data columns (first column is gene name).
    df['mean_tpm'] = df.iloc[:, 1:df.shape[1]].mean(axis=1)
    df['std_tpm'] = df.iloc[:, 1:df.shape[1]].mean(axis=1)
    return df


# Load 10A rna-seq.
H1all = pd.read_csv('/media/ngs/data/H1_diff_ucsd_reference/summary_table_diff.csv',sep=',')

wt = calc_mean(H1all[["name","SRR488684","SRR488685"]])
tropho= calc_mean(H1all[["name","SRR486237","SRR486238"]])
mes = calc_mean(H1all[["name","SRR486239","SRR486240"]])
npc = calc_mean(H1all[["name","SRR486241","SRR486242"]])

# Perform averaging.

# Cut-off
max_dist = 50000
genes_sel = genes.loc[genes['distance_to_repeat'] < max_dist]

# Reduce by unique.
genes_unique = genes_sel['names'].unique()
idx = []
for gene in tqdm(genes_unique):
    # Get index of gene.
    idx.append(wt.index[wt['name'] == gene].values[0])

def get_ratio_pd(genes, control, experiment):
    ratio=pd.DataFrame({'name':genes,'control':control['mean_tpm'],'experiment':experiment['mean_tpm'],'ratio':experiment['mean_tpm'] / control['mean_tpm']})
    return ratio

ratio = get_ratio_pd(genes_unique, wt.loc[idx], mes.loc[idx])

# Get rid of entries where neither control or experiment are expressed.
SEL = (ratio['control']==0) & (ratio['experiment']==0)
ratio = ratio.loc[ ~SEL ]

# Now set ratio for exclusively expressed genes.
ratio.loc[ratio['ratio']==np.inf] =  10 ** 4
ratio.loc[ratio['ratio']==0] = 10 ** -4

plotting.histplot(np.log10(ratio['ratio']), 'log10(wt / mes)', bin_count=200, xlims=[-4.1, 4.1])


