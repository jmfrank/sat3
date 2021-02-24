import os, pysam, math
import pandas as pd
from gen import gen

from utilities import extract_matching_strings, load_repeat_masker_data, get_fasta_faidx, generate_rna_dataframe
from utilities import nearest_distance, filter_pandas, filter_similar_transcripts
from utilities import load_bin_bam

# Define file system.
genomes_dir ='/media/ngs/data/genomes/chm13_v1.0/'
#genomes_dir ='/Users/franklin/Dropbox/TEAD_paper/data/genomics/chm13_v1.0/'
genome_fa = 'chm13.draft_v1.0.fasta'
genome = gen(genomes_dir, genome_fa)

# Add repeat bin file.
genome.add_field('binned_repeats','binned_CATTCC.bed')

# Add rna mapped to genome.
genome.add_field('rna','RNA_aligned.bam')

# add gene list.
genome.add_field('genes', 'all_matches_1M.dat')

# read gene list.
genes=pd.read_csv(genome.base_dir+genome.genes, delim_whitespace=True)
genes['names'] = ''

# retrieve the gene name from accession number.
names = pd.read_csv('acc_vs_name.txt', delim_whitespace=True, header=None)
names.columns = ['acc','names']

for index, rep_row in genes.iterrows():
    this_acc = rep_row['transcript_name']
    # find this acc in the names df.
    genes.loc[index, 'names'] = names.loc[names['acc']== this_acc, 'names'].values[0]


# Need to filter transcripts so we only have unique reads?