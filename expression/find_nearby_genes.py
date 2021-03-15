# Find genetic elements nearby CATTCC repeat clusters.

import os, pysam, math, sys, tqdm
import pandas as pd
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from gen import gen
import plotting

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

# maximum distance between end of repeat and beginning of transcript:
max_distance=10**6

# Filter binned data base on simple threshold. Counts/ kb.
motif_threshold = 15
repeat_data = load_bin_bam(genomes_dir+genome.binned_repeats)
repeat_data = filter_pandas(repeat_data, count=[2, math.inf])
plotting.histplot(repeat_data['count'].values,'TEAD motifs/2kb', bin_count=150)

# Looking at histograms of the data, there's basically three densities of TEAD motifs. Very low, and two high. Let's look
# at data in each of these density ranges.
low_density=[0,70]
med_density=[70,150]
hi_density =[150,math.inf]

# Make a new column for classification.
repeat_data['class']=''
# Label data by class.
low_data = (repeat_data['count'] > low_density[0]) & (repeat_data['count'] <= low_density[1])
repeat_data.loc[low_data,'class']='low'
med_data = (repeat_data['count'] > med_density[0]) & (repeat_data['count'] <= med_density[1])
repeat_data.loc[med_data,'class']='medium'
high_data = (repeat_data['count'] > hi_density[0]) & (repeat_data['count'] <= hi_density[1])
repeat_data.loc[high_data,'class']='high'

# Histogram by chromosome?
plotting.histplot(repeat_data.loc[repeat_data['class'] == 'high', 'contig'], 'contig', title='high', label_rotation=45)

# load genome.
fasta = get_fasta_faidx( genome.base_dir+genome.genome_fa )
# load aligned rna bam. **make sure to create bam file index: "samtools index **.bam"
rna_bam = pysam.AlignmentFile(genome.base_dir+genome.rna)

#Create a dataframe for total output.
column_names = ['contig','ref_start','ref_end','repeat_type','motif_counts','transcript_name', 'start_position', 'length', \
                'distance_to_repeat', 'gene_loc']
L=list()

# Loop over repeatmasker entries that match the sat3 repeats, get their center location in the assembly. Find the closest transcript(s).
for index, rep_row in tqdm(repeat_data.iterrows()):
    sat_start=rep_row['start']
    sat_end = rep_row['stop']
    # Search for RNAs within 'max_distance'.
    rna_search_range = [max(sat_start - max_distance, 0), min(sat_end + max_distance, fasta.get_reference_length(rep_row['contig']))]
    rna_matches = rna_bam.fetch(reference=rep_row['contig'], start=rna_search_range[0], end=rna_search_range[1])

    # iterate over reads, append to matches dataframe
    mini_list=list()
    positions=list()
    for read in rna_matches:
        distance2repeat, gene_loc = nearest_distance( [sat_start, sat_end], [read.reference_start, read.reference_end])
        mini_list.append([rep_row['contig'], rep_row['start'], rep_row['stop'], rep_row['class'], rep_row['count'], read.query_name, read.reference_start, read.query_length, distance2repeat, gene_loc])
        # keep track of positions for filtering out later.
        positions.append(read.positions)

    # Keep all matches.
    L.extend(mini_list)

# Convert list to df.
DATA=pd.DataFrame(L,columns=column_names)
# Save DATA to file.
DATA.to_csv(genome.base_dir+'/all_matches_'+str(max_distance)+'.dat', sep='\t', index=False)

# Test reading/fetching from sam file.
#read = sat_bam.fetch("NC_000004.12",1,50000)
#for x in read:
#    print(str(x))




