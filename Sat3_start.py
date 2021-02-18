

# Scratch script for getting started with genomics.

import os, pysam
import pandas as pd

from utilities import extract_matching_strings, load_repeat_masker_data, get_fasta_faidx, generate_rna_dataframe
from utilities import nearest_distance, filter_pandas, filter_similar_transcripts

# Define file system.
genomes_dir ='/media/ngs/data/genomes/'

refseq_assembly_genome_dir ='GRCh38.p13/'
refseq_assembly_genome = genomes_dir+refseq_assembly_genome_dir+'GCF_000001405.39_GRCh38.p13_genomic.fna'
refseq_assembly_aligned_transcripts = genomes_dir+refseq_assembly_genome_dir+'RefSeq_transcripts_alignments/'+'GCF_000001405.39_GRCh38.p13_modelrefseq_alns.bam'
bamfile = pysam.AlignmentFile(refseq_assembly_aligned_transcripts)

# Let's look at the recently assembled genome
HG002_CCS_dir = 'HG002_CCS_canu_paternal/'
HG002_CCS_genome = genomes_dir + HG002_CCS_dir + 'GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.fna'
HG002_CCS_aligned_bam = genomes_dir + HG002_CCS_dir + 'aligned_to_GRCh38.p13.bam'
sat_bam = pysam.AlignmentFile(HG002_CCS_aligned_bam)
HG002_CCS_rna_aligned_bam = genomes_dir + HG002_CCS_dir + 'RNA_aligned_to_HG002.bam'
rna_bam = pysam.AlignmentFile(HG002_CCS_rna_aligned_bam)

# Repeat masker file.
HG002_CCS_rm = genomes_dir + HG002_CCS_dir + 'GCA_004796285.1_HG002_CCS_canu_paternal_1.0_rm.out'
repeats=['ATTCC', 'GGAAT', 'CCTTA', 'TAAGG', 'HSATII']
repeats=['ATTCC', 'GGAAT', 'CCTTA', 'TAAGG','HSATII']
# maximum distance between end of repeat and beginning of transcript:
max_distance=10**6

# Find lines in repeat masker that match repeat base strings.
sat3_matches = genomes_dir + HG002_CCS_dir + 'sat3_repeats.out'

extract_matching_strings(HG002_CCS_rm, sat3_matches, repeats)


# Filter sat3 data based on parameters.
repeat_data = load_repeat_masker_data(sat3_matches)
repeat_data = filter_pandas(repeat_data, score=[1000, float('inf')],length=[15, float('inf')])

# open up the genome for quick random access. Need to check if index was generated.
fasta = get_fasta_faidx( HG002_CCS_genome )

# Create a dataframe for total output.
column_names = ['contig','ref_start','ref_end','repeat_type','transcript_name', 'start_position', 'length', 'distance_to_repeat', 'gene_loc'] #, 'seq']
L=list()
# Loop over repeatmasker entries that match the sat3 repeats, get their center location in the assembly. Find the closest transcript(s).
for index, rep_row in repeat_data.iterrows():

    sat_start=rep_row['query_start']
    sat_end = rep_row['query_end']
    # Search for RNAs within 'max_distance'.
    rna_search_range = [max(sat_start - max_distance, 0), min(sat_end + max_distance, fasta.get_reference_length(rep_row['query_name']))]
    rna_matches = rna_bam.fetch(reference=rep_row['query_name'], start=rna_search_range[0], end=rna_search_range[1])


    # iterate over reads, append to matches dataframe
    mini_list=list()
    positions=list()
    for read in rna_matches:

        distance2repeat, gene_loc = nearest_distance( [sat_start, sat_end], [read.reference_start, read.reference_end])
        mini_list.append([rep_row['query_name'], rep_row['query_start'], rep_row['query_end'], rep_row['repeat_name'], read.query_name, read.reference_start, read.query_length, distance2repeat, gene_loc])
        # keep track of positions for filtering out later.
        positions.append(read.positions)

    # Do a quick filter to remove similar transcripts. First just by location, then by sequence?
    #keep = filter_similar_transcripts(positions)
    #keep_list = list()
    #for i in keep:
    #    keep_list.append( mini_list[i] )

    # concat to big list.
    #L.extend(keep_list)
    # Keep all matches.
    L.extend(mini_list)

# Convert list to df.
DATA=pd.DataFrame(L,columns=column_names)
# Save DATA to file.
DATA.to_csv('all_matches_10k.dat', sep='\t', index=False)

# Test reading/fetching from sam file.
#read = sat_bam.fetch("NC_000004.12",1,50000)
#for x in read:
#    print(str(x))


