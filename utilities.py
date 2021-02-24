# Utilities for analyzing genomic datas.
# Scratch script for getting started with genomics.

import os
import pysam
import pandas as pd
import numpy as np

# Get access to fasta. Auto-creates fai index.
def get_fasta_faidx( in_file):

    print(in_file)

    F = pysam.FastaFile(in_file)

    return F

# Find repeats that match substring(s)
def extract_matching_strings(in_file, out_file, strings):

    out_f = open(out_file,'w')

    with open(in_file) as f:
        line = f.readline()
        while line:
            if any(x in line for x in strings):
                out_f.write(line)
            line = f.readline()

    out_f.close()

# Generate an empty dataframe to describe transcripts.
def generate_rna_dataframe():
    column_names = ['transcript_name','start_position','length','distance_to_repeat','gene_loq','seq']
    data = pd.DataFrame(columns=column_names)

# Load repeat masker-style data.
def load_repeat_masker_data(in_file):
    data = pd.read_csv( in_file, header=None, delim_whitespace=True, usecols=range(15))
    data.columns = ['score', '%sub','%del','%ins','query_name','query_start','query_end','num_past_match','strand', \
                    'repeat_name','class','bases_before_match','match_start','match_end', 'linkageID']
    # calculate lengths.
    data["length"] = data["query_end"] - data["query_start"]

    return data

# Load simple bam file. *** might need to change column definition.
def load_bin_bam(in_file):
    data = pd.read_csv( in_file, header=None, delim_whitespace=True)
    data.columns = ['contig','start','stop','count']

    # Replace periods with 0.
    data.loc[data['count']=='.','count']=0

    # Convert columns to integers.
    data = data.astype({'contig': str, 'start': int, 'stop': int, 'count': int})

    return data

# Filter dataframe by column values. All filters are specified [min,max] range.
def filter_pandas(data, **filters):

    # Now iterate over keyword args **filters and filter the data as suggested.
    for column in filters:

        # Check if column exists.
        if column not in data.columns:
            print('Column: '+column+' not found. Skipping')
            continue

        vals = filters[column]
        sel = (data[column] >= vals[0]) & (data[column] <= vals[1])
        data = data[sel]

    return data

def index_bam_file( bam_file):

    #bam = pysam.AlignmentFile(bam_file,"rb")
    pysam.index(bam_file)

def convert_sam_2_bam( sam_file ):

    sam = pysam.AlignmentFile(sam_file,"r")
    bam_file = os.path.splitext(sam_file)[0]+'.bam'
    bam = pysam.AlignmentFile(bam_file, "wb", header=sam.header)

    for read in sam.fetch():
        bam.write(read)

    sam.close()
    bam.close()
    print("Converted sam to bam")
    pysam.index(bam_file)


def average(lst):
    return sum(lst) / len(lst)


# Calculate nearest distance between two genome elements.
def nearest_distance( seq_ref, seq_query ):

    # First see if they overlap.
    ref_range = set(range( seq_ref[0], seq_ref[1]+1))
    query_range = set(range( seq_query[0], seq_query[1] + 1))
    overlap = ref_range & query_range

    if len( overlap ) > 0:

        gene_loc = 'internal'

        # Distance is difference in mean position.
        distance = abs( average(seq_query) - average((seq_ref)))

    else:

        # Now check if query starts before or after reference.
        if seq_query[0] > seq_ref[1]:

            distance = seq_query[0] - seq_ref[1]
            gene_loc = 'downstream'
        elif seq_query[1] < seq_ref[0]:

            distance = seq_ref[0] - seq_query[1]
            gene_loc = 'upstream'



    return distance, gene_loc

# Filter similar transcripts by either start/end, coverage, or by sequence? Input is an iterable from pysam fetch function.
def filter_similar_transcripts( indices, cutoff=0.2):

    # overlap matrix.
    overlap = np.zeros([len(indices),len(indices)])
    # Now loop over the list and remove entries that are similar by some fraction, 'cutoff'
    for i in range(len(indices)-1):
        for j in range(i+1, len(indices)):

            # Compute overlap of indices.
            overlap[i,j] = len( set(indices[i]) & set(indices[j]) ) / len(indices[i])

    # All lengths
    lengths=np.empty(len(indices))
    for i in range(len(indices)):
        lengths[i] = len(indices[i])

    # Now iterate and remove re-dundent transcripts.
    non_redundant = []
    redundant = []
    for i in range(len(indices)):

        # Skip if already marked as similar.
        if i in redundant or i in non_redundant:
            continue

        sel = np.where( overlap[i,:] >= cutoff )[0]

        # add current idx.
        sel = np.append(sel,i)

        # Find maximum length of these.
        max_idx = sel[ np.argmax( lengths[sel] )]

        non_redundant.append( max_idx )
        redundant.extend( sel[np.where( sel != max_idx )].tolist())

    return non_redundant



