
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

# Read the processed repeat-transcript matches
DATA = pd.read_csv('all_matches_10k.dat', sep='\t')

# Now choose some particular candidates, extract the whole region.
index = 0
match = DATA.loc[index]

fasta = pysam.FastaFile(HG002_CCS_genome)
read = fasta.fetch(reference=match['contig'], start=match['ref_start'], end=match['ref_end'])
read



# Can we find an arbitrary transcript? Looking for ANKRD30A
transcript = "NM_052997.3"
#transcript = "NM_003430.4" # ZNF91
A = pysam.IndexedReads(rna_bam)
A.build()
#RNA = A.find('XM_011535254.3')
RNA = A.find(transcript)
for r in RNA:
    print(r)
    print(r.query_name)
    print(r.reference_name)

