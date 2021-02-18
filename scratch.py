import pysam

F = pysam.AlignmentFile("/media/ngs/data/genomes/HG002_CCS_canu_paternal/RNA_aligned_to_HG002.bam")
ref_contig='SRHB01000005.1'
reads = F.fetch(reference=ref_contig)

# Get first of iterator.
r = next(reads)

# Indexing reads.
A = pysam.IndexedReads(F)
A.build()
D = A.find('XM_011535254.3')

# Get the total coverage from fasta?
fasta = pysam.FastaFile("/media/ngs/data/genomes/HG002_CCS_canu_paternal/GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.fna")

# Find introns?
introns = F.find_introns(D)

# Loop over introns, spit them out into multi-fasta.
intron_fasta=open('intron_fasta.fasta','w+')
count=0
for I in list(introns.keys()):
    # I is a tuple. Extract genome region.
    this_intron_seq = fasta.fetch(reference=ref_contig, start=I[0], end=I[1])
    # create line in fasta.
    intron_fasta.write( ">"+str(count)+" \n"+ this_intron_seq +" \n")
    count = count + 1

intron_fasta.close()

P = next(D)


total_gene = fasta.fetch(reference=ref_contig, start=P.reference_start, end=P.reference_end)

save_fasta = open('temp_fasta.fasta','w+')
DNA = ">SRHB01000005.1, start = 3079415, end=3210794 \n" + total_gene

save_fasta.write(DNA)
save_fasta.close()

import pysam

# Look for a transcript in fasta.
rna_fasta = pysam.FastaFile("/media/ngs/data/genomes/GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna")