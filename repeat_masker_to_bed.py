# Converting a repeat masker output to BED file format for visualization?

import os, pysam
import pandas as pd
def rm2bed(in_file, out_file):

    out_f = open(out_file,'w')

    with open(in_file) as f:
        line = f.readline()

        while line:
            S = line.split()

            cols = len(S)
            # custom list of desired outputs values.
            outL=list()
            outL.append(S[4])
            outL.append(S[5])
            outL.append(S[6])
            outL.append(S[9])
            out_f.write("\t".join(outL) + '\n')

            line = f.readline()

    out_f.close()


dir='/media/ngs/data/genomes/HG002_CCS_canu_paternal/'

#rm_file =dir+ 'sat3_repeats.out'
#rm_file_out = dir+ 'sat3_repeats.bed'

# Need to read line by line cause different number of columns due to *
#rm2bed(rm_file, rm_file_out)

rm_file = dir + 'HSATII_repeats.out'
rm_file_out = dir + 'HSATII_repeats.bed'
rm2bed(rm_file, rm_file_out)

rm_file = dir + 'ATTCC-related_repeats.out'
rm_file_out = dir + 'ATTCC-related_repeats.bed'
rm2bed(rm_file, rm_file_out)

