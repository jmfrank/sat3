# Converting a repeat masker output to BED file format for visualization?

import os, pysam
import pandas as pd

def meme2bed(in_file, out_file):

    out_f = open(out_file,'w')

    with open(in_file) as f:
        line = f.readline()

        #Line count
        L=0
        while line:
            if L > 0:
                S = line.split()

                cols = len(S)

                # check that format is correct.

                # custom list of desired outputs values.
                outL=list()
                # Chromosome
                outL.append(S[1])
                # Start
                outL.append(S[2])
                # Stop
                outL.append(S[3])
                # Feature (0,1,2,3)
                outL.append("0")
                # Size [0,1]
                outL.append("1")
                # Color
                outL.append("##FF0000")
                # Chr copy (1,2)
                outL.append("1")

                out_f.write("\t".join(outL) + '\n')

            L = L+1
            line = f.readline()

    out_f.close()


dir='/media/ngs/data/genomes/chm13_v1.0/'

meme_file = dir + 'CATTCC_fimo.out'
bed_file = dir + 'CATTCC_fimo.bed'
meme2bed(meme_file, bed_file)
