

import os
from subprocess import Popen, PIPE
from Bio import Entrez, SeqIO
import time
import numpy as np
import pandas as pd


# Output from satellite analysis.
acc_list = '/media/ngs/py_projects/sat3/GI.csv'

df = pd.read_csv(acc_list)
gi_list=open('gi_list.csv', 'w+')
# Get GIs at int.
I = []
for i, l in enumerate(df.GI):

    I.append(int(l.strip("['").strip("']")))
    gi_list.write( l.strip("['").strip("']") + "\n")

gi_list.close()

# Output from satellite analysis.
dat_file = '/media/ngs/py_projects/sat3/all_matches.dat'

df = pd.read_csv(dat_file, sep='\t')
#print(df)

# Collect the transcript name.
names=df['transcript_name'].unique()
acc_list=open('acc_list.csv','w+')
for i, l in enumerate(names):

    acc_list.write( l + "\n")

acc_list.close()