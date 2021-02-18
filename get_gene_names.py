

# Extract transcript names and then query them at NCBI.


import os
from subprocess import Popen, PIPE
from Bio import Entrez, SeqIO
import time
import numpy as np
import pandas as pd
from eutils import Client, QueryService


# Output from satellite analysis.
dat_file = '/media/ngs/py_projects/sat3/all_matches.dat'

df = pd.read_csv(dat_file, sep='\t')
#print(df)

# Collect the transcript name.
names=df['transcript_name'].unique()
print("Number of unique transcripts:" + str( len(names)))
# API key
key='c31652edceedd7db5b686afc53ff92f47e09'

#os.system("esearch -db 'nuccore' -query %s | efetch -format fasta" %names[0])
#A = Popen(["esearch -db nuccore -query %s" %names[0]], shell=True, stdout=PIPE)
#B = Popen(["efetch -format xml"], shell=True, stdin=A.stdout)
#C = Popen(["xtract -pattern "])
#out, err = B.communicate()

#define email for entrez login
db           = "nuccore"
Entrez.email = "mfranklin@health.ucsd.edu"
Entrez.api_key=key
batchSize    = 1
retmax       = 10**9
delay=5
# Loop over batches.
L = []

for start in range(0, len(names), batchSize):

    these_names = names[start:start+batchSize]

    # Accessions comma separated
    query = " ".join(these_names)
    print(query)

    handle = Entrez.esearch(db=db, term=query, retmax=retmax)
    giList = Entrez.read(handle)['IdList']

    L.append([these_names, giList])

df=pd.DataFrame(L,columns=['accession','GI'])
df.to_csv('GI.csv')
