# generate a large table with all replicates of KIT vs PLPPR3

import os, sys, pysam, math
import pandas as pd
import numpy as np
from tqdm import tqdm

import plotting

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

base_dir = '/media/ngs/data/10A_rna_seq_rRNA_depl_GSE150003/'

# Sra ids.
sras=["SRR11714175", "SRR11714176", "SRR11714177", "SRR11714178","SRR11714179","SRR11714180","SRR11714181","SRR11714182"]

# array for est_counts.
est_counts=np.empty([163975,len(sras)])

for i, sra in enumerate(sras):
    this_file = base_dir + sra + '/abundance.tsv'
    this_data = pd.read_csv(this_file, sep='\t')
    est_counts[:,i]=this_data['tpm']

out_data=this_data[['target_id', 'length']]

# Now add est_counts columns.
out_data[sras]=est_counts

# Get acc vs name. Same order as kallisto data.
acc_vs_name=pd.read_csv('~/Dropbox/code_bank/sat3/acc_vs_name.txt',header=None,delim_whitespace=True)
acc_vs_name.columns=['accession', 'name']
out_data.insert(1, 'name', acc_vs_name['name'])

#Now we need to reduce the data by looping over unique names and summing est_counts.
unique_names = out_data['name'].unique()
reduced = np.empty([unique_names.shape[0],8])

for i, n in enumerate(tqdm(unique_names)):
    # Sum over transcript variants.
    reduced[i,] = out_data.iloc[(out_data['name']==n).values, 3:3+len(sras)+1].sum().values

out_data_reduced=pd.DataFrame(unique_names)
out_data_reduced.columns=['name']
out_data_reduced[sras]=reduced

out_data_reduced.to_csv(base_dir+'summary_table.csv',index=False)

