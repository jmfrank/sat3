# generate a large table with all replicates of KIT vs PLPPR3

import os, sys, pysam, math
import pandas as pd
import numpy as np
from tqdm import tqdm

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from gen import gen

base_dir = '/media/ngs/data/GSE144085_sperm_stem_cells/'

# Sra ids.
kit_sras=["SRR10953071", "SRR10953072", "SRR10953073", "SRR10953074"]
PLPPR3_sras=["SRR10953067","SRR10953068", "SRR10953069", "SRR10953070"]

# array for est_counts.
est_counts=np.empty([163975,4])

# First loop over KIT+
for i, sra in enumerate(kit_sras):
    this_file = base_dir + sra + '/abundance.tsv'
    this_data = pd.read_csv(this_file, sep='\t')
    est_counts[:,i]=this_data['tpm']


out_data=this_data[['target_id', 'length']]

# Now add est_counts columns.
out_data[kit_sras]=est_counts

# Now loop over PLPPR3
est_counts=np.empty([163975,4])
for i, sra in enumerate(PLPPR3_sras):
    this_file = base_dir + sra + '/abundance.tsv'
    this_data = pd.read_csv(this_file, sep='\t')
    est_counts[:,i]=this_data['tpm']

# Now add est_counts columns.
out_data[PLPPR3_sras]=est_counts

# Get acc vs name. Same order as kallisto data.
acc_vs_name=pd.read_csv('~/Dropbox/code_bank/sat3/acc_vs_name.txt',header=None,delim_whitespace=True)
acc_vs_name.columns=['accession', 'name']
out_data.insert(1, 'name', acc_vs_name['name'])

#Now we need to reduce the data by looping over unique names and summing est_counts.
unique_names = out_data['name'].unique()
reduced = np.empty([unique_names.shape[0],8])

for i, n in enumerate(tqdm(unique_names)):
    # Sum over transcript variants.
    reduced[i,] = out_data.iloc[(out_data['name']==n).values, 3:11].sum().values

out_data_reduced=pd.DataFrame(unique_names)
out_data_reduced.columns=['name']
out_data_reduced[kit_sras+PLPPR3_sras]=reduced

out_data_reduced.to_csv(base_dir+'summary_table.csv',index=False)