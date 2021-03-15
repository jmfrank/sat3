# generate a large table with all replicates of KIT vs PLPPR3

import os, sys, pysam, math, glob
import pandas as pd
import numpy as np
from tqdm import tqdm

import plotting
import utilities

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

base_dir = '/media/ngs/data/H1_diff_ucsd_reference/'

# Sra ids.
dirs = glob.glob(base_dir+'SRR*/')
sras = [x.split('/')[5] for x in dirs]

# array for est_counts.
est_counts=np.empty([163975, len(sras)])

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

#Now we need to reduce the data by looping over unique names and taking sum of transcripts.
unique_names = out_data['name'].unique()
reduced = np.empty([unique_names.shape[0], len(sras)])

# Create index mapping the rows(idx) of out_data to specific gene names as dictionary:
#name_dict={}
#for i, n in enumerate(tqdm(unique_names)):
#    name_dict[n] = out_data[ out_data['name']==n].index.values
#utilities.save_obj(name_dict,'/home/matt/Dropbox/code_bank/sat3/name_array.pkl')
name_dict=utilities.load_obj('/home/matt/Dropbox/code_bank/sat3/name_array.pkl')

for i, n in enumerate(tqdm(name_dict)):
    # Sum over transcript variants.
    reduced[i,] = out_data.iloc[name_dict[n], 3:3+len(sras)+1].sum().values


out_data_reduced=pd.DataFrame(unique_names)
out_data_reduced.columns=['name']
out_data_reduced[sras]=reduced

out_data_reduced.to_csv(base_dir+'summary_table_diff.csv', index=False)