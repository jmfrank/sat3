# Look at tissue-specific expression of genes nearby Sat3 repeats.

import os, pysam, math
from tqdm import tqdm
import pandas as pd
from gen import gen
import plotting
import numpy as np


# Full file.
rna_counts = pd.read_csv('/media/ngs/data/H1_hesc_seq_GSE127201/summary_table.csv',sep=',')

# Split into whitehead and crispri studies.
white_head = rna_counts[['name','SRR574820','SRR574821']]
white_head.to_csv('/media/ngs/data/H1_hesc_seq_GSE127201/summary_table_whitehead.csv',index=False)

crispri = rna_counts[['name','SRR8633134','SRR8633135']]
crispri.to_csv('/media/ngs/data/H1_hesc_seq_GSE127201/summary_table_CRISPRi_endoderm.csv',index=False)
