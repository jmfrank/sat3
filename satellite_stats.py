import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from utilities import load_repeat_masker_data ,filter_pandas

# Define file system.
genomes_dir ='/media/ngs/data/genomes/'

# Let's look at the recently assembled genome
HG002_CCS_dir = 'HG002_CCS_canu_paternal/'
sat3_matches = genomes_dir + HG002_CCS_dir + 'sat3_repeats.out'


repeat_data = load_repeat_masker_data(sat3_matches)
#repeat_data = filter_pandas(repeat_data, length=[80000,float('inf')])

fig, axs = plt.subplots(1, 3, tight_layout=True)
fig.show()

axs[0].hist(repeat_data['length'], bins=50)
axs[0].set_yscale('log')
axs[0].set_xlabel('Repeat length')
axs[1].hist(repeat_data['score'], bins=50)
axs[1].set_yscale('log')
axs[1].set_xlabel('Repeat score')

# Color repeat by type.
names = repeat_data.repeat_name.unique()
L = repeat_data.shape[0]
labels=[]
for name in repeat_data.repeat_name:
    for i, n in enumerate(names):
        if name == n:
            labels.append(i)

axs[2].scatter(repeat_data['length'],repeat_data['score'], c=labels)
axs[2].set_xlabel('Repeat length')
axs[2].set_ylabel('Repeat score')