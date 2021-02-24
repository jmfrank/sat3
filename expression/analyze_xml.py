
import pandas as pd


dat_file = '/Users/franklin/Dropbox/TEAD_paper/data/genomics/ncbi_expression/PRJEB4337_GRCh38.p7_108_expression.xml.gz_fixed.xml.dat'

data = pd.read_csv( dat_file, sep='\t')


# First organize data into tissues. Use metadata_ key to associate SRA number to tissue.
meta = data.loc[data['idis_metadata'].str.contains('metadata_9606_')]
# First, find unique tissues.
tissues = meta['sra_id'].unique()
# Loop over tissues. get associated sra numbers for each tissue.
tissue_dict = {}
for i in tissues:
        # Find associated sras.
        tissue_dict[i] = meta.loc[ meta['sra_id'] == tissues[0], 'source_name']

# Now, can we add standard deviation to data knowing tissue ids?

