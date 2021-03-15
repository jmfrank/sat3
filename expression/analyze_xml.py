# Parses NLB xml data from RNAseq collection into more friendly format. Further, calculates the STD of RPKM across individual SRAs. Outputs a simple dataframe(csv)
# with rows for each RPKM X TISSUE X GENE. Also, adds a column with gene symbol. This will be useful for other analysis.
import pandas as pd
import numpy as np
from tqdm import tqdm
#dat_file = '/Users/franklin/Dropbox/TEAD_paper/data/genomics/ncbi_expression/PRJEB4337_GRCh38.p7_108_expression.xml.gz_fixed.xml.dat'
#dat_file = '/media/ngs/data/ncbi_expression/PRJEB4337_GRCh38.p7_108_expression.xml.gz_fixed.xml.dat'
dat_file = '/media/ngs/data/development_project_PRJNA270632/PRJNA270632_GRCh38.p7_108_expression.xml.gz_fixed.xml.dat'
data = pd.read_csv(dat_file, sep='\t')

# First organize data into tissues. Use metadata_ key to associate SRA number to tissue.
meta = data.loc[data['idis_metadata'].str.contains('metadata_9606_')]
# First, find unique tissues.
tissues = meta['sra_id'].unique()
# Loop over tissues. get associated sra numbers for each tissue.
tissue_dict = {}
for i in tissues:
    # Find associated sras.
    tissue_dict[i] = meta.loc[meta['sra_id'] == i, 'source_name'].astype('str').values

# Find unique gene ids.
gene_ids = data['gene'].unique()
# Get rid of 'nan'
gene_ids = gene_ids[~np.isnan(gene_ids)].astype('int')

# output as tuple.
summary_data = []

# Now, can we add standard deviation to data knowing tissue ids. Loop over gene id, then tissue.
for g in tqdm(gene_ids):
    # Collect relevant data.
    g_data = data.loc[data['gene'] == g]
    # Loop for tissues.
    for t in tissue_dict:
        sras = tissue_dict[t]
        vals = np.empty(sras.shape)
        # Loop over sras.
        for i, s in enumerate(sras):
            vals[i] = g_data.loc[g_data['source_name'] == s, 'full_rpkm'].values
        summary_data.append((g, t, vals.mean(), vals.std()))

# Now make Dataframe from tuple.
out_data = pd.DataFrame(summary_data, columns=['gene', 'tissue', 'mean', 'std'])

# Add gene names?
all_names_file = '/media/ngs/data/ncbi_expression/gene_resultSHORT.txt'
names_data = pd.read_csv(all_names_file, sep='\t')
# Get homo sapiens only.
names_data = names_data.loc[names_data['tax_id'] == 9606]

# New column in out_data.
out_data['symbol'] = ''
# Fill in symbol from names_data. Loop over unique names so we don't repeat for each tissue.
for g in tqdm(gene_ids):
    # See if there's a match.
    entry = names_data.loc[names_data['GeneID'] == g]
    if entry.shape[0] == 0:
        out_data.loc[out_data['gene'] == g, 'symbol'] = 'MISSING'
        print('missing geneID: ' + str(g))
        continue
    # Find all matching entries in out_data.
    out_data.loc[out_data['gene'] == g, 'symbol'] = entry['Symbol'].values[0]

# Save summary.
out_data.to_csv('/media/ngs/data/development_project_PRJNA270632/summary_PRJNA270632.csv', sep='\t', index=False)
