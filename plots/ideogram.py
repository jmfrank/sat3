"""
Demonstrates plotting chromosome ideograms and genes (or any features, really)
using matplotlib.
1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
saved as "ideogram.txt". Lines look like this::
    #chrom  chromStart  chromEnd  name    gieStain
    chr1    0           2300000   p36.33  gneg
    chr1    2300000     5300000   p36.32  gpos25
    chr1    5300000     7100000   p36.31  gneg
2) Assumes another file, "ucsc_genes.txt", which is a BED format file
   downloaded from UCSC's Table Browser. This script will work with any
   BED-format file.
"""

import matplotlib, pandas
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection

# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print(chrom)
        # Perform colormap as needed.
        norm_vals = (group['count'] / ideo['count'].max())
        # Log-norm.
        #norm_vals = np.log10(group['count'])/ np.log10(ideo['count'].max())
        #norm_vals[ norm_vals == -np.inf] = 0
        COLORS = get_color_map(norm_vals,'hot')
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values

        # Create rectangle that surrounds chromosome.
        rec = matplotlib.patches.Rectangle([xranges[0,0], yrange[0]], xranges[-1,0], height, angle=0.0,facecolor=None, fill=False)
        yield (BrokenBarHCollection(
            xranges, yrange, facecolors=COLORS, **kwargs), rec)


    if del_width:
        del df['width']

# Colormap.
def get_color_map( vals, cmap_name):
    # Colormap.
    cmap = matplotlib.cm.get_cmap(cmap_name)
    COLORS = cmap(vals)
    COLORS[:,3]=vals
    COLORS[:,0]=1
    COLORS[:,1]=0
    COLORS[:,2]=0
    return COLORS

# Height of each ideogram
chrom_height = 1

# Spacing between consecutive ideograms
chrom_spacing = 1

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 0.4

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.1

# Width, height (in inches)
figsize = (6, 8)

# Decide which chromosomes to use
chromosome_list = ['chr%s' % i for i in range(1, 23)] + ['M', 'chrX', 'chrY']

# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase - gene_height - gene_padding
    ybase += chrom_height + chrom_spacing

# Read in BED data.
bed_file = '/media/ngs/data/genomes/chm13_v1.0/binned_CATTCC_counts.bed'

ideo = pandas.read_table(
    bed_file,
    names=['chrom', 'start', 'end', 'count']
)
# Replace periods with 0's.
ideo['count'].replace({".": "0"}, inplace=True)
# Now convert column to integer.
ideo['count'] = ideo['count'].astype('int16')

# Filter out chromosomes not in our list
ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

# Add a new column for width
ideo['width'] = ideo.end - ideo.start

fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)
# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
    ax.add_collection(collection[0])
    ax.add_patch(collection[1])


# Axes tweaking
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list)
ax.axis('tight')
plt.show()