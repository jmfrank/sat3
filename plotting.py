import matplotlib.pyplot as plt
import numpy as np

from matplotlib import colors
from matplotlib.ticker import PercentFormatter


# Plot an arbitrary histogram.
def histplot( data, xlabel, bin_count=50, xlims=None, title=None, label_rotation=None ):


    fig, axs = plt.subplots(tight_layout=True)


    if xlims is not None:
        data = data[(data<=xlims[1]) & (data>=xlims[0])]

    axs.hist(data, bins=bin_count)
    axs.set_xlabel(xlabel)
    axs.set_ylabel('Counts')

    # xlimits.
    if xlims is not None:
        axs.set_xlim(xlims)

    # Rotate ticks. NEEDS FIX!
    if label_rotation is not None:
        axs.set_xticklabels(data.unique(), rotation=label_rotation, ha="right")


    axs.set_title(title)
    fig.show()

# Tissue expression bar chart. Send data as pd series?
def tissue_bar( data ):
    cats = data['tissue'].astype('string').values

    x_pos = np.arange(len(cats))


    fig, ax = plt.subplots()
    ax.bar(x_pos, data['mean'], yerr=data['std'],align='center')
    ax.set_ylabel('Mean RPKM')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(cats,rotation=45, ha="right")
    fig.show()

# Bar chart from dictionary
def bar( data ):

    cats = list(data.keys())

    x_pos = np.arange(len(cats))


    fig, ax = plt.subplots()
    ax.bar(x_pos, list(data.values()), align='center')
    ax.set_ylabel('counts')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(cats,rotation=45, ha="right")
    fig.show()