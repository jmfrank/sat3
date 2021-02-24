import matplotlib.pyplot as plt

from matplotlib import colors
from matplotlib.ticker import PercentFormatter


# Plot an arbitrary histogram.
def histplot( data, xlabel, bin_count=50, title=None ):

    fig, axs = plt.subplots(tight_layout=True)

    axs.hist(data, bins=bin_count)
    axs.set_xlabel(xlabel)
    axs.set_ylabel('Counts')
    axs.set_title(title)
    fig.show()