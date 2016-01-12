#!/usr/bin/env python
"""
Abstract: Plot categorized gramox data. Currently, this script plots only 4
          data groups. Special thanks to Peter Kerpedjiev, whose blog post
          'Creating a Grouped Bar Chart in Matplotlib' helped clarify a lot.
          Here is the link to it -
          http://emptypipes.org/2013/11/09/matplotlib-multicategory-barchart/

Date: 09/24/2015

Author: Akshay Paropkari
"""

import sys
import argparse
from itertools import cycle
importerror = []
try:
    import numpy as np
except ImportError:
    importerror.append("numpy")
try:
    import pandas as pd
except ImportError:
    importerror.append("pandas")
try:
    import matplotlib.pyplot as plt
except ImportError:
    importerror.append("matplotlib")
try:
    from palettable import colorbrewer
except ImportError:
    importerror.append("palettable")
if len(importerror) != 0:
    for module in importerror:
        print "Please install {} package.".format(module)
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description='Plot categorized gramox '
                                     'data.')
    parser.add_argument('gramox_fnh',
                        help='File with categorized gramox count data.')
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.gramox_fnh):
            pass
    except IOError as ioe:
        err_msg = '\nError opening gramox count file: {}\n'
        sys.exit(err_msg.format(ioe))

# Read gramox data
    cat_gramox_data = pd.read_csv(args.gramox_fnh, sep='\t')
    gdata = np.array(cat_gramox_data)

# Data to plot
    conditions = np.unique(gdata[:, 0])
    categories = np.unique(gdata[:, 1])
    width = (1 - 0.3) / (len(conditions))

# Get colors for plot
    bmap = colorbrewer.qualitative.Accent_8
    brewer_colors = cycle(bmap.hex_colors)
    colors = [brewer_colors.next() for i in range(len(conditions))]

# Create plot
    fig = plt.figure(figsize={30, 15})
    ax = fig.add_subplot(111)
    for i, cond in enumerate(conditions):
        ind = range(1, len(categories) + 1)
        vals = gdata[gdata[:, 0] == cond][:, 2].astype(np.float)
        pos = [j - (1 - 0.3) / 2. + i * width for j in ind]
        ax.bar(pos, vals, width=width, label=cond, color=colors[i])
        for x, y in zip(pos, vals):
            ax.text(x+(width/len(conditions)), y+1, int(y), fontsize=14)
    ax.set_ylabel("Counts", size=24)
    ax.set_ylim(top=max(gdata[:, 2]) + 15)
    ax.set_xticks(ind)
    ax.set_xticklabels(categories, size=22)
    ax.legend(loc='best', fontsize=22)
    plt.show()

if __name__ == '__main__':
    main()
