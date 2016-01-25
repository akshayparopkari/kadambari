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
    import seaborn as sns
except ImportError:
    importerror.append("seaborn")
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
    cat_gramox_data.columns = ["category", "gramox", "counts",
                               "total_otus", "otu_pct"]

# Get colors for plot
    bmap = colorbrewer.qualitative.Accent_8
    brewer_colors = cycle(bmap.hex_colors)
    colors = [brewer_colors.next()
              for i in range(len(cat_gramox_data.category))]

# Create plot
    plt.figure(figsize=(20, 15))
    sns.set_style("ticks")
    ax = sns.barplot(x="gramox", y="otu_pct", hue="category",
                     data=cat_gramox_data,
                     hue_order=["Control", "Smokers", "Pregnant",
                                "Pregnant Smokers"],
                     palette=colors,
                     saturation=1)
    ax.set_xlabel("")
    ax.set_xticklabels(np.unique(cat_gramox_data.gramox), size=18, alpha=1)
    ax.set_ylabel("OTU Overlap Percentage (%)", labelpad=20, size=18, alpha=1)
    ax.set_yticklabels(range(0, int(max(cat_gramox_data.otu_pct) + 10), 5),
                       alpha=1, size=18)
    ax.legend(loc=0, fontsize=18, frameon=True)
    ax.grid()
    plt.show()

if __name__ == '__main__':
    main()
