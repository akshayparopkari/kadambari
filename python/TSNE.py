#!/usr/bin/env python
"""
:Abstract: Visualize oral microbiome using T-SNE algorithm. From scikit-learn's help page,
           t-SNE [1] is a tool to visualize high-dimensional data. It converts
           similarities between data points to joint probabilities and tries to minimize
           the Kullback-Leibler divergence between the joint probabilities of the low
           dimensional embedding and the high-dimensional data.
:Date: 08/29/2016
:Author: Akshay Paropkari
"""

import sys
import argparse
from phylotoast import util, graph_util as gu
errors = []
try:
    import matplotlib as mpl
    mpl.rc("font", family="Arial")  # define font for figure text
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d  # for annotating 3D plot
except ImportError as ie:
    errors.append(ie)
try:
    import pandas as pd
except ImportError as ie:
    errors.append(ie)
try:
    from sklearn.manifold import TSNE, MDS
except ImportError as ie:
    errors.append(ie)
if len(errors) != 0:
    for item in errors:
        print("Import Error:", item)
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description="This script calculates and returns TSNE"
                                                 " plots based on distance matrices "
                                                 "(for e.g. unifrac distance matrix).")
    parser.add_argument("-dm", "--dist_matrix_file", required=True,
                        help="Input distance matrix file. [REQUIRED]")
    parser.add_argument("-m", "--map_fp", required=True,
                        help="Metadata mapping file. [REQUIRED]")
    parser.add_argument("-g", "--group_by", required=True,
                        help="A column name in the mapping file containing\
                              categorical values that will be used to identify \
                              groups. Each sample ID must have a group entry. \
                              Default is no categories and all the data will be \
                              treated as a single group. [REQUIRED]")
    parser.add_argument("-c", "--color_by", required=True,
                        help="A column name in the mapping file containing\
                              hexadecimal (#FF0000) color values that will\
                              be used to color the groups. Each sample ID must\
                              have a color entry. [REQUIRED]")
    parser.add_argument("--plot_title", default=None,
                        help="Plot title. Default is no title.")
    parser.add_argument("-o", "--out_fp", default="",
                        help="The path and file name to save the LDA plot under.\
                              If specified, the figure will be saved directly\
                              instead of opening a window in which the plot \
                              can be viewed before saving")
    parser.add_argument("-d", "--dimensions", default=2, type=int, choices=[2, 3],
                        help="Choose whether to plot 2D or 3D.")
    parser.add_argument("--z_angles", type=float, nargs=2, default=[45., 30.],
                        help="Specify the azimuth and elevation angles for a 3D plot.")
    parser.add_argument("--figsize", default=[14, 8], type=int, nargs=2,
                        help="Specify the 'width height' in inches for LDA plots."
                             "By default, figure size is 14x8 inches.")
    parser.add_argument("--font_size", default=12, type=int,
                        help="Sets the font size for text elements in the plot.")
    parser.add_argument("--label_padding", default=15, type=int,
                        help="Sets the spacing in points between the each axis and its \
                             label.")
    parser.add_argument("--annotate", action="store_true",
                        help="If specified, each data point will be labeled with "
                             "its sample ID. Default is False.")
    parser.add_argument("--ggplot2_style", action="store_true",
                        help="Apply ggplot2 styling to the figure. Default is False.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Read in the distance data
    try:
        dm_data = pd.read_csv(args.dist_matrix_file, sep="\t", index_col=0)
    except IOError as ioe:
        sys.exit("\nError reading in distance matrix file: {}.".format(ioe))

    # mapping and colors info for plotting
    try:
        header, map_data = util.parse_map_file(args.map_fp)
    except IOError as ioe:
        sys.exit("\nError reading mapping file: {}.".format(ioe))
    y = [map_data[sid][header.index(args.group_by)] for sid in dm_data.index]
    cond_colors = util.color_mapping(map_data, header, args.group_by, args.color_by)

    # Prep input data for t-SNE
    X = dm_data[range(dm_data.shape[1])].values
    X_tsne = TSNE(n_components=3, metric="precomputed").fit_transform(X)

    # Plot t-SNE result
    fig = plt.figure(figsize=(14, 8))
    for cond, sid, xy in zip(y, dm_data.index, X_tsne):
        plt.scatter(x=xy[0], y=xy[1], s=150, c=cond_colors[cond], alpha=0.85,
                    edgecolors="k")
        if args.annotate:
            plt.annotate(s=sid, xy=(xy[0], xy[1]), xytext=(12, 12),
                         textcoords="offset points", ha="center", va="center",
                         alpha=1, style="italic")
    if args.plot_title is not None:
        plt.title(args.plot_title, fontsize=16, weight="bold")
    l = [plt.scatter([], [], c=cond_colors[cond], s=150, edgecolors="k")
         for cond in cond_colors]
    plt.legend(l, ["{}".format(cond) for cond in cond_colors], loc="best",
               scatterpoints=3, frameon=True, framealpha=1, fontsize=14)
    plt.xlabel("t-SNE 1", fontsize=16)
    plt.ylabel("t-SNE 2", fontsize=16)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.grid()
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
