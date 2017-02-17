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
import itertools
from phylotoast import util, graph_util as gu
errors = []
try:
    import matplotlib as mpl
    mpl.rc("font", family="Arial")  # define font for figure text
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d  # for annotating 3D plot
except ImportError as ie:
    errors.append("matplotlib")
try:
    from palettable.colorbrewer.qualitative import Set1_9
except ImportError as ie:
    errors.append("palettable")
try:
    import pandas as pd
except ImportError as ie:
    errors.append("pandas")
try:
    import numpy as np
except ImportError as ie:
    errors.append("numpy")
try:
    from sklearn.metrics.pairwise import pairwise_distances
    from sklearn.manifold import TSNE
except ImportError as ie:
    errors.append("scikit-learn")
if len(errors) > 0:
    for item in errors:
        print("Import Error: {}".format(item))
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
                        help="A column name in the mapping file containing categorical "
                             "values that will be used to identify groups. Each sample ID"
                             " must have a group entry. Default is no categories and all "
                             "the data will be treated as a single group. [REQUIRED]")
    parser.add_argument("-c", "--color_by",
                        help="A column name in the mapping file containing hexadecimal "
                             "(#FF0000) color values that will be used to color the "
                             "groups. Each sample ID must have a color entry.")
    parser.add_argument("-p", "--perplexity", default=25, type=int,
                        help="According to Laurens van der Maaten, perplexity is a "
                             "measure for information that is defined as 2 to the power "
                             "of the Shannon entropy. In t-SNE, the perplexity may be "
                             "viewed as a knob that sets the number of effective nearest "
                             "neighbors. It is comparable with the number of nearest "
                             "neighbors k that is employed in many manifold learners. "
                             "Please see http://lvdmaaten.github.io/tsne/ for more "
                             "details. Default value is set to 25. Ideally this value "
                             "should be less than the number of samples in distance "
                             "matrix file.")
    parser.add_argument("-t", "--title", default="", help="Title of the plot.")
    parser.add_argument("-s", "--point_size", default=100, type=int,
                        help="Specify the size of the circles representing each of the "
                             "samples in the plot")
    parser.add_argument("--plot_title", default=None,
                        help="Plot title. Default is no title.")
    parser.add_argument("-o", "--out_fp", default="",
                        help="The path and file name to save the LDA plot under.If "
                             "specified, the figure will be saved directly instead of "
                             "opening a window in which the plot can be viewed before "
                             "saving.")
    parser.add_argument("--figsize", default=[14, 8], type=int, nargs=2,
                        help="Specify the 'width height' in inches for LDA plots. By "
                             "default, figure size is 14x8 inches.")
    parser.add_argument("--font_size", default=12, type=int,
                        help="Sets the font size for text elements in the plot.")
    parser.add_argument("--label_padding", default=10, type=int,
                        help="Sets the spacing in points between the each axis and its "
                             "label.")
    parser.add_argument("--annotate", action="store_true",
                        help="If specified, each data point will be labeled with its "
                             "sample ID. Default is False.")
    parser.add_argument("--ggplot2_style", action="store_true",
                        help="Apply ggplot2 styling to the figure. Default is False.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Read in the distance data
    try:
        dm_data = pd.read_csv(args.dist_matrix_file, sep="\t", index_col=0)
        dm_data_sids = dm_data.index
        dm_data = pairwise_distances(dm_data[range(dm_data.shape[1])].values,
                                     metric="precomputed")
    except IOError as ioe:
        sys.exit("\nError reading in distance matrix file: {}.".format(ioe))

    # Mapping and colors info for plotting
    try:
        header, map_data = util.parse_map_file(args.map_fp)
    except IOError as ioe:
        sys.exit("\nError reading mapping file: {}.".format(ioe))
    y = [map_data[sid][header.index(args.group_by)] for sid in dm_data_sids]

    # Get colors for all categories
    if not args.color_by:
        categories = set(y)
        bcolors = itertools.cycle(Set1_9.hex_colors)
        cond_colors = {c: bcolors.next() for c in categories}
    else:
        cond_colors = util.color_mapping(map_data, header, args.group_by, args.color_by)

    # Prep input data for t-SNE
    X_tsne = TSNE(n_components=3, perplexity=args.perplexity, metric="precomputed",
                  method="exact", verbose=2, random_state=0, angle=0.8)
    X_new = X_tsne.fit_transform(dm_data)
    print("KL divergence after optimization: {}\n".format(X_tsne.kl_divergence_))
    x_min, x_max = np.min(X_new, 0), np.max(X_new, 0)
    X_new = (X_new - x_min) / (x_max - x_min)

    # Plot t-SNE result
    fig = plt.figure(figsize=(14, 8))
    for cond, sid, xy in zip(y, dm_data_sids, X_new):
        ax = fig.add_subplot(111)
        ax.scatter(x=xy[0], y=xy[1], s=args.point_size, c=cond_colors[cond],
                   alpha=0.9, edgecolors="k")
        if args.annotate:
            ax.annotate(s=sid, xy=(xy[0], xy[1]), xytext=(12, 12),
                        textcoords="offset points", ha="center", va="center",
                        alpha=1, style="italic")
    if args.plot_title is not None:
        ax.set_title(args.plot_title, fontsize=16, weight="bold")
    l = [plt.scatter([], [], c=cond_colors[cond], s=150, edgecolors="k")
         for cond in cond_colors]
    plt.legend(l, ["{}".format(cond) for cond in cond_colors], loc="best",
               scatterpoints=3, frameon=True, framealpha=1, fontsize=14)
    ax.set_xlabel("t-SNE 1", fontsize=14)
    ax.set_ylabel("t-SNE 2", fontsize=14)
    plt.tight_layout()
    if args.ggplot2_style:
        gu.ggplot2_style(ax)
        fc = "0.8"
    else:
        fc = "none"

    # save or display result
    if args.out_fp:
        plt.savefig(args.out_fp, facecolor=fc, edgecolor="none", dpi=300, pad_inches=0.1,
                    bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
