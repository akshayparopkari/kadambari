#! /usr/bin/env python
"""
Abstract: Create MA plots from DESeq results.
Author: Akshay Paropkari
Date: 08/24/2016
"""
import sys
import math
import argparse
from random import uniform
from collections import Counter
importerrors = []
try:
    import pandas as pd
except:
    importerrors.append("pandas")
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc("font", family="Arial")  # define font for figure text
    mpl.rc("xtick", labelsize=14)  # increase X axis ticksize
    mpl.rc("ytick", labelsize=14)  # increase Y axis ticksize
except:
    importerrors.append("matplotlib")
if len(importerrors) > 0:
    for err in importerrors:
        print("Please install module: {}".format(err))
    sys.exit()


# Get data for plotting
def get_plot_data(deseq_res):
    """
    From DESeq2 results, obtain formatted data for MA plots.

    :type deseq_res: file handle or file path
    :param deseq_res: DESEq results file
    """
    in_data = pd.read_csv(deseq_res, sep="\t")
    x_min, x_max = min(in_data["avgLogExpr"]), max(in_data["avgLogExpr"])
    coord = {"x": [], "y": [], "color": [], "label": []}
    for rows in in_data.iterrows():
        row = rows[1]
        try:
            # Check if row represents unique or shared gene
            # Unique genes will fail the assertion test
            assert math.isnan(row["avgLogExpr"]) == False
        except AssertionError:
            # If avgLogExpr not given (unique genes), pick a value
            # randomly between min and max avgLogExpr value
            coord["x"].append(uniform(x_min, x_max))
            coord["color"].append("#0000FF")
            coord["y"].append(row["rLogFC"])
            coord["label"].append("Unique Genes")
        else:
            # Common genes
            coord["x"].append(row["avgLogExpr"])
            coord["y"].append(row["rLogFC"])
            if row["rLogFC"] > 2:
                coord["color"].append("#FF0000")
                coord["label"].append("Shared & >2 Log$_2$Fold Genes")
            elif row["rLogFC"] < -2:
                coord["color"].append("#FF0000")
                coord["label"].append("Shared & >2 Log$_2$Fold Genes")
            else:
                coord["color"].append("#808080")
                coord["label"].append("Shared Genes")
    return coord


def handle_program_options():
    """Command line arguments."""
    parser = argparse.ArgumentParser(description="Create MA plots from DESeq results. "
                                     "There are no checks for "
                                     "significant results, all logfold changes are "
                                     "plotted. Shared genes, unique genes and shared "
                                     "genes greater than 2 log fold difference are "
                                     "highlighted in output.")
    parser.add_argument("deseq_result",
                        help="Tab-separated DESEq results file. [REQUIRED]")
    parser.add_argument("-t", "--plot_title", default=None,
                        help="Plot title. Default is no title.")
    parser.add_argument("-o", "--out_filepath", default="",
                        help="The path and file name to save the MA plot under.If "
                             "specified, the figure will be saved directly instead of "
                             "opening a window in which the plot can be viewed before "
                             "saving.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Get input data
    plot_data = get_plot_data(args.deseq_result)
    counts = Counter(plot_data["label"])
    print("\nUnique Genes: {}\nShared Genes: {}\nShared & >2 LFC Genes: {}\n".
          format(counts["Unique Genes"], counts["Shared Genes"],
          counts["Shared & >2 Log$_2$Fold Genes"]))

    # Plot input data
    plt.figure(figsize=(12, 10))
    plt.scatter(x=plot_data["x"], y=plot_data["y"], c=plot_data["color"],
                s=20, alpha=1, edgecolors="none")
    colors = {"Unique Genes": "#0000FF", "Shared Genes": "#808080",
              "More than 2 Log$_2$Fold": "#FF0000"}

    # Manually add horizontal legend on top
    c = [plt.scatter([], [], c=colors[s1], edgecolors="none", s=125, alpha=1)
         for s1 in ["Unique Genes", "Shared Genes", "More than 2 Log$_2$Fold"]]
    plt.legend(c, ["{}".format(s2) for s2 in ["Unique Genes", "Shared Genes",
                                              "Shared & >2 Log$_2$Fold Genes"]],
               fontsize=14, scatterpoints=3, loc=2, ncol=3, mode="expand",
               bbox_to_anchor=(0., 0.96, 1., 0.1), borderaxespad=0., frameon=True)

    plt.xlabel("avgLogExpr", fontsize=14)
    plt.ylabel("LogFC", fontsize=14)
    plt.axhline(linewidth=2, color='g', alpha=1)
    if args.plot_title:
        plt.title(args.plot_title, fontsize="xx-large", weight="bold", x=0.5, y=1.075)

    # Save plot file (optional)
    if args.out_filepath:
        plt.savefig(args.out_filepath, edgecolor="none", dpi=300, bbox_inches="tight",
                    pad_inches=0.1)
    else:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
