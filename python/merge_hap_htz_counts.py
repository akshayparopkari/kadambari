#!/usr/bin/env python
"""
:Abstract: Merge output of heterozygosity count calculations and plot all female
           individuals for all autosomes.
:Date: 05/23/2018
:Author: Akshay Paropkari
"""

import sys
import gzip
import argparse
from os import walk
from time import strftime
from os.path import join
err = []
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc("font", family="Arial")
    mpl.rc("xtick", labelsize=11)  # set X axis ticksize
    mpl.rc("ytick", labelsize=11)  # set Y axis ticksize
except ImportError:
    err.append("matplotlib")
try:
    from palettable.colorbrewer.sequential import Blues_9    # Europe
    from palettable.colorbrewer.sequential import Greens_9   # East Asia
    from palettable.colorbrewer.sequential import Purples_9  # South Asia
    from palettable.colorbrewer.sequential import Reds_9     # Americas
    from palettable.colorbrewer.sequential import YlOrBr_9   # Sub-Saharan Africa
except ImportError:
    err.append("palettable")
try:
    import pandas as pd
except ImportError:
    err.append("pandas")
try:
    assert len(err) == 0
except AssertionError:
    for error in err:
        sys.exit("Please install {}".format(error))


def handle_program_options():
    parser = argparse.ArgumentParser(description="Merge output of heterozygosity count "
                                     "calculations and plot all female individuals for "
                                     "all autosomes.")
    parser.add_argument("input_directory", help="Path to the folder which only contains"
                        " the output of heterozygosity counts. [REQUIRED]")
    parser.add_argument("metadata", help="Metadata mapping file for 1000 genomes project")
    parser.add_argument("-p", "--plot_data", action = "store_true",
                        help="Supply this parameter to plot all autosome heterozygosity "
                        "counts for each population.")
    parser.add_argument("-s", "--savefile", help="Save the plot as an SVG file. Provide "
                        "file path and file name with extension.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Iterate through all files and get combined heterozygosity results for autosomes
    for root, folders, files in walk(args.input_directory):
        all_data = pd.DataFrame()
        for f in files:
            if "chrX" not in f:
                try:
                    input = pd.read_csv(join(root, f), sep="\t", index_col=None,
                                        usecols= ["sample", "variant_counts"])
                except Exception as ee:
                    sys.exit("Error while reading heterozygosity counts files\n{}".
                             format(ee))
                else:
                    all_data = pd.concat([all_data, input], ignore_index=True,
                                         verify_integrity=True)
            else:
                all_data = pd.read_csv(join(root, f), sep="\t", index_col=None)
                all_data = all_data.query("variant_counts != 0").dropna()
                break

    try:
        assert "gender" not in all_data.columns.values
    except AssertionError:
        all_data = all_data[["sample", "variant_counts", "pop", "super_pop", "gender"]]
        all_data.to_csv("fem_chrX_htz_counts.txt", sep="\t", index=False)
    else:
        all_data = all_data.groupby(["sample"]).sum()
        md_data = pd.read_csv(args.metadata, sep="\t", index_col=False,
                              usecols=[0, 1, 2, 3])
        all_data = pd.merge(all_data, md_data, on="sample")
        all_data.to_csv("all_autosome_htz_counts.txt", sep="\t", index=False)

    # Plot heterozygosity counts
    if args.plot_data:
        colors = Blues_9.hex_colors[4:] + Greens_9.hex_colors[4:] +\
                 Purples_9.hex_colors[4:] + Reds_9.hex_colors[2:6] +\
                 YlOrBr_9.hex_colors[2:]
        x_order = tuple(["FIN", "GBR", "CEU", "IBS", "TSI", "CHS", "CDX", "CHB",
                         "JPT", "KHV", "GIH", "STU", "PJL", "ITU", "BEB", "PEL",
                         "MXL", "CLM", "PUR", "ASW", "ACB", "GWD", "YRI", "LWK",
                         "ESN", "MSL"])
        assigned_colors = {a: b for a, b in zip(x_order, colors)}
        x_order = all_data.groupby(["pop"]).mean().sort_values(by="variant_counts")
        with mpl.style.context("seaborn-white"):
            plt.figure(figsize=(12, 6))
            for entry in x_order.itertuples():
                plt.scatter(entry[0], entry[1], marker="_", s=300, c="k")
            for entry in all_data.itertuples():
                try:
                    assert entry[5] == "female"
                except AssertionError:
                    continue
                else:
                    plt.scatter(entry[3], entry[2], marker=".", edgecolors="face",
                                c=assigned_colors[entry[3]])
            plt.title("Heterozygosity in female individuals' chromosome X", fontsize=12)
            plt.grid(axis="y", linestyle=":")
            plt.tight_layout()
            if args.savefile:
                plt.savefig(args.savefile, dpi=300, format="svg",
                            bbox_inches="tight", pad_inches=0.25)
            else:
                plt.show()


if __name__ == "__main__":
    sys.exit(main())
