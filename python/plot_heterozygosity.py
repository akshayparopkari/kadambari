#!/usr/bin/env python
"""
:Abstract: Plot number of non-reference sites per genome from 1000 genome project.
:Date: 05/14/2018
:Author: Akshay Paropkari
"""

import sys
import argparse
from re import findall
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
    parser = argparse.ArgumentParser(description="Plot number of heterozygous sites per "
                                     "genomes for 1000 genome project.")
    parser.add_argument("main_file", help="Input file of non-reference counts. Output "
                        "file fom count_variants.py")
    parser.add_argument("-s", "--savefile",
                        help="Save the plot as an SVG file. Provide file path and file "
                        "name with extension.")
    parser.add_argument("-g", "--gender", action="store_true",
                        help="Set this parameter if you need counts for males. Default is"
                        " False, which will not run a haploid check (required for males)")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Plot variant calls
    colors = Blues_9.hex_colors[4:] + Greens_9.hex_colors[4:] + Purples_9.hex_colors[4:] \
             + Reds_9.hex_colors[2:6] + YlOrBr_9.hex_colors[2:]
    x_order = tuple(["FIN", "GBR", "CEU", "IBS", "TSI", "CHS", "CDX", "CHB",
                     "JPT", "KHV", "GIH", "STU", "PJL", "ITU", "BEB", "PEL",
                     "MXL", "CLM", "PUR", "ASW", "ACB", "GWD", "YRI", "LWK",
                     "ESN", "MSL"])
    assigned_colors = {a: b for a, b in zip(x_order, colors)}
    try:
        count_data_df = pd.read_csv(args.main_file, sep="\t", index_col=False)
    except Exception as ee:
        sys.exit("ERROR\n{}".format(ee))
    if args.gender:
        count_data_df = count_data_df.query("gender == 'male'")
    if not args.gender:
        count_data_df = count_data_df.query("gender == 'female'")
    count_sorted_df = count_data_df.sort_values(by=["pop", "htz_counts"])
    x_order = count_sorted_df.groupby(["pop"]).mean().sort_values(by="htz_counts")
    with mpl.style.context("seaborn-white"):
        plt.figure(figsize=(12, 6))
        for entry in x_order.itertuples():
            plt.scatter(entry[0], entry[1], marker="_", s=300,
                        c="k")
        for entry in count_sorted_df.itertuples():
            plt.scatter(entry[2], entry[5], marker=".",
                        c=assigned_colors[entry[2]], edgecolors="face")
        chr_num =  findall(r"(\d+)", args.main_file.split("_")[0])[0]
        if args.gender:
            plt.title("Heterozygosity in male chromosome X samples for all populations",
                      fontsize=12)
        else:
            plt.title("Heterozygosity in chromosome {}".format(chr_num),
                      fontsize=12)
        plt.grid(axis="both", linestyle=":")
        plt.tight_layout()
        if args.savefile:
            plt.savefig(args.savefile, dpi=300, format="svg", bbox_inches="tight")
        else:
            plt.show()


if __name__ == "__main__":
    sys.exit(main())
