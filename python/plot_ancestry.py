#!/usr/bin/env python
"""
:Abstract: Plot number of non-reference sites per genome in ChrX of 1000 genome project.
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
    parser = argparse.ArgumentParser(description="Plot number of non-reference sites per "
                                     "genome in ChrX of 1000 genome project.")
    parser.add_argument("main_file", help="Input file of non-reference counts. Output "
                        "file fom count_variants.py")
    parser.add_argument("savefile", help="Save the plot as an SVG file. Provide file path"
                        " and file name.")
    parser.add_argument("-pd", "--plot_ancestry", nargs="+",
                        help="Supply local ancestry information file from 1000 genomes "
                        "project.")
    parser.add_argument("-g", "--gender", action="store_true",
                        help="Set this parameter if you need plot for males. Default is"
                        " False, which will not run a haploid check (required for males)")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Plot number of variants vs ancestry (NAT, AFR, EUR) per genome
    try:
        vc_data = pd.read_csv(args.main_file, sep="\t", index_col=False, header=0,
                              names=["ID", "pop", "super_pop",
                                     "gender", "variant_counts"])
    except Exception as ee:
        sys.exit("Error reading local ancestry information file\n{}".format(ee))
    else:
        all_data = pd.DataFrame()
        for infile in args.plot_ancestry:
            try:
                lai_data = pd.read_csv(infile, sep="\t", index_col=False)
            except Exception as ex:
                sys.exit("Error reading variant count data file: {}".format(ex))
            else:
                merged_data = pd.merge(lai_data, vc_data, on="ID")
                all_data = pd.concat([all_data, merged_data], ignore_index=True,
                                     verify_integrity=True)
        colors = Blues_9.hex_colors[4:] + Greens_9.hex_colors[4:] + \
                 Purples_9.hex_colors[4:] + Reds_9.hex_colors[2:6] + \
                 YlOrBr_9.hex_colors[2:]
        x_order = tuple(["FIN", "GBR", "CEU", "IBS", "TSI", "CHS", "CDX", "CHB",
                         "JPT", "KHV", "GIH", "STU", "PJL", "ITU", "BEB", "PEL",
                         "MXL", "CLM", "PUR", "ASW", "ACB", "GWD", "YRI", "LWK",
                         "ESN", "MSL"])
        assigned_colors = {a: b for a, b in zip(x_order, colors)}
        fig, axarr = plt.subplots(3, figsize=(10, 8))
        for entry in all_data.itertuples():
            try:
                axarr[0].scatter(entry[4], entry[9], marker=".",
                                 c=assigned_colors[entry[6]])
                axarr[1].scatter(entry[2], entry[9], marker=".",
                                 c=assigned_colors[entry[6]])
                axarr[2].scatter(entry[3], entry[9], marker=".",
                                 c=assigned_colors[entry[6]])
            except Exception:
                continue
        axarr[0].set_title("Native American Ancestry", fontdict={"fontsize": 18})
        axarr[1].set_title("African Ancestry", fontdict={"fontsize": 18})
        axarr[2].set_title("European Ancestry", fontdict={"fontsize": 18})
        legend_dict = {"Europe": Blues_9.hex_colors[5],
                       "East Asia": Greens_9.hex_colors[4:][3],
                       "South Asia": Purples_9.hex_colors[4:][3],
                       "Sub-Saharan Africa": YlOrBr_9.hex_colors[3],
                       "Americas": Reds_9.hex_colors[7]}
        l = [plt.scatter([], [], c=legend_dict[lab], s=100)
             for lab in legend_dict.keys()]
        for ax in axarr:
            ax.set_ylabel("Frequecy of heterozygosity per genome in chromosome {}".
                          format(findall(r"(\d+)", args.main_file.split("_")[0])[0]),
                          fontdict={"fontsize": 12})
            ax.grid(axis="x", linestyle=":")
            ax.legend(l, ["{}".format(lab) for lab in legend_dict.keys()],
                      fontsize=14, scatterpoints=3, loc="upper left",
                      bbox_to_anchor=(1, 1))
        plt.tight_layout()
        try:
            plt.savefig("{}.svg".format(args.savefile), dpi=300, format="svg",
                        bbox_inches="tight")
        except Exception as ee:
            sys.exit("Error while saving the image file\n{}".format(ee))


if __name__ == "__main__":
    sys.exit(main())
