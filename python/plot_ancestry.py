#!/usr/bin/env python
"""
:Abstract: Plot number of non-reference sites vs ancestry proportion per genome in autosomes of 1000 genome project.
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
    from palettable.cartocolors.qualitative import Safe_6    # 6 populations
except ImportError:
    err.append("palettable")
try:
    import pandas as pd
except ImportError:
    err.append("pandas")
try:
    from numpy import polyfit, poly1d
except ImportError:
    err.append("numpy")
try:
    assert len(err) == 0
except AssertionError:
    for error in err:
        sys.exit("Please install {}".format(error))


def handle_program_options():
    parser = argparse.ArgumentParser(description="Plot number of non-reference sites per "
                                     "genome vs ancestry proportion in autosomes of 1000 "
                                     "genome project.")
    parser.add_argument("main_file", nargs="+", help="Input file of non-reference counts."
                        " Output file fom count_variants.py")
    parser.add_argument("-s", "--savefile", action="store_true",
                        help="Save the plot as an SVG file. Provide file path and file "
                        "name.")
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
    for vcfile in args.main_file:
        try:
            vc_data = pd.read_csv(vcfile, sep="\t", index_col=False, header=0,
                                  names=["ID", "pop", "super_pop",
                                         "gender", "variant_counts"])
        except Exception as ee:
            sys.exit("Error reading variant count data file\n{}".format(ee))
        else:
            all_data = pd.DataFrame()
            for infile in args.plot_ancestry:
                try:
                    lai_data = pd.read_csv(infile, sep="\t", index_col=False)
                except Exception as ex:
                    sys.exit("Error reading local ancestry information file: {}".
                             format(ex))
                else:
                    merged_data = pd.merge(lai_data, vc_data, on="ID")
                    all_data = pd.concat([all_data, merged_data], ignore_index=True,
                                         verify_integrity=True)
            colors = Safe_6.hex_colors
            pop_labels = tuple(["PEL", "MXL", "CLM", "PUR", "ASW", "ACB"])
            assigned_colors = {a:b for a, b in zip(pop_labels, colors)}
            male_data = all_data.query("gender == 'male'")       # male only data
            female_data = all_data.query("gender == 'female'")   # female only data
            for subset_data in [male_data, female_data]:
                fig, axarr = plt.subplots(3, figsize=(10, 15))
                for entry in subset_data.itertuples():
                    try:
                        axarr[0].scatter(entry[4], entry[9], marker=".",
                                         c=assigned_colors[entry[6]])
                        axarr[1].scatter(entry[2], entry[9], marker=".",
                                         c=assigned_colors[entry[6]])
                        axarr[2].scatter(entry[3], entry[9], marker=".",
                                         c=assigned_colors[entry[6]])
                    except Exception:
                        continue
                for anc, ax in zip(["NAT", "AFR", "EUR"], axarr):
                    x = subset_data[anc].values
                    y = subset_data["variant_counts"]
                    z = polyfit(x, y, 3)
                    p = poly1d(z)
                    ax.plot(sorted(x), p(sorted(x)), "k-")
                axarr[0].set_title("Native American Ancestry", fontdict={"fontsize": 18})
                axarr[1].set_title("African Ancestry", fontdict={"fontsize": 18})
                axarr[2].set_title("European Ancestry", fontdict={"fontsize": 18})
                for ax in axarr:
                    ax.set_ylabel("Frequency of heterozygosity in chromosome {}".
                                  format(findall(r"(\d+)", vcfile.split("_")[1])[0]),
                                  fontdict={"fontsize": 12})
                    ax.grid(axis="x", linestyle=":")
                    l = [plt.scatter([], [], c=assigned_colors[lab], s=100)
                         for lab in pop_labels]
                    ax.legend(l, ["{}".format(lab) for lab in pop_labels], fontsize=18,
                              bbox_to_anchor=(1,1), loc="upper left", scatterpoints=3)
                plt.tight_layout()
                try:
                    if subset_data["gender"].values[0] == "male":
                        plt.savefig("ancestry_figs/{}_male.svg".
                                    format(vcfile.split("_")[1].split("/")[1]),
                                    dpi=300, format="svg", bbox_inches="tight")
                    else:
                        plt.savefig("ancestry_figs/{}_female.svg".
                                    format(vcfile.split("_")[1].split("/")[1]),
                                    dpi=300, format="svg", bbox_inches="tight")
                except Exception as ee:
                    sys.exit("Error while saving the image file\n{}".format(ee))


if __name__ == "__main__":
    sys.exit(main())
