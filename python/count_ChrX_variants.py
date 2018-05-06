#!/usr/bin/env python
"""
:Abstract: Calculate number of variant sites per genome in ChrX of 1000 genome project.
:Date: 05/03/2018
:Author: Akshay Paropkari
"""

import sys
import gzip
import argparse
from collections import defaultdict
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
    parser = argparse.ArgumentParser(description="Calculate number of variant sites per "
                                     "genome in ChrX of 1000 genome project.")
    parser.add_argument("-vcf", "--chr21_vcf_file", help="Path to input ChrX VCF file")
    parser.add_argument("-md", "--map_fp",
                        help="Metadata mapping file corresponsding to ChrX VCF")
    parser.add_argument("-mf", "--main_file", help="Input file of variant counts.")
    parser.add_argument("-s", "--savefile",
                        help="Save the plot as an SVG file. Provide file path and file "
                        "name with extension.")
    parser.add_argument("-o", "--output_file",
                        help="Save consolidated data a tab-separated file. Provide file "
                        "path and file name with extension.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Get metadata for each genome
    if args.map_fp:
        with open(args.map_fp, "r") as md:
            md_data = dict()
            for line in md.readlines()[1:]:
                line = line.split()
                md_data[line[0]] = line[1]

    variants = ["1|0", "0|1"]

    # Get variant site data for all genomes
    if args.chr21_vcf_file:
        with gzip.open(args.chr21_vcf_file, "rb") as vcff:
            for line in vcff:
                try:
                    line = line.decode()
                    assert line.startswith("#CHROM")
                except Exception:
                    try:
                        assert line.startswith("X")
                        line = line.strip().split("\t")
                        assert line[3] in ["A", "T", "C", "G"]
                        assert line[4] in ["A", "T", "C", "G"]
                    except Exception:
                        continue
                    else:
                        for i, entry in enumerate(line[9:]):
                            try:
                                assert entry in variants
                            except Exception:
                                continue
                            else:
                                variant_data[genome_order[i]] += 1
                else:
                    line = line.strip().split("\t")
                    genome_order = line[9:]
                    variant_data = {col: 0 for col in genome_order}

        # Consolidate data
        all_data = defaultdict(list)
        for sample in md_data.keys():
            all_data[sample] = [md_data[sample], variant_data[sample]]
        all_data_df = pd.DataFrame.from_dict(all_data, orient="index")
        all_data_df.columns = ["population", "variant counts"]
        if args.output_file:
            all_data_df.to_csv(args.output_file, sep="\t")

    # The code above is run on MERCED, and below is using output of variant counts and
    # plotting the data.

    # Plot variant calls
    if args.main_file:
        colors = Blues_9.hex_colors[4:] + Greens_9.hex_colors[4:] +\
                 Purples_9.hex_colors[4:] + Reds_9.hex_colors[2:6] +\
                 YlOrBr_9.hex_colors[2:]
        x_order = tuple(["FIN", "GBR", "CEU", "IBS", "TSI", "CHS", "CDX", "CHB",
                         "JPT", "KHV", "GIH", "STU", "PJL", "ITU", "BEB", "PEL",
                         "MXL", "CLM", "PUR", "ASW", "ACB", "GWD", "YRI", "LWK",
                         "ESN", "MSL"])
        assigned_colors = {a: b for a, b in zip(x_order, colors)}
        count_data_df = pd.read_csv(args.main_file, sep="\t", index_col=0)
        admixed_pop = ["CEU", "GIH", "PEL", "MXL", "CLM", "PUR", "ASW", "ACB"]
        count_subset_df = count_data_df.query("population in {}".format(admixed_pop))
        count_subset_df = count_subset_df.sort_values(by=["population", "variant counts"])
        x_order = count_subset_df.groupby(["population"]).mean().sort_values(by="variant counts")
        with mpl.style.context("seaborn-white"):
            plt.figure(figsize=(7, 5))
            for entry in x_order.itertuples():
                plt.scatter(entry[0], entry[1], marker="_", s=300,
                            c="k")
            for entry in count_subset_df.itertuples():
                plt.scatter(entry[1], entry[2], marker=".",
                            c=assigned_colors[entry[1]], edgecolors="face")
            plt.xlabel("Admixed Populations", fontsize=12)
            plt.ylabel("Number of Variant Sites", fontsize=12)
            plt.tight_layout()
            plt.grid(axis="y", linestyle=":")
            if args.savefile:
                plt.savefig(args.savefile, dpi=300, format="svg",
                            bbox_inches="tight", pad_inches=0.25)
            else:
                plt.show()


if __name__ == "__main__":
    sys.exit(main())
