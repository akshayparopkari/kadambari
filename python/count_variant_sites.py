#!/usr/bin/env python
"""
:Abstract: Calculate number of variant sites per genome in Chr21 of 1000 genome project.
:Date: 02/12/2018
:Author: Akshay Paropkari
"""

import sys
import gzip
import argparse
import itertools
from collections import defaultdict
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc("font", family="Arial")
mpl.rc("xtick", labelsize=9.5)  # set X axis ticksize
mpl.rc("ytick", labelsize=11)  # set Y axis ticksize


def handle_program_options():
    parser = argparse.ArgumentParser(description="Calculate number of variant sites per "
                                     "genome in Chr21 of 1000 genome project.")
    parser.add_argument("-vcf", "--chr21_vcf_file", help="Path to input Chr21 VCF file")
    parser.add_argument("-md", "--map_fp",
                        help="Metadata mapping file corresponsding to Chr21 VCF")
    parser.add_argument("-mf", "--main_file", help="Output file of variant counts.")
    parser.add_argument("-s", "--savefile",
                        help="Save the plot to this file. PDF preferred.")
    parser.add_argument("-o", "--output_file", action="store_true",
                        help="Save consolidated data to this file")
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

    variants = ["1|0", "0|1", "1|1"]

    # Get variant site data for all genomes
    if args.chr21_vcf_file:
        with gzip.open(args.chr21_vcf_file, "rb") as vcff:
            for line in vcff.readlines():
                try:
                    line = line.decode()
                    assert line.startswith("#CHROM")
                except Exception:
                    try:
                        assert line.startswith("21")
                    except Exception:
                        continue
                    else:
                        line = line.strip().split()
                        for i, entry in enumerate(line[9:]):
                            try:
                                assert entry in variants
                            except Exception:
                                continue
                            else:
                                variant_data[genome_order[i]] += 1
                else:
                    line = line.strip().split()
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
        count_data_df = pd.read_csv(args.main_file, sep="\t", index_col=0)
        count_data_df = count_data_df.sort_values(by=["population", "variant counts"])
        x_order = count_data_df.groupby(["population"]).mean().sort_values(by="variant counts")
        fig = plt.figure(figsize=(12, 8))
        for entry in x_order.itertuples():
            plt.scatter(entry[0], entry[1], marker="_", s=300, c="r")
        for entry in count_data_df.itertuples():
            plt.scatter(entry[1], entry[2], marker="+", c="b")
        plt.xlabel("Individual", fontsize=15)
        plt.ylabel("Variant sites per genome", fontsize=15)
        if args.savefile:
            plt.savefig(args.savefile, dpi=300, format="svg",
                        bbox_inches="tight", pad_inches=0.25)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    sys.exit(main())
