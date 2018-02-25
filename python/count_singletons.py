#!/usr/bin/env python
"""
:Abstract: Calculate number of singletons in Chr21 of 1000 genome project.
:Date: 02/21/2018
:Author: Akshay Paropkari
"""

import sys
import gzip
import argparse
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
    parser = argparse.ArgumentParser(description="Calculate number of singletons per "
                                     "population in Chr21 of 1000 genome project.")
    parser.add_argument("-vcf", "--chr21_vcf_file", help="Path to input Chr21 VCF file")
    parser.add_argument("-md", "--map_fp",
                        help="Metadata mapping file corresponsding to Chr21 VCF")
    parser.add_argument("-mf", "--main_file", help="File of singleton counts.")
    parser.add_argument("-s", "--savefile",
                        help="Save the plot to this file. PDF preferred.")
    parser.add_argument("-o", "--output_file",
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

    singletons = ["1|0", "0|1"]

    # Get singleton data for all genomes
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
                        for singleton in singletons:
                            try:
                                assert line[3] in ["A", "T", "C", "G"]
                                assert line[4] in ["A", "T", "C", "G"]
                                assert line[9:].count(singleton) == 1
                            except AssertionError:
                                continue
                            else:
                                individual = line[9:].index(singleton)
                                singleton_data[genome_order[individual]] += 1
                else:
                    line = line.strip().split()
                    genome_order = line[9:]
                    singleton_data = {col: 0 for col in genome_order}

    # Get normalized singleton data
    if args.output_file:
        with open(args.output_file, "w") as outf:
            outf.write("sampleID\tpopulation\tsingleton_count\n")
            for sid in md_data.keys():
                outf.write("{0}\t{1}\t{2}\n".format(sid, md_data[sid],
                                                    singleton_data[sid]))

    # Plot the data
    if args.savefile:
        try:
            singleton_counts = pd.read_csv(args.main_file, sep="\t", index_col=0)
        except Exception:
            sys.exit("Error reading in singleton count file. Please check the --main_file"
                     " parameter.")
        else:
            colors = Blues_9.hex_colors[4:] + Greens_9.hex_colors[4:] +\
                     Purples_9.hex_colors[4:] + Reds_9.hex_colors[5:] +\
                     YlOrBr_9.hex_colors[2:]
            x_order = tuple(["FIN", "GBR", "CEU", "IBS", "TSI", "CHS", "CDX", "CHB",
                             "JPT", "KHV", "GIH", "STU", "PJL", "ITU", "BEB", "PEL",
                             "MXL", "CLM", "PUR", "ASW", "ACB", "GWD", "YRI", "LWK",
                             "ESN", "MSL"])
            assigned_colors = {a:b for a, b in zip(x_order, colors)}
            sorted_counts = singleton_counts.groupby("population").median().\
                            reindex(x_order, axis=0)
            plt.figure(figsize=(20, 8))
            for entry in sorted_counts.itertuples():
                plt.scatter(entry[0], entry[1], marker="_", s=300, c="k")
            for entry in singleton_counts.itertuples():
                plt.scatter(entry[1], entry[2], marker="+", c=assigned_colors[entry[1]])
            plt.xlabel("Populations", fontdict={"fontsize": 12})
            plt.ylabel("Number of Singletons in Chr21", fontdict={"fontsize": 12})
            plt.savefig(args.savefile, dpi=300, format="svg",
                       bbox_inches="tight", pad_inches=0.25)


if __name__ == "__main__":
    sys.exit(main())
