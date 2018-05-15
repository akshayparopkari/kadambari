#!/usr/bin/env python
"""
:Abstract: Calculate number of variant sites per genome in ChrX of 1000 genome project.
:Date: 05/03/2018
:Author: Akshay Paropkari
"""

import sys
import gzip
import argparse
from time import strftime
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


def get_variant_counts(line, males, linenum):
    """Iterate through dataframe of genotype entries and count 0|1 or 1|0."""
    print("{}: Processing line {}".format(strftime("%d %b %Y %H:%M:%S"), linenum))
    try:
        line = line.strip().split("\t")
        if males:
            try:
                assert 60001 <= int(line[1]) <= 2699520
            except AssertionError:
                try:
                    154931044 <= int(line[1]) <= 155260560
                except AssertionError:
                    pass
                else:
                    print("{}: Skipped line because coordinate in ignored region".
                          format(strftime("%d %b %Y %H:%M:%S"), linenum))
                    return None
            else:
                print("{}: Skipped line because coordinate in ignored region".
                          format(strftime("%d %b %Y %H:%M:%S"), linenum))
                return None
        assert line[3] in ["A", "T", "C", "G"]
        assert line[4] in ["A", "T", "C", "G"]
    except AssertionError:
        print("{}: Skipped line since position isn't biallelic {}".
              format(strftime("%d %b %Y %H:%M:%S"), linenum))
        return None
    else:
        variant_indices = [i for i, entry in enumerate(line[9:])
                           if entry in ["1|0", "0|1"]]
        print("{}: Completed line {}".format(strftime("%d %b %Y %H:%M:%S"), linenum))
        return variant_indices


def handle_program_options():
    parser = argparse.ArgumentParser(description="Calculate number of variant sites per "
                                     "genome in ChrX of 1000 genome project.")
    parser.add_argument("-vcf", "--vcf_file", help="Path to input ChrX VCF file")
    parser.add_argument("-md", "--map_fp",
                        help="Metadata mapping file corresponsding to ChrX VCF")
    parser.add_argument("-g", "--gender", action="store_true",
                        help="Set this parameter if you need counts for males. Default is"
                        " False, which will not run a haploid check (required for males)")
    parser.add_argument("-o", "--output_file",
                        help="Save consolidated data a tab-separated file. Provide file "
                        "path and file name with extension.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    #  Obtain non-reference site counts for all individuals
    if args.vcf_file:
        with gzip.open(args.vcf_file, "rb") as vcff:
            for line in vcff:
                for i, line in enumerate(vcff):
                    try:
                        line = line.decode()
                        assert line.startswith("#CHROM")
                    except Exception:
                        try:
                            assert line.startswith("X")
                        except AssertionError:
                            continue
                        else:
                            input = [line, args.gender, i]
                            res = get_variant_counts(*input)
                            try:
                                assert res is not None
                            except AssertionError:
                                continue
                            else:
                                for entry in res:
                                    variant_data[genome_order[entry]] += 1
                    else:
                        line = line.strip().split("\t")
                        genome_order = line[9:]
                        variant_data = {col: 0 for col in genome_order}

        # Consolidate data
        # Get metadata for each genome
        if args.map_fp:
            md_data = pd.read_csv(args.map_fp, sep="\t", index_col=False,
                                  usecols=[0, 1, 2, 3])
            md_data["variant_counts"] = md_data["sample"].map(variant_data)
            if args.output_file:
                md_data.to_csv(args.output_file, sep="\t", index=False)
        else:
            sys.exit("Please supply --map_fp parameter with the metadata file.")


if __name__ == "__main__":
    sys.exit(main())
