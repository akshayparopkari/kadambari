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
from re import findall
try:
    import pandas as pd
except ImportError:
    sys.exit("Please install pandas")


def get_variant_counts(line, males, linenum):
    """Iterate through dataframe of genotype entries and count 0|1 or 1|0."""
    print("{}: Processing line {}".format(strftime("%d %b %Y %H:%M:%S"), linenum))
    try:
        line = line.strip().split("\t")
        if males:
            try:
                assert 60001 <= int(line[1]) <= 2699520
            except AssertionError:  # fail assertion
                try:
                    154931044 <= int(line[1]) <= 155260560
                except AssertionError:  # fail assertion
                    try:
                        assert line[3] in ["A", "T", "C", "G"]
                        assert line[4] in ["A", "T", "C", "G"]
                    except AssertionError:
                        print("{}: Skipped line because coordinate in ignored region".
                              format(strftime("%d %b %Y %H:%M:%S"), linenum))
                        return None
                    else:
                        variant_indices = [i for i, entry in enumerate(line[9:])
                                           if int(entry) == 1]
                        print("{}: Completed line {}".
                              format(strftime("%d %b %Y %H:%M:%S"), linenum))
                        return variant_indices
                else:
                    print("{}: Skipped line because {} in ignored region".
                          format(strftime("%d %b %Y %H:%M:%S"), line[1], linenum))
                    return None
            else:
                print("{}: Skipped line because coordinate in ignored region".
                      format(strftime("%d %b %Y %H:%M:%S"), linenum))
                return None
        assert line[3] in ["A", "T", "C", "G"]
        assert line[4] in ["A", "T", "C", "G"]
    except AssertionError:
        print("{}: Skipped line since position isn't biallelic entry {}".
              format(strftime("%d %b %Y %H:%M:%S"), linenum))
        return None
    else:
        variant_indices = [i for i, entry in enumerate(line[9:])
                           if entry in ["1|0", "0|1"]]
        print("{}: Completed line {}".format(strftime("%d %b %Y %H:%M:%S"), linenum))
        return variant_indices


def get_diploid_htz_counts(line, linenum):
    """Iterate through dataframe of genotype entries and count 0|1 or 1|0."""
    line = line.strip().split("\t")
    print("{}: Processing position {} in line {}".
          format(strftime("%d %b %Y %H:%M:%S"), line[1], linenum))
    try:
        assert 60001 <= int(line[1]) <= 2699520
    except AssertionError:
        try:
            assert 154931044 <= int(line[1]) <= 155260560
        except AssertionError:
            print("{}: Skipped line {} as coordinate {} is not in range of interest".
                  format(strftime("%d %b %Y %H:%M:%S"), linenum, line[1]))
            return None
        else:
            try:
                assert line[3] in ["A", "T", "C", "G"]
                assert line[4] in ["A", "T", "C", "G"]
            except AssertionError as ae:
                print("{}: Skipped line {} as it is not biallelic entry".
                      format(strftime("%d %b %Y %H:%M:%S"), linenum))
                return None
            else:  # entry is biallelic and not in ignored region
                variant_indices = [i for i, entry in enumerate(line[9:])
                                   if entry in ["0|1", "1|0"]]
                print("{}: Processed line {}".
                      format(strftime("%d %b %Y %H:%M:%S"), linenum))
                return variant_indices
    else:
        try:
            assert line[3] in ["A", "T", "C", "G"]
            assert line[4] in ["A", "T", "C", "G"]
        except AssertionError as ae:
            print("{}: Skipped line {} as it is not biallelic entry".
                  format(strftime("%d %b %Y %H:%M:%S"), linenum))
            return None
        else:  # entry is biallelic and not in ignored region
            variant_indices = [i for i, entry in enumerate(line[9:])
                               if entry in ["0|1", "1|0"]]
            print("{}: Processed line {}".
                  format(strftime("%d %b %Y %H:%M:%S"), linenum))
            return variant_indices


def get_haploid_htz_counts(line, linenum):
    """Iterate through dataframe of genotype entries and count 0|1 or 1|0."""
    line = line.strip().split("\t")
    print("{}: Processing position {} in line {}".
          format(strftime("%d %b %Y %H:%M:%S"), line[1], linenum))
    try:
        assert 60001 <= int(line[1]) <= 2699520
    except AssertionError:
        try:
            assert 154931044 <= int(line[1]) <= 155260560
        except AssertionError:    # when coordinates not in the ignore range
            try:
                assert line[3] in ["A", "T", "C", "G"]
                assert line[4] in ["A", "T", "C", "G"]
            except AssertionError:
                print("{}: Skipped line {} as it is not biallelic entry".
                  format(strftime("%d %b %Y %H:%M:%S"), linenum))
                return None
            else:  # entry is biallelic and not in ignored region
                variant_indices = [i for i, entry in enumerate(line[9:])
                                   if entry in ["0|1", "1|0"]]
                print("{}: Processed line {}".
                      format(strftime("%d %b %Y %H:%M:%S"), linenum))
                return variant_indices
        else:
            print("{}: Skipped line {} as coordinate {} is not in range of interest".
                  format(strftime("%d %b %Y %H:%M:%S"), linenum, line[1]))
            return None
    else:
        print("{}: Skipped line {} as coordinate {} is not in range of interest".
              format(strftime("%d %b %Y %H:%M:%S"), linenum, line[1]))
        return None


def get_autosome_htz_counts(line, linenum):
    """Iterate through dataframe of genotype entries and count 0|1 or 1|0."""
    line = line.strip().split("\t")
    print("{}: Processing position {} in line {}".
          format(strftime("%d %b %Y %H:%M:%S"), line[1], linenum))
    try:
        assert line[3] in ["A", "T", "C", "G"]
        assert line[4] in ["A", "T", "C", "G"]
    except AssertionError as ae:
        print("{}: Skipped line {} as it is not biallelic entry".
              format(strftime("%d %b %Y %H:%M:%S"), linenum))
        return None
    else:  # entry is biallelic and not in ignored region
        variant_indices = [i for i, entry in enumerate(line[9:])
                           if entry in ["0|1", "1|0"]]
        print("{}: Processed line {}".
              format(strftime("%d %b %Y %H:%M:%S"), linenum))
        return variant_indices


def handle_program_options():
    parser = argparse.ArgumentParser(description="Calculate number of variant sites per "
                                     "genome in ChrX of 1000 genome project.")
    parser.add_argument("-vcf", "--vcf_file", help="Path to input ChrX VCF file")
    parser.add_argument("-md", "--map_fp",
                        help="Metadata mapping file corresponsding to ChrX VCF")
    parser.add_argument("-i", "--include", action="store_true",
                        help="Supply this parameter to include the 'diploid' coordinates "
                        "of ChrX in calculations. By default, these coordinates will be "
                        "excluded.")
    parser.add_argument("-a", "--autosomes", action="store_true",
                        help="Supply this parameter to get heterozygosity counts for all"
                        " autosomes.")
    parser.add_argument("-o", "--output_file",
                        help="Save consolidated data a tab-separated file. Provide file "
                        "path and file name with extension.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    #  Obtain non-reference site counts for all individuals
    if args.vcf_file:
        with gzip.open(args.vcf_file, "rb") as vcff:
            chr_name = args.vcf_file.split(".")[1]
            try:
                chr_num = findall(r"(\d+)", chr_name)[0]
            except IndexError:
                chr_num = "X"    # when processing ChrX vcf
            for line in vcff:
                for i, line in enumerate(vcff):
                    try:
                        line = line.decode()
                        assert line.startswith("#CHROM")
                    except Exception:
                        try:
                            assert line.startswith("{}".format(chr_num))
                        except AssertionError:
                            continue
                        else:
                            if args.include:
                                res = get_diploid_htz_counts(line, i)
                            elif args.autosomes:
                                res = get_autosome_htz_counts(line, i)
                            else:
                                res = get_haploid_htz_counts(line, i)
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
            md_data["htz_counts"] = md_data["sample"].map(variant_data)
            if args.output_file:
                md_data.to_csv(args.output_file, sep="\t", index=False)
        else:
            sys.exit("Please supply --map_fp parameter with the metadata file.")


if __name__ == "__main__":
    sys.exit(main())
