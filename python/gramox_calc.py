#!/usr/bin/env python
"""
Abstract: Calculate gram stain and oxygen requirement (gramox) for all OTU from BIOM
          format file.
Author: Akshay Paropkari
Date: 06/29/2015
"""
import sys
import argparse
importerrors = []
try:
    import biom
except ImportError:
    importerrors.append("biom")
try:
    from phylotoast.otu_calc import otu_name
except ImportError:
    importerrors.append("phylotoast")
if len(importerrors) > 0:
    for err in importerrors:
        print("Please install missing package: {}".format(err))
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description="Get gram stain and oxygen requirement "
                                     "(gramox) data for s-OTUs from a given BIOM file. "
                                     "Gram Positives are 1, Gram Negatives are 0. "
                                     "Obligate Anaerobes are 0 and rest are classified as"
                                     " Aerobes with 1.")
    parser.add_argument("master_fnh", help="Path to master gramox data file [REQUIRED]")
    parser.add_argument("biom_file", help="BIOM file/OTU table. This file is used to "
                        "collect the list of OTUs for processing their gramox data "
                        "[REQUIRED]")
    parser.add_argument("out_fnh", help="Path to output file [REQUIRED]")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Read gramox data from master file
    try:
        with open(args.master_fnh, "rU") as bgof:
            bgodata = {line.strip().split("\t")[1]: "\t".join(line.strip().split("\t")[2:4])
                       for line in bgof.readlines()}
    except Exception as err:
        sys.exit("\nError parsing master gramox file: {}. Please check master gramox data"
                 " file and re-run this script.".format(err))

    # Read relative abundance data
    try:
        biomf = biom.load_table(args.biom_file)
    except Exception as err:
        sys.exit("\nError opening BIOM file: {}\n".format(err))
    else:
        otus = [otu_name(biomf.metadata(otuid, "observation")["taxonomy"])
                for otuid in biomf.ids("observation")]

    # Write classified gramox data to tsv file
    print("\nWriting out results to file. For OTUs with missing gramox data, please "
          "manually input the relevant information or fill in NA for unverified "
          "information. This is required before running pct_abd_gramox.py script.\n")
    with open(args.out_fnh, "w") as gramoxout:
        gramoxout.write("#OTU\tGram Status\tOxygen Requirement\n")
        for otu in otus:
            if otu in bgodata.keys():
                gramoxout.write("{}\t{}\n".format(otu, bgodata[otu]))
            else:
                gramoxout.write("{}\n".format(otu))


if __name__ == "__main__":
    sys.exit(main())
