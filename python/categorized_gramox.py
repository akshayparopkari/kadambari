#!/usr/bin/env python
"""
Abstract: Get gram stain and oxygen requirement (gramox) data by category.

Date: 03/11/2015

Author: Akshay Paropkari
"""

import sys
import argparse
from collections import defaultdict


def categorize_otus(catdata):
    """
    Categorize OTU"s based on their gram stain and oxygen requirement.

    :type catdata: defaultdict(list)
    :param catdata: OTU data obtained from overall gramox data file.
    """
    gramox = ["|Gram_Negative_Anaerobes", "|Gram_Negative_Aerobes",
              "|Gram_Positive_Anaerobes", "|Gram_Positive_Aerobes",
              "|Unknown"]
    cat_gramox = defaultdict(list)

    Gram_Negative_Anaerobe = ["0", "0"]
    Gram_Negative_Aerobe = ["0", "1"]
    Gram_Positive_Anaerobe = ["1", "0"]
    Gram_Positive_Aerobe = ["1", "1"]

    for k1, v1 in catdata.iteritems():
        for k2, v2 in v1.iteritems():
            if v2 == Gram_Negative_Anaerobe:
                cat_gramox[k1+gramox[0]].append(k1)
            elif v2 == Gram_Negative_Aerobe:
                cat_gramox[k1+gramox[1]].append(k2)
            elif v2 == Gram_Positive_Anaerobe:
                cat_gramox[k1+gramox[2]].append(k2)
            elif v2 == Gram_Positive_Aerobe:
                cat_gramox[k1+gramox[3]].append(k2)
            elif "NA" in v2:
                cat_gramox[k1+gramox[4]].append(k2)
    return cat_gramox


def handle_program_options():
    parser = argparse.ArgumentParser(
        description="Get categorized gramox data for OTU\"s from annotated "
                    "relative abundance file.")
    parser.add_argument("rel_abd_fnh",
                        help="Path to relative abundance data file. This file"
                             " is the transposed output from "
                             "biom_relative_abundance and its second column "
                             "must be the category column.")
    parser.add_argument("in_gramox_fnh",
                        help="Path to completed gramox file. All entries "
                             "must be present, unidentified entries shold be "
                             "labeled 'NA'.")
    parser.add_argument("out_gramox_fnh",
                        help="Categorized gramox data will be written to "
                             "this file.")
    parser.add_argument("otu_dist_fnh",
                        help="Categorized OTU will be written into this file.")
    parser.add_argument("sample_otu_fnh",
                        help="Per sample OTU will be written into this file.")
    parser.add_argument("categories", nargs="+",
                        help="Provide category names. For eg. 'A' 'B' 'C'")
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.rel_abd_fnh):
            pass
    except IOError as ioe:
        err_msg = "\nError opening relative abundance file: {}\n"
        sys.exit(err_msg.format(ioe))

    try:
        with open(args.in_gramox_fnh):
            pass
    except IOError as ioe:
        err_msg = "\nError opening completed gramox file: {}\n"
        sys.exit(err_msg.format(ioe))

# Read gramox data.
    with open(args.in_gramox_fnh, "rU") as gdf:
        gramox_data = {line.strip().split("\t")[0]: line.strip().split("\t")[1:]
                       for line in gdf.readlines()[1:]}

# Read category names.
    categories = args.categories
    cat = defaultdict()
    for c in categories:
        cat[c] = {}

# Get the bacteria based on categories and sampleid.
    sample_otu = defaultdict(list)
    with open(args.rel_abd_fnh, "rU") as bacat:
        bacteria = [entry
                    for entry in bacat.readline().strip().split("\t")[2:]]
        for line in bacat:
            line = line.strip().split("\t")
            for i, abd in enumerate(line[2:]):
                if float(abd) > 0:
                    sample_otu[line[0]].append(bacteria[i])
            for c in categories:
                if line[1] == c:
                    for i, data in enumerate(line[2:]):
                        if float(data) > 0:
                            cat[c][bacteria[i]] = gramox_data[bacteria[i]]

# Obtain numbers for all gramox categories.
    final_data = categorize_otus(cat)
    ord_keys = sorted(final_data.keys())

# Write to output file.
    with open(args.out_gramox_fnh, "w") as outf:
        outf.write("Category\tGramox\tCounts\n")
        for k in ord_keys:
            category = k.split("|")[0]
            gramox = k.split("|")[1]
            outf.write("{}\t{}\t{}\n".format(" ".join(category.split("_")),
                                             " ".join(gramox.split("_")),
                                             len(final_data[k])))

# Write out OTU per sampleid
    with open(args.sample_otu_fnh, "w") as sdf:
        sdf.write("SAMPLEID\tOTU\tGRAM STAIN\tOXYGEN REQUIREMENT\n")
        for k, v in sample_otu.iteritems():
            for item in v:
                sdf.write("{}\t{}\t{}\t{}\n".format(k, item,
                          gramox_data[item][0], gramox_data[item][1]))

# Write out categorized OTU
    with open(args.otu_dist_fnh, "w") as sdf:
        sdf.write("CATEGORY\tOTU\tGRAM STAIN\tOXYGEN REQUIREMENT\n")
        for k, v in cat.iteritems():
            for k1, v1 in v.iteritems():
                sdf.write("{}\t{}\t{}\t{}\n".format(k, k1, v1[0], v1[1]))

if __name__ == "__main__":
    main()
