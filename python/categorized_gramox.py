#!/usr/bin/env python
"""
Abstract: Get gram stain and oxygen requirement (gramox) data by category and 
          sampleID.

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


def write_output(fnh1, fnh2, output, classify="category"):
    """
    Write the output to file.

    :type fnh: str
    :param fnh: Name of the output file.
    :type output: defaultdict(dict)
    :param output: Dict containing output data keyed on categories or sampleIDs
                   and their values is a dict keyed on gram stain and oxygen
                   requirement acronym with OTU's falling in that particular
                   group as the dict value.
    :type classify: str
    :param classify: Specify if the output written is classified by category or
                     by sampleIDs. The dafult option is 'category'.
    """
    if classify == "sampleID":
        total_otus = {sid: float(len(output[sid])) for sid in output.keys()}
        final_sid_data = categorize_otus(output)
        ord_sid_keys = sorted(final_sid_data.keys())
        with open(fnh1, "w") as outf:
            outf.write("SampleID\tGramox\tCounts\tTotal OTUs in SampleID\t\
                        OTU Percent in Gramox\n")
            for k in ord_sid_keys:
                sampleid = k.split("|")[0]
                gramox = k.split("|")[1]
                tot = total_otus[sampleid]
                pct = (len(final_sid_data[k]) / tot) * 100
                outf.write("{}\t{}\t{}\t{}\t{}\n".format(
                    " ".join(sampleid.split("_")), " ".join(gramox.split("_")),
                    len(final_sid_data[k]), total_otus[k.split("|")[0]], pct))
        with open(fnh2, "w") as out2f:
            out2f.write("SAMPLEID\tOTU\tGRAM STAIN\tOXYGEN REQUIREMENT\n")
            for k, v in output.iteritems():
                for k1, v1 in v.iteritems():
                    out2f.write("{}\t{}\t{}\t{}\n".format(k, k1, v1[0], v1[1]))
    else:
        total_otus = {sid: len(output[sid]) for sid in output.keys()}
        final_cat_data = categorize_otus(output)
        ord_cat_keys = sorted(final_cat_data.keys())
        with open(fnh1, "w") as outf:
            outf.write("Category\tGramox\tCounts\tTotal OTUs in SampleID\t\
                        OTU Percent in Gramox\n")
            for k in ord_cat_keys:
                category = k.split("|")[0]
                gramox = k.split("|")[1]
                tot = float(total_otus[k.split("|")[0]])
                pct = (len(final_cat_data[k]) / tot) * 100
                outf.write("{}\t{}\t{}\t{}\t{}\n".format(
                    " ".join(category.split("_")), " ".join(gramox.split("_")),
                    len(final_cat_data[k]), total_otus[k.split("|")[0]], pct))
        with open(fnh2, "w") as sdf:
            sdf.write("CATEGORY\tOTU\tGRAM STAIN\tOXYGEN REQUIREMENT\n")
            for k, v in output.iteritems():
                for k1, v1 in v.iteritems():
                    sdf.write("{}\t{}\t{}\t{}\n".format(k, k1, v1[0], v1[1]))


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
    parser.add_argument("cat_gramox_fnh",
                        help="Categorized gramox data will be written to "
                             "this file.")
    parser.add_argument("sid_gramox_fnh",
                        help="Sample-wise gramox data will be written to "
                             "this file.")
    parser.add_argument("cat_otu_fnh",
                        help="OTU's in each category will be written to this\
                             file.")
    parser.add_argument("sid_otu_fnh",
                        help="OTU's in each sample will be written to this\
                             file.")
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
        gramox_data = {line.strip().split("\t")[0]: line.strip().split("\t")
                       [1:] for line in gdf.readlines()[1:]}

# Read category names.
    categories = args.categories
    cat = defaultdict()
    for c in categories:
        cat[c] = {}

# Get the bacteria based on categories and sampleid.
    sid_otu = defaultdict(dict)
    with open(args.rel_abd_fnh, "rU") as bacat:
        bacteria = [entry
                    for entry in bacat.readline().strip().split("\t")[2:]]
        for line in bacat:
            line = line.strip().split("\t")
            for i, abd in enumerate(line[2:]):
                if float(abd) > 0:
                    sid_otu[line[0]][bacteria[i]] = gramox_data[bacteria[i]]
            for c in categories:
                if line[1] == c:
                    for i, data in enumerate(line[2:]):
                        if float(data) > 0:
                            cat[c][bacteria[i]] = gramox_data[bacteria[i]]

# Write to output file.
    write_output(args.sid_gramox_fnh, args.sid_otu_fnh, sid_otu, "sampleID")
    write_output(args.cat_gramox_fnh, args.cat_otu_fnh, cat)

if __name__ == "__main__":
    main()
