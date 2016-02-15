#!/usr/bin/env python
"""
Abstract: Filter biom file on both 'sample' and 'observation' axes, given a
          list of sampieIDs to retain.
Author: Akshay Paropkari
Date: 02/15/2016
"""
import sys
import argparse
importerrors = []
try:
    import biom
    from biom.util import biom_open as bo
except ImportError:
    importerrors.append("biom-format")
try:
    import pandas as pd
except ImportError:
    importerrors.append("pandas")
if len(importerrors) > 0:
    for err in importerrors:
        print "Please install missing module: {}".format(err)
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description="Filter biom file on both \
                                                 'sample' and 'observation' \
                                                 axes, given a list of \
                                                 sampieIDs to retain.")
    parser.add_argument("input_biom_fnh", help="BIOM file path.")
    parser.add_argument("output_biom_fnh", default="filtered_otu_table.biom",
                        help="Filtered biom output file.")
    parser.add_argument("mapping_fnh",
                        help="Mapping file with sampleIDs to retain in it. The"
                             " '#SampleID' column will be used to get the list"
                             " of all ids to retain.")
    parser.add_argument("filter_otuids_fnh", default="filter_otuids.txt",
                        help="OTUIDs to filter from .tre file using QIIME "
                             "'filter_tree.py' script. Input file to be used "
                             "with '-n' option.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Error check input file
    try:
        with open(args.input_biom_fnh):
            pass
    except IOError as ioe:
        sys.exit("\nError in BIOM file path: {}\n".format(ioe))
    try:
        with open(args.mapping_fnh):
            pass
    except IOError as ioe:
        sys.exit("\nError in mapping file path: {}\n".format(ioe))

    # Read input files
    biomf = biom.load_table(args.input_biom_fnh)
    mapf = pd.read_csv(args.mapping_fnh, sep="\t")

    # Get filtered biom data
    keep_sampleIDs = list(mapf["#SampleID"])
    sid_filtered_biomf = biomf.filter(keep_sampleIDs, inplace=False)
    print "{} sampleIDs retained from original biom file.".format(len(keep_sampleIDs))
    obs_abd_sum = sid_filtered_biomf.sum(axis="observation")
    otuids = sid_filtered_biomf.ids("observation")
    abd_sum = {a: b for a, b in zip(otuids, obs_abd_sum)}

    # Run checks for redundant otus
    redundant_otuids = [otu for otu, abd in abd_sum.iteritems() if abd == 0]
    otuid_filtered_biomf = sid_filtered_biomf.filter(redundant_otuids,
                                                     "observation",
                                                     invert=True,
                                                     inplace=False)
    print "{} otuIDs filtered out of the original biom file.\n".format(len(redundant_otuids))

    # Write out files
    with bo(args.output_biom_fnh, "w") as rth:
        otuid_filtered_biomf.to_hdf5(rth, "Filtered OTU Table.")
    with open(args.filter_otuids_fnh, "w") as yui:
        for otuid in redundant_otuids:
            yui.write("{}\n".format(otuid))

if __name__ == "__main__":
    main()
