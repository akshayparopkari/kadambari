#!/usr/env/bin python

"""
Abstract: Given a file of gene or OTU IDs per group, this script will output
          groups present in each gene or OTU IDs.
"""

import sys
import csv
import argparse
from collections import defaultdict


def handle_program_options():
    parser = argparse.ArgumentParser(
        description=("Given a file of gene or OTU IDs per group, this script "
                     "will output groups present in each gene or OTU IDs.")
    )
    parser.add_argument("in_fnh",
                        help="Input file containing identifiers as columns "
                             "and gene or OTU IDs as rows, one per line.")
    parser.add_argument("out_fnh",
                        help="Output file handle name. First column contains "
                             "gene or OTU IDs and subsequenct columns contain "
                             "identifiers/groups they are present in.")
    return parser.parse_args()


def main():
    args = handle_program_options()
    # Get input
    with open(args.in_fnh, "rU") as tykyuo:
        reader = csv.DictReader(tykyuo, delimiter="\t")
        data = [line for line in reader]

    # Collect groups from input data - column names are assumed as group names
    groups = data[1].keys()

    # Get master list of gene or OTU IDs
    all_ids = set()
    for d in data:
        for grp in groups:
            if d[grp] != "":
                all_ids.add(d[grp])

    # Compare master list to all groups
    id_groups = defaultdict(list)
    for d in data:
        for i in sorted(all_ids):
            for grp in groups:
                if i == d[grp]:
                    id_groups[i].append(grp)

    # Write to output
    with open(args.out_fnh, "w") as uyiouy:
        for k, v in id_groups.iteritems():
            uyiouy.write("{}\t{}\n".format(k, "\t".join(sorted(v))))

if __name__ == "__main__":
    sys.exit(main())
