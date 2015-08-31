#!/usr/bin/env python
"""
Abstract: Fill up the empty cells in gramox data with 'NA' due to
          unavailability of the information. Other cells remain unchanged.

Date: 02/27/2015

Author: Akshay Paropkari
"""
import argparse
import sys


def handle_program_options():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Fill up the empty cells in gramox data with "NA" due to'
                    '\nunavailability of the information. Other cells remain'
                    '\nunchanged. Gramox data format:'
                    '\n#OTU\tGram Status\tOxygen Requirement\tSource'
        )
    parser.add_argument('in_fnh',
                        help='Path to input gramox data file.')
    parser.add_argument('out_fnh',
                        help='Path to output gramox file.')
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.in_fnh):
            pass
    except IOError as ioe:
        err_msg = '\nError opening input data file: {}\n'
        sys.exit(err_msg.format(ioe))

    with open(args.in_fnh, 'rU') as f:
        gramoxdata = [line.strip().split('\t') for line in f]

# Filling up the missing values with 'NA'.
    for entry in gramoxdata:
        for i, item in enumerate(entry):
            if item == '':
                entry[i] = 'NA'    # For missing gram status, oxygen req data.
        if len(entry) == 1:
            length = range(3 - len(entry))
            for l in reversed(length):
                l += 1
                entry.insert(l, 'NA')    # For unidentified OTUs with no info.

    with open(args.out_fnh, 'w') as ogf:
        for entry in gramoxdata:
            entry = '\t'.join(entry)
            ogf.write('{}\n'.format(entry))

if __name__ == '__main__':
    main()
