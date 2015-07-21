#!usr/bin/env/python
"""
Abstract: Calculate gram stain and oxygen requirement (gramox) for all OTU
          from their relative abundance file and write the output to file.

Date: 06/29/2015

Author: Akshay Paropkari
"""
import argparse
import sys
import pandas as pd


def handle_program_options():
    parser = argparse.ArgumentParser(
        description='Get gramox data for OTU\'s from relative abundance file.'
    )
    parser.add_argument('master_fnh',
                        help='Path to master gramox data file.')
    parser.add_argument('rel_abd_fnh',
                        help='Path to relative abundance data file.')
    parser.add_argument('out_fnh',
                        help='Path to output file.')
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.master_fnh):
            pass
    except IOError as ioe:
        err_msg = '\nError opening master gramox file: {}\n'
        sys.exit(err_msg.format(ioe))

    try:
        with open(args.rel_abd_fnh):
            pass
    except IOError as ioe:
        err_msg = '\nError opening relative abundance file: {}\n'
        sys.exit(err_msg.format(ioe))

# Read relative abundance data
    rel_abd_data = pd.read_csv(args.rel_abd_fnh, sep='\t')
    otus = rel_abd_data['#OTU ID']

# Read gramox data from master file
    with open(args.master_fnh, 'rU') as bgof:
        bgodata = {line.strip().split('\t')[0]: '\t'.join(line.strip().split('\t')[1:4])
                   for line in bgof.readlines()}

# Write classified gramox data to tsv file
    with open(args.out_fnh, 'w') as gramoxout:
        gramoxout.write('#OTU\tGram Status\tOxygen Requirement\tSource\n')
        for otu in otus:
            if otu in bgodata.keys():
                gramoxout.write('{}\t{}\n'.format(otu, bgodata[otu]))
            else:
                gramoxout.write('{}\n'.format(otu))

if __name__ == '__main__':
    main()
