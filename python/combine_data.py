#!/usr/bin/env python
"""
Abstract: Script to merge multiple raw abundance files. This is a slight
          variation of primer_average.py of PhyloToAST.

Date: 01/05/15

Author: Akshay Paropkari
"""
import sys
import argparse


def handle_program_options():
    parser = argparse.ArgumentParser(
        description='Combine raw abundance files, which were generated for '
        'different categories in order to identity unique OTU\'s and their '
        'abundances for each sample of each category, into a single file.'
    )
    parser.add_argument('in1_fnh',
                        help='First input file handle name.')
    parser.add_argument('in2_fnh',
                        help='Second input file handle name.')
    parser.add_argument('out_fnh',
                        help='Output file handle name.')
    return parser.parse_args()


def main():
    args = handle_program_options()

# Checking the necessary input arguments for validity.
    try:
        with open(args.in1_fnh):
            pass
    except IOError as ioe:
        sys.exit('\nError with first input file:{}\n'.format(ioe))

    try:
        with open(args.in2_fnh):
            pass
    except IOError as ioe:
        sys.exit('\nError with second input file:{}\n'.format(ioe))

    try:
        with open(args.out_fnh, 'w'):
            pass
    except IOError as ioe:
        sys.exit('{}: Path to output file is invalid.'.format(ioe))

# Creating dictionaries keyed on OTU and rel abd as for each sample as their
# values. Each dictionary represents data from a single input file.
    data_1 = {}
    with open(args.in1_fnh, 'rU') as f:
        data1 = [line.strip().split('\t') for line in f]
    len_data1 = len(data1[1][1:])
    for item in data1:
        data_1[item[0]] = '\t'.join(item[1:])
    header1 = data1[0][0]+'\t' + data_1[data1[0][0]]
    data1_zeros = '\t'.join(['0'] * len_data1)

    data_2 = {}
    with open(args.in2_fnh, 'rU') as f:
        data2 = [line.strip().split('\t') for line in f]
    len_data2 = len(data2[1][1:])
    for item in data2:
        data_2[item[0]] = '\t'.join(item[1:])
    header2 = data_2[data2[0][0]]
    data2_zeros = '\t'.join(['0'] * len_data2)

# Creating dict of shared OTU and their rel abd for both input files.
    otuid1 = set([item[0] for item in data1[1:]])
    otuid2 = set([item[0] for item in data2[2:]])
    shared = list(otuid1 & otuid2)

    shared_data1 = {}
    for k, v in data_1.iteritems():
        if k in shared:
            shared_data1[k] = v

    shared_data2 = {}
    for k, v in data_2.iteritems():
        if k in shared:
            shared_data2[k] = v

# Creating dict of exclusive OTU and their rel abd for each file.
    only_data1_otu = list(otuid1.difference(otuid2))
    only_data2_otu = list(otuid2.difference(otuid1))

    only_data1 = {}
    for k, v in data_1.iteritems():
        if k in only_data1_otu:
            only_data1[k] = v
    only_data2 = {}
    for k, v in data_2.iteritems():
        if k in only_data2_otu:
            only_data2[k] = v

# Writing to output file.
    with open(args.out_fnh, 'w') as f:
        f.write('{}\t{}\n'.format(header1, header2))
        for k1, v1 in shared_data1.iteritems():
            for k2, v2 in shared_data2.iteritems():
                if k1 == k2:
                    f.write('{}\t{}\t{}\n'.format(k1, v1, v2))
        for k, v in only_data1.iteritems():
            f.write('{}\t{}\t{}\n'.format(k, v, data2_zeros))
        for k, v in only_data2.iteritems():
            f.write('{}\t{}\t{}\n'.format(k, data1_zeros, v))

if __name__ == '__main__':
    main()
