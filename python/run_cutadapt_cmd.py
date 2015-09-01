#!/usr/bin/env python

'''
Abstract: Iterate through all folders and execute a cutadapt command.

Date: 08/31/2015

Author: Akshay Paropkari
'''

import os
import sys
import shlex
import argparse
import subprocess as sp
from collections import Counter


def count_common_seqs(bc, tc):
    '''
    Count the most common sequences found in your fastq files.

    :type bc: int
    :param bc: Number of bases from start to include in count results.
    :type tc: int
    :param tc: Number of most common count results to include in output.
    '''
    flash_out = 'out.extendedFrags.fastq'

# Get list of reads
    with open(flash_out, 'r') as fastq:
        seq = [line.strip()[:bc]
               for line in fastq.readlines()[1::4]]

# Count of most common sequence
    cnt = Counter(seq)
    return cnt.most_common(tc)


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder in a file and run '
                            'cutadapt command on sample fastq files.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    parser.add_argument('-c', '--adapter_counts', default=False,
                        help='Path to file containing one adapter sequences '
                             'per line.')
    parser.add_argument('-nb', '--base_counts', type=int, default=15,
                        help='Number of bases from start to include in count'
                             ' result. Default is 15 bases.')
    parser.add_argument('-t', '--top_common_count', type=int, default=20,
                        help='Number of most common count results to '
                             'include in output. Default is top 20 most '
                             'common ("-nb" bases long) sequences.')
    return parser.parse_args()


def main():
    args = prog_options()

# Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            print root

# Read adapter sequences from file
            with open(args.adapter_counts, 'r') as acf:
                adapters = [line.strip()
                            for line in acf.readlines()]



# Option to get counts of adapters in FLASh output merged file
            if args.adapter_counts:
                print count_common_seqs(args.base_counts,
                                        args.top_common_count)
                print


if __name__ == '__main__':
    sys.exit(main())
