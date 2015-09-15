#!/usr/bin/env python

'''
Abstract: Iterate through all folders and get counts of most common sequences
          found in reads.

Date: 09/03/2015

Author: Akshay Paropkari
'''

import os
import sys
import subprocess as sp
from collections import Counter
import argparse
try:
    import shlex
except ImportError as ie:
    sys.exit('Please install {} module before executing this script.'
             .format(ie))


def count_common_seqs(bc, tc, fo):
    '''
    Count the most common sequences found in your fastq files.

    :type bc: int
    :param bc: Number of bases from start to include in count results.
    :type tc: int
    :param tc: Number of most common count results to include in output.
    :type fo: str
    :param fo: Name of the fastq file being processed.
    '''

    # Get list of reads
    with open(fo, 'r') as fastq:
        seq = [line.strip()[:bc]
               for line in fastq.readlines()[1::4]]

    # Get read count from fastq file
    qcs = 'count_seqs.py -i {}'.format(fo)
    kwargs = shlex.split(qcs)
    out = sp.check_output(kwargs)
    read_count = out.split(':')[0]

    # Count of most common sequence
    result = []
    cnt = Counter(seq)
    for item in cnt.most_common(tc):
        seq_pct = 100*(item[1]/float(read_count))
        result.append((read_count, item[0], item[1], seq_pct))
    return result


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder in a file and get '
                            'counts of most common sequences found in reads.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    parser.add_argument('-nb', '--base_counts', type=int, default=15,
                        help='Number of bases from start to include in count'
                             ' result. Default is 15 bases.')
    parser.add_argument('-tcc', '--top_common_count', type=int, default=20,
                        help='Number of most common count results to '
                             'include in output. Default is top 20 most '
                             'common ("-nb" bases long) sequences.')
    parser.add_argument('out_fnh',
                        help='Output file path to save most common count '
                              'results.')
    return parser.parse_args()


def main():
    args = prog_options()

    # Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            print root[24:27]

    # Get relevant file name
            for file in files:
                if file.endswith('.fastq'):
                    trim_out = file

    # Option to get counts of adapters in FLASh output merged file
            most_common_seqs = count_common_seqs(args.base_counts,
                                                 args.top_common_count,
                                                 trim_out)
            with open(args.out_fnh, 'w') as outf:
                outf.write('{}\n'.format(trim_out))
                outf.write('Read Count\tCommon Sequence\t'
                           'Common Sequence Counts\tPercent in reads')
                for d in most_common_seqs:
                    outf.write('{}\t{}\t{}\t{}'
                               .format(d[0], d[1], d[2], d[3]))
            return

if __name__ == '__main__':
    sys.exit(main())
