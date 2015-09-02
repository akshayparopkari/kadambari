#!/usr/bin/env python

'''
Abstract: Iterate through all folders and execute a cutadapt command.

Date: 08/31/2015

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
    sys.exit('Please install {} module before executing this script.'\
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


def trim_adapter_seqs(ac, fo):
    '''
    Trim adapter sequences ligated either at 5'-end or 3'-end of the reads.

    :type ac: file name
    :param bc: File containing one adapter sequence per line.
    :type fo: str
    :param fo: Name of the fastq file being processed.
    '''

# Checking the necessary input arguments for validity.
    try:
        with open(ac):
            pass
    except IOError as ioe:
        sys.exit('\nError with adapter list file: {}\n'.format(ioe))

    # Read adapter sequences from file
    with open(ac, 'r') as acf:
        for line in acf.readlines():
            adapters = line.strip().split('\r')

# Generate cutadapt command to run
    s = ''
    for adapter in adapters:
        if args.ligation_end == 5:
            s += ' -g ' + adapter
        else:
            s += ' -a ' + adapter
    cmd = 'cutadapt{} -o trimmed.fastq {}'.format(s, flash_out)

# Run cutadapt command
    kwargs = shlex.split(cmd)
    out = sp.check_output(kwargs)
#             print out
    with open('trimmed.txt', 'w') as outlog:
        outlog.write('{}'.format(out))
    return out


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder in a file and run '
                            'cutadapt command on sample fastq files.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    parser.add_argument('-t', '--trim_adapters', default=False,
                        help='Option to trim adapter sequences.')
    parser.add_argument('-le', '--ligation_end', type=int, choices=[5, 3],
                        help="Specify whether adapter sequence was ligated "
                             "to 5'-end (5) or 3'-end (3) of the reads.")
    parser.add_argument('-al', '--adapter_list', default=None,
                        help='Path to file containing one adapter sequences '
                             'per line.')
    parser.add_argument('-c', '--counts', default=False,
                        help='Option to get most common sequence counts.')
    parser.add_argument('-nb', '--base_counts', type=int, default=15,
                        help='Number of bases from start to include in count'
                             ' result. Default is 15 bases.')
    parser.add_argument('-tcc', '--top_common_count', type=int, default=20,
                        help='Number of most common count results to '
                             'include in output. Default is top 20 most '
                             'common ("-nb" bases long) sequences.')
    parser.add_argument('-mco', '--most_common_out',
                         help='Output file path to save most common count '
                              'results.')
    return parser.parse_args()


def main():
    args = prog_options()

    # Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            sample = root[24:27]
            print sample

            for file in files:
                if file.endswith('extendedFrags.fastq'):
                    flash_out = file
                elif file.endswith('trimmed.fastq'):
                    trim_out = file

    # Option to trim adapters from sequences
            if args.trim_adapters:
                print trim_adapter_seqs(args.adapter_list, flash_out)
                print

    # Option to get counts of adapters in FLASh output merged file
            if args.counts:
                most_common_seqs = count_common_seqs(args.base_counts,
                                                     args.top_common_count,
                                                     trim_out)
                with open(args.most_common_out, 'w') as outf:
                    outf.write('Read Count\tCommon Sequence\t'
                               'Common Sequence Counts\tPercent in reads')
                    for d in most_common_seqs:
                        outf.write('{}\t{}\t{}\t{}'\
                                   .format(d[0], d[1], d[2], d[3]))
    return

if __name__ == '__main__':
    sys.exit(main())
