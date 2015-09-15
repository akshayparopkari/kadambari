#!/usr/bin/env python

'''
Abstract: Iterate through all folders and execute a cutadapt command.

Date: 08/31/2015

Author: Akshay Paropkari
'''

import os
import sys
import subprocess as sp
import argparse
try:
    import shlex
except ImportError as ie:
    sys.exit('Please install {} module before executing this script.'
             .format(ie))


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder in a file and run '
                            'cutadapt command on sample fastq files.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    parser.add_argument('adapter_list',
                        help='Path to file containing one adapter sequences '
                             'per line.')
    parser.add_argument('nos', choices=['1', '2'],
                        help='Number of samples in the fastq files.')
    return parser.parse_args()


def main():
    args = prog_options()

# Checking the necessary input arguments for validity.
    try:
        with open(args.adapter_list):
            pass
    except IOError as ioe:
        sys.exit('\nError with adapter list file: {}\n'.format(ioe))

    # Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            print root[29:33]

            for file in files:
                if file.endswith('extendedFrags.fastq'):
                    mergedfile = file

    # Read adapter sequences from file
            with open(args.adapter_list, 'r') as acf:
                for line in acf.readlines():
                    adapters = line.strip().split('\r')

    # Map seq to their names
            f_names = ['P1_V1-V3', 'P1_V4-V5', 'P2_V1-V3', 'P2_V4-V5']
            fwd_adapters = adapters[0::2]
            sp_seq = {f: s for f, s in zip(f_names, fwd_adapters)}

    # Generate and run cutadapt command to run
            for k, v in sp_seq.iteritems():
                if args.nos == '2':
                    cmd = 'cutadapt -g {} -o {}.fastq --discard-untrimmed {}'\
                          .format(v, k, mergedfile)
                    kwargs = shlex.split(cmd)
                    print kwargs
#                     out = sp.check_output(kwargs)
#                     print out
#                     with open(k+'_trimlog.txt', 'w') as outlog:
#                         outlog.write('{}'.format(out))
                elif args.nos == '1' and k in ['P1_V1-V3', 'P1_V4-V5']:
                    cmd = 'cutadapt -g {} -o {}.fastq --discard-untrimmed {}'\
                          .format(v, k, mergedfile)
                    kwargs = shlex.split(cmd)
                    out = sp.check_output(kwargs)
                    print out
                    with open(k+'_trimlog.txt', 'w') as outlog:
                        outlog.write('{}'.format(out))
    return

if __name__ == '__main__':
    sys.exit(main())
