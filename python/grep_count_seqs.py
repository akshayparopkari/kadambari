#!/usr/bin/env python

'''
Abstract: Iterate through all folders and grep count all adapter sequences.

Date: 09/03/2015

Author: Akshay Paropkari
'''

import os
import sys
import subprocess as sp
import argparse
try:
    import shlex
except ImportError as ie:
    sys.exit('Please install {} module before executing this script.'\
             .format(ie))


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder (sample) in a '
                            'directory and grep count all adapter sequences '
                            'in zipped fastq file.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    parser.add_argument('adapter_list',
                        help='File containing one sequence per line. This is '
                             'the sequence that will be counted by grep.')
    return parser.parse_args()


def main():
    args = prog_options()

    # Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)

            for file in files:
                if file.endswith('R1_001.fastq.gz'):
                    R1 = file
                elif file.endswith('R2_001.fastq.gz'):
                    R2 = file

    # Read adapter sequences from file
            with open(args.adapter_list, 'r') as acf:
                for line in acf.readlines():
                    adapters = line.strip().split('\r')

    # Dict mapping adapter seq to primer name
            adapter_name = {'TCGATCG' : '27F_A',
                            'CGGACTTGATGTACGA' : '519R_A',
                            'CGAGCAATCCACTC' : '515F_C',
                            'GATTAGCTGC' : '806R_C',
                            'ATCTGTCATG' : '27F_B',
                            'TCAGTAGCTACGC' : '519R_B',
                            'GATCAGTCGTCTCACTC' : '515F_D',
                            'ATCAGCA' : '806R_D'}

    # Grep adapter sequences
            for f in [R1, R2]:
                for adapter in adapters:
                    cmd = 'zgrep -c "{}" {}'.format(adapter, f)
                    kwargs = shlex.split(cmd)
                    print '{}__{}'.format(f, adapter_name[adapter])
                    sp.call(kwargs)
                print
            print

    return

if __name__ == '__main__':
    sys.exit(main())
