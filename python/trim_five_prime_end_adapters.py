#!/usr/bin/env python

'''
Abstract: Iterate through all folders and perform adapter trimming on all
          reads. This script uses a sequencing mapping file to map wells to
          sampleids and trim 5'-end adapters of reads. Additionally, reads with
          lengths <200 nt are discarded. For more information on Skewer, visit
          http://www.biomedcentral.com/1471-2105/15/182

Date: 08/31/2015

Author: Akshay Paropkari
'''

import os
import sys
import argparse
import subprocess as sp
from csv import DictReader
try:
    import shlex
except ImportError as ie:
    sys.exit('Please install {} module before executing this script.'
             .format(ie))


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder in a file and trim '
                            '5\'-end adapters sample fastq files via Skewer. '
                            'Using adapter-to-primer mapping (adapters) and '
                            'well-to-sampleid mapping (ws_map), trim command '
                            'and output files are generated. This scripts '
                            'discards reads with length <200 nt.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    parser.add_argument('ap_map',
                        help='Path to tab-delimited file containing '
                             'adapter-primer pair per line. Reverse primer '
                             'adapters must be reverse complemented entries.')
    parser.add_argument('ws_map',
                        help='Path to tab-separated sample ID to adapter '
                             'mapping file. A "-" must be used for empty '
                             'cells. Format of the file: '
                             'well->sampleid->adapters')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Optionally, you may specify number of threads to'
                             ' be used for trimming. Default is 4.')
    return parser.parse_args()


def main():
    args = prog_options()

# Checking the necessary input arguments for validity.
    try:
        with open(args.ap_map):
            pass
    except IOError as ioe:
        sys.exit('\nError with adapter file: {}\n'.format(ioe))

    try:
        with open(args.ws_map):
            pass
    except IOError as ioe:
        sys.exit('\nError with the well-to-sampleid mapping file: {}\n'
                 .format(ioe))

# Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            well = root[24:27]
            print '\n', well

# Read adapter sequences from file
            with open(args.ap_map, 'rU') as apf:
                adapters = {line.strip().split('\t')[0]: line.strip().split('\t')[1]
                            for line in apf.readlines()[1:]}

# Read mapping file data
            with open(args.ws_map, 'rU') as csvf:
                map_data = [line
                            for line in DictReader(csvf, delimiter='\t')]

# Read FLASh merged output file
            for file in files:
                if file.endswith('extendedFrags.fastq'):
                    mergedfile = file

# Calculate skewer command for each sample
            for sample in map_data:
                if sample['well'] == well:
                    # Plate 1 trimming
                    cmd1 = 'skewer -t {} -b -m head -l 200 -x {} -o {}_27F {}'\
                           .format(args.threads, sample['27F_P1'],
                                   sample['P1_sample'], mergedfile)       # P1
                    kwargs1 = shlex.split(cmd1)
                    print 'Well: {} | Primer: {} | SampleID: {}'\
                          .format(sample['well'], adapters[sample['27F_P1']],
                                  sample['P1_sample'])                    # P1
                    print kwargs1, '\n'
                    out1 = sp.check_output(kwargs1)
                    print out1

                    cmd2 = 'skewer -t {} -b -m head -l 200 -x {} -o {}_515F {}'\
                           .format(args.threads, sample['515F_P1'],
                                   sample['P1_sample'], mergedfile)      # P1
                    kwargs2 = shlex.split(cmd2)
                    print 'Well: {} | SampleID: {} | Primer: {}'\
                          .format(sample['well'], adapters[sample['515F_P1']],
                                  sample['P1_sample'])                    # P1
                    print kwargs2, '\n'
                    out2 = sp.check_output(kwargs2)
                    print out2

                    # Plate 2 trimming
                    if sample['27F_P2'] != '-':
                        cmd3 = 'skewer -t {} -b -m head -l 200 -x {} \
                               -o {}_27F {}'.format(args.threads,
                                                    sample['27F_P2'],
                                                    sample['P2_sample'],
                                                    mergedfile)           # P2
                        kwargs3 = shlex.split(cmd3)
                        print 'Well: {} | Primer: {} | SampleID: {}'\
                              .format(sample['well'],
                                      adapters[sample['27F_P2']],
                                      sample['P2_sample'])                # P2
                        print kwargs3, '\n'
                        out3 = sp.check_output(kwargs3)
                        print out3

                    if sample['515F_P2'] != '-':
                        cmd4 = 'skewer -t {} -b -m head -l 200 -x {} \
                                -o {}_515F {}'.format(args.threads,
                                                      sample['515F_P2'],
                                                      sample['P2_sample'],
                                                      mergedfile)         # P2
                        kwargs4 = shlex.split(cmd4)
                        print 'Well: {} | Primer: {} | SampleID: {}'\
                              .format(sample['well'],
                                      adapters[sample['515F_P2']],
                                      sample['P2_sample'])                # P2
                        print kwargs4, '\n'
                        out4 = sp.check_output(kwargs4)
                        print out4
    return

if __name__ == '__main__':
    sys.exit(main())
