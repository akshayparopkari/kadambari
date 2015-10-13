#!/usr/bin/env python

'''
Abstract: Iterate through all folders and perform adapter trimming on all
          reads. This script uses a sequencing mapping file to map wells to
          sampleids and trim 3'-end adapters of reads. Additionally, reads with
          lengths <200 nt are discarded. For more information on Skewer, visit
          http://www.biomedcentral.com/1471-2105/15/182

Date: 10/7/2015

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
                            '3\'-end adapters sample fastq files via Skewer. '
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
                             'cells.Please note that reverse adapters must be '
                             'reverse complemented wherever necessary. Format '
                             'of the file: well->sampleid->adapters.')
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
            well = root
            print '\n', well

# Read adapter sequences from file
            with open(args.ap_map, 'rU') as apf:
                adapters = {line.strip().split('\t')[1]: line.strip().split('\t')[0]
                            for line in apf.readlines()[1:]}

# Read mapping file data
            with open(args.ws_map, 'rU') as csvf:
                map_data = [line
                            for line in DictReader(csvf, delimiter='\t')]

# Calculate skewer command for each sample
            for sample in map_data:
                if sample['well'] == well.split('-')[0].split('/')[5]:
                    fname1 = sample['sampleid']+'_27F-assigned-01.fastq'
                    fname2 = sample['sampleid']+'_515F-assigned-01.fastq'

                    # Plate 1 trimming
                    for file in files:
                        if file == fname1:
                            if sample['519R_A'] != '-':
                                cmd1 = 'skewer -t {} -b -m tail -l 200 -x {} \
                                       -o {}_519R {}'\
                                       .format(args.threads, sample['519R_A'],
                                               file.split('_')[0], file)
                                kwargs1 = shlex.split(cmd1)
                                print 'Well: {} | Primer: {} | SampleID: {}'\
                                      .format(sample['well'],
                                              adapters[sample['519R_A']],
                                              file.split('_')[0])
                                print kwargs1, '\n'
                                out1 = sp.check_output(kwargs1)
                                print out1
                            elif sample['519R_B'] != '-':
                                cmd3 = 'skewer -t {} -b -m tail -l 200 -x {} \
                                       -o {}_519R {}'\
                                       .format(args.threads, sample['519R_B'],
                                               file.split('_')[0], file)
                                kwargs3 = shlex.split(cmd3)
                                print 'Well: {} | Primer: {} | SampleID: {}'\
                                      .format(sample['well'],
                                              adapters[sample['519R_B']],
                                              file.split('_')[0])
                                print kwargs3, '\n'
                                out3 = sp.check_output(kwargs3)
                                print out3

                        elif file == fname2:
                            if sample['806R_C'] != '-':
                                cmd2 = 'skewer -t {} -b -m tail -l 200 -x {} \
                                       -o {}_806R {}'\
                                       .format(args.threads, sample['806R_C'],
                                               file.split('_')[0], file)
                                kwargs2 = shlex.split(cmd2)
                                print 'Well: {} | SampleID: {} | Primer: {}'\
                                      .format(sample['well'],
                                              adapters[sample['806R_C']],
                                              file.split('_')[0])
                                print kwargs2, '\n'
                                out2 = sp.check_output(kwargs2)
                                print out2
                            elif sample['806R_D'] != '-':
                                cmd4 = 'skewer -t {} -b -m tail -l 200 -x {} \
                                       -o {}_806R {}'\
                                       .format(args.threads, sample['806R_D'],
                                               file.split('_')[0], file)
                                kwargs4 = shlex.split(cmd4)
                                print 'Well: {} | Primer: {} | SampleID: {}'\
                                      .format(sample['well'],
                                              adapters[sample['806R_D']],
                                              file.split('_')[0])
                                print kwargs4, '\n'
                                out4 = sp.check_output(kwargs4)
                                print out4
    return

if __name__ == '__main__':
    sys.exit(main())
