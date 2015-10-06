#!/usr/bin/env python

'''
Abstract: Iterate through all folders and execute quality filtering on
          merged-trimmed fastq files.

Date: 10/06/2015

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
                            'fastx quality filtering on each merged-trimmed '
                            'fastq file. The output stats are written to '
                            'STDOUT and saved in text file.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    parser.add_argument('-q', '--min_qual_score', type=int, default='20',
                        help='Specify the minimum quality score to keep.'
                             'Default is Q20.')
    parser.add_argument('-p', '--min_base_pct', type=int, default='90',
                        help='Specify the minimum percent of bases with "-q" '
                             'quality scores. Default is 90 percent.')
    return parser.parse_args()


def main():
    args = prog_options()

    # Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            print root[24:27]

            for file in files:
                if file.endswith('R-assigned-01.fastq'):
                    fileout = raw_input('For {} input file, please enter '
                                        'output file name: '.format(file))
                    cmd = 'fastq_quality_filter -v -q {} -p {} -i {} '\
                          '-o {}.fastq'.format(args.min_qual_score,
                                               args.min_base_pct, file,
                                               fileout)
                    kwargs = shlex.split(cmd)
                    print kwargs, '\n'
                    out = sp.check_output(kwargs)
                    print out
                    with open(fileout+'.txt', 'w') as outf:
                        outf.write('{}\n{}'.format(fileout, out))
    return

if __name__ == '__main__':
    sys.exit(main())
