#!/usr/bin/env python

'''
Abstract: Iterate through all folders and execute quality filtering on
          merged-trimmed fastq files. filtering is achieved through Fastx-
          Toolkit. For more information, visit:
          http://hannonlab.cshl.edu/fastx_toolkit/index.html

Date: 10/06/2015

Author: Akshay Paropkari
'''


import sys
import argparse
from os import walk
import subprocess as sp
from os.path import join, relpath
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
                        help='Specify the minimum quality score to keep. '
                             'Default is Q20.')
    parser.add_argument('-p', '--min_base_pct', type=int, default='90',
                        help='Specify the minimum percent of bases with "-q" '
                             'quality scores. Default is 90 percent.')
    return parser.parse_args()


def main():
    args = prog_options()

# Iterate through all directories and access their files
    for root, dirs, files in walk(args.sample_dir):
        if root != args.sample_dir:
            print root

# Map primer names to their 16S regions
            gene_region = {'519R': 'V1-V3', '806R': 'V4-V5'}

# Iterate through all files and perform quality filtering
            for file in files:
                if file.endswith('R-assigned-01.fastq'):
                    fname = file.split('_')[0] + '_' +\
                            gene_region[file.split('_')[1][:4]]
                    cmd = 'fastq_quality_filter -v -q {} -p {} -i {} '\
                          '-o {}_qual_fil.fastq'\
                          .format(args.min_qual_score, args.min_base_pct,
                                  relpath(join(root, file)),
                                  relpath(join(root, fname)))
                    kwargs = shlex.split(cmd)
                    print kwargs, '\n'
                    out = sp.check_output(kwargs)
                    print out
                    with open(relpath(join(root, fname+'_qual_fil.txt')), 'w')as fo:
                        fo.write('{}\n{}'.format(fname, out))
    return

if __name__ == '__main__':
    sys.exit(main())
