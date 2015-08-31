#!/usr/bin/env python

'''
Abstract: Iterate through all folders and execute a shell/python command.

Date: 08/06/2015

Author: Akshay Paropkari
'''

import os
import subprocess as sp
import shlex
import argparse


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder in a file and run '
                            'FLASh command on sample fastq files.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    return parser.parse_args()


def main():
    args = prog_options()

    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            print root

            for file in files:
                if file.endswith('R1_001.fastq.gz'):
                    R1 = file
                else:
                    R2 = file
            in1 = 'flash {} {} -M 300 --cap-mismatch-quals'.format(R1, R2)
            kwargs = shlex.split(in1)
            out = sp.check_output(kwargs)
            print out
            with open('flash.txt', 'w') as outlog:
                outlog.write('{}'.format(out))
    return


if __name__ == '__main__':
    main()
