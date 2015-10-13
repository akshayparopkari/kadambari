#!/usr/bin/env python


import os
import subprocess as sp
import argparse
try:
    import shlex
except ImportError as ie:
    sys.exit('Please install {} module before executing this script.'\
             .format(ie))


def prog_options():
    parser = argparse.ArgumentParser(
                description='Iterate through each folder in a file and run '
                            'some command.')
    parser.add_argument('sample_dir',
                        help='Directory containing sample folders.')
    return parser.parse_args()


def main():
    args = prog_options()

    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir and 'Data' not in root:
            os.chdir(root)
            print 'root: ', root
            print 'dir: ', dirs
            print 'files: ', files

#             if dirs == 'BaseCalls':
#             for file in files:
#                 cmd = 'mv {} ../../../{}'.format(file, file)
#                 kwargs = shlex.split(cmd)
#                 print 'kwargs', kwargs, '\n'
#                 sp.call(kwargs)

#             for file in files:
#                 if not file.endswith('fastq.gz') and not file.startswith('out') and not file.startswith('flash'):
#                     in1 = 'rm {}'.format(file)
#                     kwargs = shlex.split(in1)
#                     print kwargs
#                     sp.call(kwargs)
    return


if __name__ == '__main__':
    main()
