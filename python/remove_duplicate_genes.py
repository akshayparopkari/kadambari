#!/usr/bin/env python
"""
Abstract: Remove duplicate gene occurrences and output into a tsv file.

Author: Akshay Paropkari
"""
import sys
import argparse


def handle_program_options():
    parser = argparse.ArgumentParser(
        description='Remove all duplicate occurences of gene functions and '
                    'their relative abundaces and output to a tsv file. The '
                    'output file will contain unique functions and their '
                    'relative abundances for each patient. Note that it is '
                    'assumed that all duplicate occurences of relative '
                    'abundances are repeated number.'
    )
    parser.add_argument('in_fnh',
                        help='Input file containing repeated gene functions '
                             'and their relative abundances for all patients.'
                        )
    parser.add_argument('out_fnh',
                        help='Path to the output file. (Required input)')
    return parser.parse_args()


def main():
    args = handle_program_options()

# Checking the necessary input arguments for validity.
    try:
        with open(args.in_fnh):
            pass
    except IOError as ioe:
        sys.exit('\nError with first input file:{}\n'.format(ioe))

    data = []
    with open(args.in_fnh, 'rU') as f:
        for line in f:
            line = line.strip().split('\t')
            data.append(line)
    header = '\t'.join(data[0])

    func = {}
    for entry in data[1:]:
        func[entry[0]] = '\t'.join(entry[1:])

    with open(args.out_fnh, 'w') as f:
        f.write('{}\n'.format(header))
        for k, v in func.iteritems():
            f.write('{}\t{}\n'.format(k, v))

if __name__ == '__main__':
    main()
