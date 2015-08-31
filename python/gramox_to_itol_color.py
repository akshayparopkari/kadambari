#!/usr/bin/env python
"""
Abstract: Python script that creates color definition file for iTol tree
          displaying gram status and oxygen requirement of all the OTU's.

Date: 12/31/2014
"""
import sys
import argparse
from argparse import RawTextHelpFormatter as rthf


def handle_program_options():
    parser = argparse.ArgumentParser(
        description='Create a output color file for iTol images. The color\n'
        'file assigns colors to each OTU based on their gram status and oxygen\n'
        'requirement. This helps in identifying the composition of OTU in a \n'
        'given sample.',
        formatter_class=rthf
        )
    parser.add_argument(
        'in_fnh',
        help='Input text file containing OTU, gram status and oxygen requirement.\n'
        'The first line must be the header line and all values must be tab\n'
        'seprateed. All unknown values must have default value of \'NA\'.\n'
        'An example of the format:\n'
        '#OTU    Gram Status    Oxygen Requirement\n'
        'Abiotrophia_defectiva    1    1\n'
        'Leptotrichia_hongkongensis    NA    1\n'
        )
    parser.add_argument('out_fnh',
                        help='Path to the output text file.')
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.in_fnh):
            pass
    except IOError as ioe:
        err_msg = '\nError opening input text file (in_fnh): {}\n'
        sys.exit(err_msg.format(ioe))

# Define colors and description for each combination of gram status and oxygen requirement.
    phen_colors = {'11': '#0868ac', '10': '#43a2ca',
                   '01': '#7bccc4', '00': '#bae4bc'}    # 5-class GnBu sequential
    phen_names = {'11': 'Gram Positive Aerobes', '10': 'Gram Positive Anaerobes',
                  '01': 'Gram Negative Aerobes', '00': 'Gram Negative Anaerobes'}

# Calculate the category for each OTU.
    with open(args.in_fnh, 'rU') as inF:
        gram_ox = {}
        for line in inF.readlines()[1:]:
            line = line.strip().split('\t')
            genus, gram, ox = line
            if 'NA' in line:
                gram_ox[genus] = 'Unknown'
            else:
                gram_ox[genus] = gram + ox

# Write output to file.
    with open(args.out_fnh, 'w') as outF:
        for otu, info in sorted(gram_ox.iteritems()):
            if info != 'Unknown':
                data = '\t'.join([otu, 'range', phen_colors[gram_ox[otu]],
                                  phen_names[gram_ox[otu]]])
            else:
                data = '\t'.join([otu, 'range', '#AAAAAA', 'Unknown'])
            outF.write(data+'\n')

if __name__ == '__main__':
    sys.exit(main())
