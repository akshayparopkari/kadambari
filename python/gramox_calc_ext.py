#!usr/bin/env python
"""
Abstract: Get gram stain and oxygen requirement (gramox) data by category.

Date: 03/11/2015

Author: Akshay Paropkari
"""
import argparse
import sys
from collections import defaultdict


def categorize_otus(catname, catabr, catdata):
    """
    Categorize OTU's based on their gram stain and oxygen requirement.

    :type catname: str
    :param catname: List of string, with all smoking categories.
    :type colorby: str
    :param colorby: Abbreviation of category name.
    :type catdata: defaultdict(list)
    :param catdata: OTU data obtained from overall gramox data file.
    """
    cats = defaultdict(list)
    gramox_cat = ['_gnan', '_gna', '_gpan', '_gpa', '_unk']

    Gram_Negative_Anaerobes = ['0', '0']
    Gram_Negative_Aerobes = ['0', '1']
    Gram_Positive_Anaerobes = ['1', '0']
    Gram_Positive_Aerobes = ['1', '1']

    for k, v in catdata[catname].iteritems():
        if v == Gram_Negative_Anaerobes:
            cats[catabr+gramox_cat[0]].append(k)
        elif v == Gram_Negative_Aerobes:
            cats[catabr+gramox_cat[1]].append(k)
        elif v == Gram_Positive_Anaerobes:
            cats[catabr+gramox_cat[2]].append(k)
        elif v == Gram_Positive_Aerobes:
            cats[catabr+gramox_cat[3]].append(k)
        elif 'NA' in v:
            cats[catabr+gramox_cat[4]].append(k)
    return [cats[catabr+gramox_cat[0]], cats[catabr+gramox_cat[1]],
            cats[catabr+gramox_cat[2]], cats[catabr+gramox_cat[3]],
            cats[catabr+gramox_cat[4]]]


def handle_program_options():
    parser = argparse.ArgumentParser(
        description='Get categorized gramox data for OTU\'s from annotated '
                    'relative abundance file.'
    )
    parser.add_argument('rel_abd_fnh',
                        help='Path to relative abundance data file. This file'
                             ' is the transposed output from '
                             'biom_relative_abundance and its second column '
                             'must be the category column.')
    parser.add_argument('in_gramox_fnh',
                        help='Path to completed gramox file. All entries '
                             'must be present, unidentified entries shold be '
                             'labeled "NA".')
    parser.add_argument('out_gramox_fnh',
                        help='Categorized gramox data will be written to '
                             'this file.')
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.rel_abd_fnh):
            pass
    except IOError as ioe:
        err_msg = '\nError opening relative abundance file: {}\n'
        sys.exit(err_msg.format(ioe))

    try:
        with open(args.gramox_fnh):
            pass
    except IOError as ioe:
        err_msg = '\nError opening completed gramox file: {}\n'
        sys.exit(err_msg.format(ioe))

# Read gramox data.
    with open(args.gramox_fnh, 'rU') as gdf:
        gramox_data = {line.strip().split('\t')[0]: line.strip().split('\t')[1:3]
                       for line in gdf.readlines()[1:]}

# Get the bacteria based on categories.
    with open(args.rel_abd_fnh, 'rU') as bacat:
        bacteria = [entry for entry in bacat.readline().strip().split('\t')[2:]]
        cat = {'Cat1': {}, 'Cat2': {},
               'Cat3': {}, 'Cat4': {}}
        for line in bacat:
            line = line.strip().split('\t')[1:]
            if line[0] == 'Cat1':
                for i, data in enumerate(line[1:]):
                    if float(data) > 0:
                        cat['Cat1'][bacteria[i]] = gramox_data[bacteria[i]]
            if line[0] == 'Cat2':
                for i, data in enumerate(line[1:]):
                    if float(data) > 0:
                        cat['Cat2'][bacteria[i]] = gramox_data[bacteria[i]]
            if line[0] == 'Cat3':
                for i, data in enumerate(line[1:]):
                    if float(data) > 0:
                        cat['Cat3'][bacteria[i]] = gramox_data[bacteria[i]]
            if line[0] == 'Cat4':
                for i, data in enumerate(line[1:]):
                    if float(data) > 0:
                        cat['Cat4'][bacteria[i]] = gramox_data[bacteria[i]]

# Obtain numbers for all gramox categories.
    cat1_data = categorize_otus('Cat1', 'c1', cat)
    cat2_data = categorize_otus('Cat2', 'c2', cat)
    cat3_data = categorize_otus('Cat3', 'c3', cat)
    cat4_data = categorize_otus('Cat4', 'c4', cat)

# Write to output file.
    with open(args.out_gramox_fnh, 'w') as outf:
        outf.write('Data for each smoking category.\n')
        outf.write('\nCat1\nGram Negative Anaerobes\t{}\n'
                   'Gram Negative Aerobes\t{}\nGram Positive Anaerobes\t{}\n'
                   'Gram Positive Aerobes\t{}\nUnknown\t{}\n'
                   .format(len(cat1_data[0]), len(cat1_data[1]),
                           len(cat1_data[2]), len(cat1_data[3]),
                           len(cat1_data[4])))
        outf.write('\nCat2\nGram Negative Anaerobes\t{}\n'
                   'Gram Negative Aerobes\t{}\nGram Positive Anaerobes\t{}\n'
                   'Gram Positive Aerobes\t{}\nUnknown\t{}\n'
                   .format(len(cat2_data[0]), len(cat2_data[1]),
                           len(cat2_data[2]), len(cat2_data[3]),
                           len(cat2_data[4])))
        outf.write('\nCat3\nGram Negative Anaerobes\t{}\n'
                   'Gram Negative Aerobes\t{}\nGram Positive Anaerobes\t{}\n'
                   'Gram Positive Aerobes\t{}\nUnknown\t{}\n'
                   .format(len(cat3_data[0]), len(cat3_data[1]),
                           len(cat3_data[2]), len(cat3_data[3]),
                           len(cat3_data[4])))
        outf.write('\nCat4\nGram Negative Anaerobes\t{}\n'
                   'Gram Negative Aerobes\t{}\nGram Positive Anaerobes\t{}\n'
                   'Gram Positive Aerobes\t{}\nUnknown\t{}\n'
                   .format(len(cat4_data[0]), len(cat4_data[1]),
                           len(cat4_data[2]), len(cat4_data[3]),
                           len(cat4_data[4])))

if __name__ == '__main__':
    main()
