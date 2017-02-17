#!/usr/bin/env python
"""
:Abstract: Given an abundance table, calculate pairwise Spearman's Rho of the columns.
:Date: 02/15/2017
:Author: Akshay Paropkari
"""

import sys
import argparse
from itertools import combinations
from phylotoast.otu_calc import otu_name
from phylotoast.util import parse_map_file, gather_categories
from phylotoast.biom_calc import arcsine_sqrt_transform as ast, relative_abundance
importerrors = set()
try:
    import biom
except:
    importerrors.add("biom")
try:
    import numpy as np
except:
    importerrors.add("numpy")
try:
    from scipy.stats import rankdata, spearmanr as Spearman, kendalltau as KendallTau
except ImportError as ie:
    importerrors.add("scipy")
try:
    from statsmodels.sandbox.stats.multicomp import multipletests as MT
except:
    importerrors.add("statsmodels")
if len(importerrors) > 0:
    for err in importerrors:
        print("Please install missing module: {}".format(ie))
    sys.exit()


def calc_corr(biom_file, abd_data, otu_combos, gather_categories, method):
    """
    This function calculates the pairwaise spearman correlation for each group.
    """
    result = []
    for combo in otu_combos:
        a = []
        b = []
        otu1 = otu_name(biom_file.metadata(combo[0], "observation")["taxonomy"])
        otu2 = otu_name(biom_file.metadata(combo[1], "observation")["taxonomy"])
        for cat in gather_categories.keys():
            for sid in gather_categories[cat].sids:
                a.append(abd_data[sid][combo[0]])
                b.append(abd_data[sid][combo[1]])
            rho_corr, rho_p_val = Spearman(rankdata(a, "ordinal"), rankdata(b, "ordinal"))
            kendalltau_corr, kendalltau_p_val = KendallTau(rankdata(a, "ordinal"),
                                                           rankdata(b, "ordinal"))
            result.append([cat, otu1, otu2, rho_corr, rho_p_val])

    # Perform multiple testing correction
    p_vals = [entry[-1] for entry in result]
    hypothesis, pvals_corrected, alphacSidak, alphacBonf = MT(p_vals, method=method)
    idx_to_retain = np.where(np.array(pvals_corrected) < 0.05)[0].tolist()
    print("{} correlations removed through multiple test correction."
          .format(len(result) - len(idx_to_retain)))
    updated_result = []
    for idx in idx_to_retain:
        updated_p_val = pvals_corrected[idx]
        updated_result.append((result[idx][0], result[idx][1], result[idx][2],
                               result[idx][3], updated_p_val,))
    return updated_result


def program_options():
    parser = argparse.ArgumentParser(
        description=("Given an abundance table, calculate pairwise Spearman's Rho of the "
                     "columns. Results undergo multiple testing correction (BH)."))
    parser.add_argument("in_biomf", help="Input abundance file containing OTU names as "
                        "columns and SampleIDs as rows. Ideally, this is the output from "
                        "biom_relative_abundance.py script.")
    parser.add_argument("map_fnh", help="Mapping file associated with input BIOM file.")
    parser.add_argument("category_column", help="Column name from mapping file which is "
                        "associated with categories.")
    parser.add_argument("out_fnh", help="Path and name to output correlation matrix file."
                        " The format for the tab-separated file will be: Category -> "
                        "Variable -> by Variable -> Correlation")
    parser.add_argument("-mt", "--multiple_test_method", default="fdr_bh",
                        choices={"bonferroni", "holm", "fdr_bh", "fdr_by"},
                        help="Method to be used for testing and adjusting p-values. By "
                        "default, it is set to Benjamini/Hochberg. Please refer to http:/"
                        "/statsmodels.sourceforge.net/stable/generated/statsmodels."
                        "sandbox.stats.multicomp.multipletests.html for additional "
                        "information.")
    return parser.parse_args()


def main():
    args = program_options()

    try:
        biomf = biom.load_table(args.in_biomf)
    except IOError as ioe:
        sys.exit("Error with input BIOM format file: {}".format(ioe))
    rel_abd = relative_abundance(biomf)
    ast_rel_abd = ast(rel_abd)

    try:
        mheader, mdata = parse_map_file(args.map_fnh)
    except IOError as ioe:
        sys.exit("Error with input mapping file: {}".format(ioe))

    # Gather sampleID categories
    sid_cat = gather_categories(mdata, mheader, [args.category_column])

    # Get pairwise combinations of OTUs
    otu_combos = combinations(biomf.ids("observation"), 2)

    # Get correlation calculation results
    corr_results = calc_corr(biomf, ast_rel_abd, otu_combos, sid_cat,
                             args.multiple_test_method)

    with open(args.out_fnh, "w") as outf:
        outf.write("Category\tVariable\tby Variable\tCorrelation\tp value\n")
        for k in corr_results:
            outf.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(k[0], k[1], k[2], k[3], k[4]))


if __name__ == "__main__":
    sys.exit(main())
