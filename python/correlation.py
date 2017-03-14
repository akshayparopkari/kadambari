#!/usr/bin/env python
"""
:Abstract: Given an abundance table, calculate pairwise Spearman's Rho of the columns.
:Date: 02/15/2017
:Author: Akshay Paropkari
"""

import sys
import argparse
from traceback import format_exc
from time import localtime, strftime
from itertools import combinations
from phylotoast.util import parse_map_file, gather_categories
from phylotoast.biom_calc import arcsine_sqrt_transform as ast, relative_abundance
importerrors = set()
try:
    import biom
except ImportError:
    importerrors.add("biom")
try:
    import numpy as np
except ImportError:
    importerrors.add("numpy")
try:
    from multiprocessing import Pool, current_process
except ImportError:
    importerrors.add("multiprocessing")
try:
    from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as FDR
except ImportError:
    importerrors.add("statsmodels")
try:
    from scipy.stats import rankdata, spearmanr as Spearman, kendalltau as KendallTau
except ImportError:
    importerrors.add("scipy")
if len(importerrors) > 0:
    for err in importerrors:
        print("Please install missing module: {}".format(err))
    sys.exit()


def calc_corr(otu_combos, gather_categories, abd_data):
    """This function calculates the pairwise spearman correlation for each group."""
    starttime = strftime("%x | %X".format(localtime))
    print("{0} | Core-{1}".format(starttime, current_process().name.split("-")[1]))
    result = []
    for combo in otu_combos:
        list1 = []
        list2 = []
        otu1 = combo[0]
        otu2 = combo[1]
        for cat in gather_categories.keys():
            for sid in gather_categories[cat].sids:
                list1.append(abd_data[sid][otu1])
                list2.append(abd_data[sid][otu2])
        try:
            assert len(list1) == len(list2)
        except AssertionError:
            sys.exit("Abundance array lengths do not match.")
        rho_corr, rho_p_val = Spearman(rankdata(list1, "ordinal"),
                                       rankdata(list2, "ordinal"))
        result.append(["Spearman", cat, otu1, otu2, rho_corr, rho_p_val])
        kendall_corr, kendall_p_val = KendallTau(rankdata(list1, "ordinal"),
                                                 rankdata(list2, "ordinal"))
        result.append(["Kendall", cat, otu1, otu2, kendall_corr, kendall_p_val])
    return result


def calc_corr_helper(args):
    """Helper function to use Pool.map() one-argument worker functions."""
    return calc_corr(*args)


def run_fdr(result_list):
    """Run BH FDR correction on correlation results."""
    try:
        assert all(isinstance(x, list) for x in result_list) is True
    except AssertionError:
        sys.exit("Data format for FDR Correction not compatible. Please check.")
    else:
        p_vals = [entry[-1] for entry in result_list]
        hypothesis, pvals_corrected = FDR(p_vals)
        idx_to_retain = np.where(np.array(pvals_corrected) < 0.05)[0].tolist()
        print("{} correlations removed through multiple test correction."
              .format(len(result_list) - len(idx_to_retain)))
    updated_result = []
    for idx in idx_to_retain:
        updated_p_val = pvals_corrected[idx]
        updated_result.append((result_list[idx][0], result_list[idx][1],
                               result_list[idx][2], result_list[idx][3],
                               result_list[idx][4], updated_p_val,))
    return updated_result


def program_options():
    """Function to house the user inputs."""
    parser = argparse.ArgumentParser(
        description=("Given an abundance table, calculate and return FDR corrected "
                     "significant pairwise Spearman's Rho."))
    parser.add_argument("in_biomf", help="Input abundance file containing OTU names as "
                        "columns and SampleIDs as rows. Ideally, this is the output from "
                        "biom_relative_abundance.py script.")
    parser.add_argument("map_fnh", help="Mapping file associated with input BIOM file.")
    parser.add_argument("category_column", help="Column name from mapping file which is "
                        "associated with categories.")
    parser.add_argument("out_fnh", help="Path and name to output correlation matrix file."
                        " The format for the tab-separated file will be: Category -> "
                        "Variable -> by Variable -> Correlation")
    return parser.parse_args()


def main():
    args = program_options()

    try:
        biomf = biom.load_table(args.in_biomf)
    except IOError as ioe:
        sys.exit("Error with input BIOM format file: {}".format(ioe))
    else:
        rel_abd = relative_abundance(biomf)
        ast_rel_abd = ast(rel_abd)
        # Get pairwise combinations of OTUs
        otu_combos = list(combinations(biomf.ids("observation"), 2))

    try:
        mheader, mdata = parse_map_file(args.map_fnh)
    except IOError as ioe:
        sys.exit("Error with input mapping file: {}".format(ioe))
    else:
        # Gather sampleID categories
        sid_cat = gather_categories(mdata, mheader, [args.category_column])

    # Create arguments for helper function to be supplied to multiprocessing pool.map()
    chunksize = 10000
    jobs = [(otu_combos[x:x+chunksize], sid_cat, ast_rel_abd,)
            for x in xrange(0, len(otu_combos), chunksize)]
    print("{0} jobs created.".format(len(jobs)))

    # Start multiprocessing jobs
    try:
        print("Starting map_async()...")
        pool = Pool()
        res = pool.map_async(calc_corr_helper, jobs)
        pool.close()
        pool.join()
    except Exception:
        sys.exit("Error while calculating correlations\n{}".format(format_exc()))
    else:
        s_rho_calc = []
        k_tau_calc = []
        for r in res.get():
            for s in r:
                if s[0] == "Spearman":
                    s_rho_calc.append(s)
                else:
                    k_tau_calc.append(s)

    # Get FDR corrected correlation results
    print("Running FDR correction on {} Spearman's Rho.".format(len(s_rho_calc)))
    fdr_corr_s_rho = run_fdr(s_rho_calc)
    print("Running FDR correction on {} Kendall Tau.".format(len(k_tau_calc)))
    fdr_corr_k_tau = run_fdr(k_tau_calc)

    # Consolidate correlation results
    k_kos = {(e[2], e[3],) for e in fdr_corr_k_tau}
    s_kos = {(f[2], f[3],) for f in fdr_corr_s_rho}
    final_kos = s_kos & k_kos
    print("{0} elements from KendallTau\n{1} elements from SpearmanRho\n{2} elements are "
          "common to both.".format(len(k_kos), len(s_kos), len(final_kos)))
    final_fdr_corr_results = [cdata[1:] for cdata in fdr_corr_s_rho
                              if (cdata[2], cdata[3],) in final_kos]

    # Write our results to file
    with open(args.out_fnh, "w") as outf:
        outf.write("Category\tVariable\tby Variable\tCorrelation\tp value\n")
        for k in final_fdr_corr_results:
            outf.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(k[0], k[1], k[2], k[3], k[4]))


if __name__ == "__main__":
    sys.exit(main())
