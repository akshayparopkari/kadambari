#!/usr/bin/env python
"""
:Abstract: Universality of microbiome: Dissimilarity-Overlap Curve Analysis. Source -
           Bashan, A., Gibson, T. E., Friedman, J., Carey, V. J., Weiss, S. T.,
           Hohmann, E. L., & Liu, Y. Y. (2016). Universality of human microbial dynamics.
           Nature, 534(7606), 259-262.
:Date: 11/08/2016
:Author: Akshay Paropkari
"""
import sys
import argparse
from random import choice
from math import sqrt, log
from operator import itemgetter
from itertools import combinations
from collections import defaultdict
from phylotoast.graph_util import ggplot2_style
from phylotoast.util import gather_categories, parse_map_file
importerrors = []
try:
    import biom
except ImportError:
    importerrors.append("biom")
try:
    import numpy as np
except ImportError:
    importerrors.append("numpy")
try:
    import pandas as pd
except ImportError:
    importerrors.append("pandas")
try:
    import matplotlib as mpl
    mpl.rc("font", family="Arial")
    mpl.rc("xtick", labelsize=12)
    mpl.rc("ytick", labelsize=12)
    import matplotlib.pyplot as plt
except ImportError:
    importerrors.append("matplotlib")
try:
    import seaborn as sns
except ImportError:
    importerrors.append("seaborn")
try:
    from scipy.stats import entropy
except ImportError:
    importerrors.append("scipy")
try:
    from sklearn.metrics import r2_score
except ImportError:
    importerrors.append("scikit-learn")
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
except ImportError:
    importerrors.append("statsmodels")
try:
    assert len(importerrors) == 0
except AssertionError:
    for err in importerrors:
        print("Please install missing module: {}".format(err))
    sys.exit()


def root_jsd(x, y):
    """
    Discrete root Jensen Shannon divergence (rJSD) calculation. Parameters x and y must
    denote abundance for same otuids.

    :type x: int or float or list
    :param x: renormalized relative abundance from Sample 1 for an otuid.
    :type y: int or float or list
    :param y: renormalized relative abundance from Sample 2 for an otuid.
    """
    try:
        assert isinstance(x, list)
        assert isinstance(y, list)
    except:
        # For individual inputs
        m = (x + y) / 2
        try:
            assert m > 0
        except AssertionError:
            return 0
        else:
            kld_x = x * log(x / m)
            kld_y = y * log(y / m)
            rJSD = sqrt(0.5 * kld_x + 0.5 * kld_y)
    else:
        # For list inputs
        m = [(a + b) / 2 for a, b in zip(x, y)]
        try:
            assert m > 0
        except AssertionError:
            return 0
        else:
            kld_x = entropy(x, m)
            kld_y = entropy(y, m)
            rJSD = sqrt(0.5 * kld_x + 0.5 * kld_y)
    return rJSD


def calc_doc(biomfile, samples=None):
    """
    Calculate the Dissimilarity-Overlap values for each input sample pair.

    :type biomfile: biom format
    :param biomfile: BIOM file which has been normalized using biom.norm() along sample
                     axis.
    :type samples: list
    :param samples: Samples for which to calculate the Dissimilarity-Overlap values.
    """
    # Create object to hold calculated values
    pairwise_calc = defaultdict(dict)

    # Get dissimilarity and shared otu
    for pair in combinations(samples, 2):
        pairwise_calc["{}_{}".format(pair[0], pair[1])] = {"Overlap": 0,
                                                           "Dissimilarity": 0}
        renorm_rel_abd = defaultdict(dict)
        shared_species = set()
        for otuid in biomfile.ids("observation"):
            # Only if present in both samples
            if biomfile.get_value_by_ids(otuid, pair[0]) > 0 and \
               biomfile.get_value_by_ids(otuid, pair[1]) > 0:

                # Dissimilarity
                count_sum = sum([biomfile.get_value_by_ids(otuid, pair[0]),
                                 biomfile.get_value_by_ids(otuid, pair[1])])
                abd1 = biomfile.get_value_by_ids(otuid, pair[0]) / count_sum
                abd2 = biomfile.get_value_by_ids(otuid, pair[1]) / count_sum
                renorm_rel_abd[pair[0]][otuid] = abd1
                renorm_rel_abd[pair[1]][otuid] = abd2
                assert renorm_rel_abd[pair[0]].keys() == renorm_rel_abd[pair[1]].keys()

                # Overlap
                shared_species.add(otuid)

        # Calculate the shared species overlap
        for shared_otuid in shared_species:
            abd1 = biomfile.get_value_by_ids(shared_otuid, pair[0])
            abd2 = biomfile.get_value_by_ids(shared_otuid, pair[1])
            shared_otuid_avg = np.mean([abd1, abd2])
            assert shared_otuid_avg == sum([abd1, abd2]) / 2
            pairwise_calc["{}_{}".format(pair[0], pair[1])]["Overlap"] += shared_otuid_avg

        # Dissimilarity calculation
        pairwise_calc["{}_{}".format(pair[0], pair[1])]["Dissimilarity"] +=\
            root_jsd(renorm_rel_abd[pair[0]].values(), renorm_rel_abd[pair[1]].values())

    # Final check for each sample pair calculations performance
    try:
        assert len(pairwise_calc) == (len(samples)*(len(samples) - 1)) / 2
    except AssertionError:
        print("Dissimilarity-overlap values were not calculated for all sample pairs.")
        return None
    return pairwise_calc


def get_doc_ci(do_values, frac, ci=False, samples=None, num_of_seqs=100):
    """
    Calculate DOC and the confidence interval via bootstrapping.

    :type do_values: dict
    :param do_values: The output from calc_doc() containing dissimilarity and overlap
                      values for each pair of input samples.
    :type frac: float
    :param frac: The fraction of data points used for predicting y-value.
    :type samples: list
    :param samples: Samples for which the Dissimilarity-Overlap values were calculated.
    :type num_of_seqs: int
    :param num_of_seqs: Number of sample lists to generate.
    """
    do_values_df = pd.DataFrame.from_dict(do_values, orient="index")

    try:
        assert ci == False
    except:
        # Get a random list of input samples with replacement
        sids_list = [[choice(samples) for _ in range(len(samples))]
                     for _ in range(num_of_seqs)]

        # Get sorted combination of random sample list
        sids_combo = [{"{}_{}".format(pair[0], pair[1])
                       for pair in combinations(sid_seq, 2)} for sid_seq in sids_list]

        # Get indices of all combinations present in do_values_df
        idx = [[combo for combo in sid_seq if combo in sorted(do_values.keys())]
               for sid_seq in sids_combo]

        lowess_regression = defaultdict(list)
        for row in idx:
            subset_df = do_values_df.ix[row]
            subset_df.sort_values("Overlap", inplace=True)
            # Apply smooth lowess regression
            lr = lowess(np.asarray(subset_df["Dissimilarity"]),
                        np.asarray(subset_df["Overlap"]), frac=frac, is_sorted=True,
                        return_sorted=False)
            for sample_pair, d in zip(subset_df.index.tolist(), lr):
                lowess_regression[sample_pair].append(d)
        try:
            assert sorted(do_values_df.index) == sorted(lowess_regression.keys())
        except AssertionError:
            sys.exit("Error: All samples pairs not included in estimating CI. Please "
                     "increase num_of_seqs parameter value.".upper())
        else:
            lowess_regression_min = [np.percentile(lowess_regression[k], 3)
                                     for k in sorted(lowess_regression.keys())]
            lowess_regression_max = [np.percentile(lowess_regression[k], 97)
                                     for k in sorted(lowess_regression.keys())]
        do_values_df["LOWESS_min"] = pd.Series(lowess_regression_min,
                                           index=do_values_df.index.values)
        do_values_df["LOWESS_max"] = pd.Series(lowess_regression_max,
                                               index=do_values_df.index.values)

    # Get overall LOWESS
    do_values_df.sort_values("Overlap", inplace=True)
    overall_lr = lowess(np.asarray(do_values_df["Dissimilarity"]),
                        np.asarray(do_values_df["Overlap"]), frac=frac, is_sorted=True,
                        return_sorted=False)
    do_values_df["LOWESS"] = pd.Series(overall_lr, index=do_values_df.index.values)
    return do_values_df


def plot_doc(do_values_df, order, ci=False, title=None, save=None):
    """
    Plot DOC diagram with confidence intervals.

    :type do_values_df: Pandas Dataframe
    :param do_values_df: Dataframe containing dissimilarity and overlap values with
                         lowess predictions.
    :type order: int
    :param order: Order of the polynomial to fit when calculating the residuals.
    :type ci: boolean
    :param ci: Set to plot the confidance intervals. Default is no plotting CI.
    :type save: boolean
    :param save: Text for plot title.
    :type save: boolean
    :param save: Set to save the plot.
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    ax.scatter(x=do_values_df["Overlap"], y=do_values_df["Dissimilarity"], s=25,
               c="#111111", label="Dissimilarity-Overlap")
    ax.plot(do_values_df["Overlap"], do_values_df["LOWESS"], "r-", lw=3,
            label="LOWESS Smoothing", alpha=0.75)
    # Plot LOWESS fit R^2
    Rsqr_ld = r2_score(y_true=do_values_df["Dissimilarity"],
                       y_pred=do_values_df["LOWESS"],
                       multioutput="uniform_average")
    ld_r2_text = r"LOWESS Fit $R^2$: {:.3f}".format(Rsqr_ld)
    ax.text(x=0.1, y=0.97, s=ld_r2_text, ha="center", va="center", fontsize=14,
            transform=ax.transAxes)
    try:
        assert ci == False
    except AssertionError:
        ymin = []
        ymax = []
        for rows in do_values_df.iterrows():
            row = rows[1]
            min_calc = abs(row["LOWESS"] - row["LOWESS_min"])
            max_calc = abs(row["LOWESS_max"] - row["LOWESS"])
            ymin.append(row["LOWESS"] - min_calc)
            ymax.append(row["LOWESS"] + max_calc)
        ax.fill_between(do_values_df["Overlap"], ymin, ymax, color="#ff0000", alpha=0.25)

    # Plot polynomial fit
#     z_orig = np.polyfit(x=do_values_df["Overlap"], y=do_values_df["Dissimilarity"],
#                         deg=order)
#     p_orig = np.poly1d(z_orig)
#     ax.plot(do_values_df["Overlap"], p_orig(do_values_df["Overlap"]), "b--", lw=3,
#             label="Polynomial Fit (degree {})".format(order), alpha=0.75)
#     polyfit_r2 = r2_score(y_true=do_values_df["Dissimilarity"],
#                           y_pred=p_orig(do_values_df["Overlap"]),
#                           multioutput="uniform_average")
#     poly_r2_text = r"Polynomial Fit $R^2$: {:.3f}".format(polyfit_r2)
#     ax.text(x=0.1, y=0.92, s=poly_r2_text, ha="center", va="center", fontsize=14,
#             transform=ax.transAxes)

    ax.legend(fontsize=14, scatterpoints=4)
    ax.set_xlabel("Overlap", fontsize=14)
    ax.set_ylabel("Dissimilarity", fontsize=14)
    ax.set_xlim([0, 1.0])
    ax.set_ylim([0, 1.0])
    if title is not None:
        ax.set_title(title, weight="bold", fontsize=15)
    ggplot2_style(ax)
    fig.tight_layout()
    if save is not None:
        fig.savefig(save, facecolor="white", edgecolor="none",
                    bbox_inches="tight", pad_inches=0.2)
    else:
        plt.show()


def plot_residplot(do_values_df, order, save=None):
    """
    Plot residual plot, which can help in determining if there is structure to the
    residuals.

    :type do_values_df: Pandas Dataframe
    :param do_values_df: Dataframe containing dissimilarity and overlap values with
                         lowess predictions.
    :type order: int
    :param order: Order of the polynomial to fit when calculating the residuals.
    :type save: boolean
    :param save: Set to save the plot.
    """
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    sns.residplot(x="Overlap", y="Dissimilarity", data=do_values_df, lowess=True,
                  order=order, ax=ax)
    ax.set_title("Residual Plot", weight="bold", fontsize=15)
    if save is not None:
        fig.savefig(save, facecolor="white", edgecolor="none", format="svg",
                    bbox_inches="tight", pad_inches=0.1)
    else:
        plt.show()


def handle_program_options():
    parser = argparse.ArgumentParser(description="Dissimilarity-Overlap Curve Analysis "
                                     "as published in http://www.nature.com/nature/journa"
                                     "l/v534/n7606/full/nature18301.html. LOWESS "
                                     "smoothing is achieved through statsmodels package.")
    parser.add_argument("input_biom_fp", help="BIOM file path. [REQUIRED]")
    parser.add_argument("map_fp", help="Input metadata mapping filepath [REQUIRED].")
    parser.add_argument("-b", "--group_by", default=None,
                        help="Any mapping categories, such as treatment type, that will "
                             "be used to group the data in the output iTol table. For "
                             "example, one category with three types will result in "
                             "three data columns in the final output. Two categories with"
                             "three types each will result in six data columns. Default "
                             "is no categories and all the data will be treated as a "
                             "single group.")
    parser.add_argument("-rp", "--residplot", default=2, type=int,
                        help="Order of the polynomial to fit when calculating the"
                             "residuals. Default is set to 2nd degree polynomial fit.")
    parser.add_argument("-f", "--frac", type=float, default=0.5,
                        help="Between 0 and 1. The fraction of the data used when "
                             "estimating each y-value.")
    parser.add_argument("-pc", "--plot_ci", action="store_true",
                        help="Boolean to specify whether to plot the confidance intervals"
                             " or not. By default, CI will not be plotted.")
    parser.add_argument("-sc", "--save_calc",
                        help="Provide a path and file name to save the dissimilarity and "
                             "overlap calculations used for plotting. By default, "
                             "calculations will not be saved.")
    parser.add_argument("-n", "--num_iterations", default=100,
                        help="Number of iterations to perform estimating confidance "
                             "intervals of non-parametric LOWESS. Default is set to 100.")
    parser.add_argument("-t", "--title", help="Title for the output figure.")
    parser.add_argument("-s", "--save_image",
                        help="Path and file name for saving the DOC curve plot. By "
                             "default, the plot will be saved in SVG format.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Read in biom file
    try:
        shared_biom = biom.load_table(args.input_biom_fp)
    except IOError as ie:
        sys.exit("\nError reading BIOM file: {}\n".format(ie))
    norm_shared_biom = shared_biom.norm(axis="sample", inplace=False)

    # Read in mapping file
    try:
        header, imap = parse_map_file(args.map_fp)
    except IOError as ioe:
        sys.exit("\nError in metadata mapping filepath: {}\n".format(ioe))

    # Samples for each group and  get DO values per category
    try:
        assert args.group_by is None
    except AssertionError:
        data_gather = gather_categories(imap, header, args.group_by.split(","))
        sample_list = [sid for cat in data_gather.keys() for sid in data_gather[cat].sids]
    else:
        sample_list = norm_shared_biom.ids()
    doc = calc_doc(norm_shared_biom, sample_list)
    try:
        assert doc is not None
    except AssertionError:
        sys.exit("Error in DOC calculations. Please check the modules.")

    # Get confidence interval
    sl_lowess_regr = get_doc_ci(doc, args.frac, args.plot_ci, sample_list,
                                num_of_seqs=args.num_iterations)

    # Plot the residual figure
    plot_residplot(sl_lowess_regr, args.residplot, save=args.save_image)

    # Plot DOC
    plot_doc(sl_lowess_regr, args.residplot, ci=args.plot_ci, title=args.title,
             save=args.save_image)

if __name__ == "__main__":
    sys.exit(main())
