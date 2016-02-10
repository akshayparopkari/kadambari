#!/usr/bin/env python
"""
Abstract: Calculate categorized gramox data and plot it.

Author: Akshay Paropkari

Date: 02/09/2016
"""

import sys
import argparse
from collections import defaultdict
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
    import matplotlib.pyplot as plt
except ImportError:
    importerrors.append("matplotlib")
try:
    import seaborn as sns
except ImportError:
    importerrors.append("seaborn")
try:
    from phylotoast import otu_calc as oc, biom_calc as bc
except ImportError:
    importerrors.append("phylotoast")
if len(importerrors) > 0:
    for err in importerrors:
        print "Please install missing module: {}".format(err)
    sys.exit()


def calc_rel_abd(biomf, sampleIDs=None):
    """
    Calculate relative abundance from a biom table either based on sampleIDs or
    otuIDs.

    :type biomf: Biom file format
    :param biomf: Biom table loaded object

    :type sampleIDs: list
    :param sampleIDs: Only calculate relative abundances for these sampleIDs.
                      Default is None.

    :return type: defaultdict(dict)
    :return: A dict keyed on sampleID with its value denoting a dict keyed on
             otuID and abundance value for that [sampleID, otuID] pair.
    """
    norm_biomf = biomf.norm(inplace=False)
    if sampleIDs is None:
        sampleIDs = norm_biomf.ids()
    otuIDs = norm_biomf.ids(axis="observation")
    rel_abd = defaultdict(dict)

    for sample in sampleIDs:
        for otu in otuIDs:
            abd = norm_biomf.get_value_by_ids(otu, sample)
            otu_tax = norm_biomf.metadata(otu, "observation")["taxonomy"]
            otu_name = oc.otu_name(otu_tax)
            rel_abd[sample][otu_name] = abd
    trans_rel_abd = bc.arcsine_sqrt_transform(rel_abd)
    return trans_rel_abd


def handle_program_options():
    parser = argparse.ArgumentParser(description="Calculate categorized gramox"
                                     " data and plot it.")
    parser.add_argument("gramox_fnh",
                        help="File with gramox data for each OTU and 'NA' if "
                             " data is unavailable or inadequate.")
    parser.add_argument("biom_file",
                        help="BIOM file/OTU table.")
    parser.add_argument("gramox_out",
                        help="File to write out gramox data for each category.")
    parser.add_argument("mapping_file",
                        help="Mapping file.")
    parser.add_argument("category", type=str,
                        help="Column name for sample categories.")
    parser.add_argument("-c", "--colors", type=str, default=None,
                        help="Column name for category colors. Default is husl"
                             " color palette from seaborn library.")
    parser.add_argument("-s", "--save_fig",
                        help="Gramox data plot will be saved into this file.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Read gramox_calc data
    gramox_data = defaultdict()
    with open(args.gramox_fnh, "rU") as rht:
        for line in rht.readlines():
            line = line.strip().split("\t")
            try:
                gramox_data[line[0]] = [line[1], line[2]]
            except IndexError:
                raise IndexError("{} gramox data incomplete".format(line[0]))

    # Read mapping file for category knowledge
    mapf = pd.read_csv(args.mapping_file, sep="\t")
    if args.colors is not None:
        cat_colors = defaultdict()
    cond = defaultdict()
    for rows in mapf.iterrows():
        row = rows[1]
        cond[row["#SampleID"]] = row[args.category]
        if args.colors is not None:
            cat_colors[row[args.category]] = row[args.colors]

    # Get relative abundances from biom file
    biomf = biom.load_table(args.biom_file)
    rel_abd = calc_rel_abd(biomf, cond.keys())
    otus = set()
    for k, v in rel_abd.iteritems():
        for o in v.keys():
            otus.add(o)
    print "Total OTUS:", len(otus)
    for a in sorted(otus):
        if a not in gramox_data.keys():
            raise RuntimeError("OTU list from biom file doesn't match with "
                               "gramox data file list.")

    # Initialize dict to collect gramox abundance data
    cat = defaultdict(dict)
    for k in cond.keys():
        cat[cond[k]]["Unknown"] = 0
        cat[cond[k]]["Gram Positive Aerobe"] = 0
        cat[cond[k]]["Gram Positive Anaerobe"] = 0
        cat[cond[k]]["Gram Negative Aerobe"] = 0
        cat[cond[k]]["Gram Negative Anaerobe"] = 0
        cat[cond[k]]["Total"] = 0

    for sid, d in rel_abd.iteritems():
        for otu, abd in d.iteritems():
            cat[cond[sid]]["Total"] += abd
            if gramox_data[otu] == ["1", "1"]:
                cat[cond[sid]]["Gram Positive Aerobe"] += abd
            elif gramox_data[otu] == ["1", "0"]:
                cat[cond[sid]]["Gram Positive Anaerobe"] += abd
            elif gramox_data[otu] == ["0", "1"]:
                cat[cond[sid]]["Gram Negative Aerobe"] += abd
            elif gramox_data[otu] == ["0", "0"]:
                cat[cond[sid]]["Gram Negative Anaerobe"] += abd
            elif "NA" in gramox_data[otu]:
                cat[cond[sid]]["Unknown"] += abd

    # Calculate and collect percent abundance data.
    pct_cat = defaultdict(dict)
    for condition, v in cat.iteritems():
        for gramox, abd in v.iteritems():
            tot_abd = cat[condition][gramox]
            total = cat[condition]["Total"]
            if gramox != "Total":
                pct_cat[condition][gramox] = tot_abd/total*100

    # Write out the categorized gramox data
    with open(args.gramox_out, "w") as kyu:
        kyu.write("condition\tgramox\tabd\ttotal_abd\tpct_abd\n")
        for cond, val in cat.iteritems():
            for gramox in sorted(val.keys()):
                if gramox != "Total":
                    kyu.write("{}\t{}\t{}\t{}\t{}\n".format(cond, gramox,
                              cat[cond][gramox], cat[cond]["Total"],
                              pct_cat[cond][gramox]))

    # Plot data
    plot_data = pd.read_csv(args.gramox_out, sep="\t")
    cat_order = sorted(np.unique(plot_data.condition))
    if args.colors is not None:
        palette_col = [cat_colors[c] for c in cat_order]
    else:
        palette_col = sns.set_palette("husl")
    fig = plt.figure(figsize=(20, 15))
    ax = sns.barplot(
        x="gramox", y="pct_abd", hue="condition", data=plot_data,
        hue_order=cat_order, palette=palette_col, saturation=1)
    for p in ax.patches:
        ax.annotate(np.round(p.get_height(), decimals=2),
                    (p.get_x()+p.get_width()/2., p.get_height()), ha='center',
                    va='bottom', xytext=(0, 5), textcoords='offset points',
                    size=50*p.get_width())
    ax.set_xlabel("")
    ax.set_xticklabels(np.unique(plot_data.gramox), size=18, alpha=1)
    ax.set_ylabel("Abundance Percentage (%)", labelpad=20, size=18, alpha=1)
    ax.set_yticklabels(range(0, int(max(plot_data.pct_abd) + 10), 5),
                       alpha=1, size=18)
    ax.legend(loc=0, fontsize=18, frameon=True)
    ax.grid()
    if args.save_fig:
        fig.savefig(args.save_fig, dpi=300, format="svg", bbox_inches="tight",
                    pad_inches=0.2, )
    plt.show()

if __name__ == "__main__":
    main()
