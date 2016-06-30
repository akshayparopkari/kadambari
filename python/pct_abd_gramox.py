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
except ImportError as ie:
    importerrors.append(ie)
try:
    import numpy as np
except ImportError as ie:
    importerrors.append(ie)
try:
    import pandas as pd
except ImportError as ie:
    importerrors.append(ie)
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except ImportError as ie:
    importerrors.append(ie)
try:
    import seaborn as sns
except ImportError as ie:
    importerrors.append(ie)
try:
    from phylotoast import otu_calc as oc, biom_calc as bc
except ImportError as ie:
    importerrors.append(ie)
if len(importerrors) > 0:
    for err in importerrors:
        print "Import Error: {}".format(err)
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
                                     " data and plot it. Script will also have"
                                     " the option to write out data files (see"
                                     " cat_gramox_out, sample_gramox_out and "
                                     "-c options).")
    parser.add_argument("gramox_fnh",
                        help="File with gramox data for each OTU and 'NA' if "
                             " data is unavailable or inadequate.")
    parser.add_argument("biom_file",
                        help="BIOM file/OTU table.")
    parser.add_argument("mapping_file",
                        help="Mapping file.")
    parser.add_argument("cat_gramox_out",
                        help="File to write out percent abundance data for "
                             "each category.")
    parser.add_argument("sample_gramox_out",
                        help="File to write out percent abundance data for "
                             "each sample.")
    parser.add_argument("category", type=str,
                        help="Column name for sample categories.")
    parser.add_argument("-a", "--annotate_bars", action="store_true",
                        help="If provided, each bar will be annotated with "
                             "percent abundance (Y-axis value) Default is no "
                             "annotation.")
    parser.add_argument("-c", "--per_sid_abd",
                        help="If provided, intermediate calculation values for"
                             " per sample percent abundance will be written to"
                             " this file.")
    parser.add_argument("-o", "--per_sid_otus",
                        help="If provided, per sample otus will be written to"
                             " this file.")
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
    cond = defaultdict()
    for rows in mapf.iterrows():
        row = rows[1]
        cond[row["#SampleID"]] = row[args.category]

    # Get relative abundances from biom file
    biomf = biom.load_table(args.biom_file)
    rel_abd = calc_rel_abd(biomf, cond.keys())
    if args.per_sid_otus:
        with open(args.per_sid_otus, "w") as fdafb:
            for sid, bacteria in rel_abd.iteritems():
                naam = [name for name, c in bacteria.iteritems() if c > 0]
                fdafb.write("{}\t{}\n".format(sid, "\t".join(naam)))
    otus = set()
    for k, v in rel_abd.iteritems():
        for o in v.keys():
            otus.add(o)
    print "Total OTUs:", len(otus)
    for a in sorted(otus):
        if a not in gramox_data.keys():
            raise RuntimeError("OTU list from biom file doesn't match with "
                               "gramox data file list.")

    # Initialize dict to collect gramox abundance data by category and sample
    cat = defaultdict(dict)
    sid_gramox_abd = defaultdict(dict)
    for k in cond.keys():
        cat[cond[k]]["Unknown"] = 0
        cat[cond[k]]["Gram Positive Facultatives"] = 0
        cat[cond[k]]["Gram Positive Anaerobe"] = 0
        cat[cond[k]]["Gram Negative Facultatives"] = 0
        cat[cond[k]]["Gram Negative Anaerobe"] = 0
        cat[cond[k]]["Total"] = 0
        sid_gramox_abd[k]["Unknown"] = 0
        sid_gramox_abd[k]["Gram Positive Facultatives"] = 0
        sid_gramox_abd[k]["Gram Positive Anaerobe"] = 0
        sid_gramox_abd[k]["Gram Negative Facultatives"] = 0
        sid_gramox_abd[k]["Gram Negative Anaerobe"] = 0
        sid_gramox_abd[k]["Total"] = 0

    for sid, d in rel_abd.iteritems():
        for otu, abd in d.iteritems():
            cat[cond[sid]]["Total"] += abd
            sid_gramox_abd[sid]["Total"] += abd
            if gramox_data[otu] == ["1", "1"]:
                cat[cond[sid]]["Gram Positive Facultatives"] += abd
                sid_gramox_abd[sid]["Gram Positive Facultatives"] += abd
            elif gramox_data[otu] == ["1", "0"]:
                cat[cond[sid]]["Gram Positive Anaerobe"] += abd
                sid_gramox_abd[sid]["Gram Positive Anaerobe"] += abd
            elif gramox_data[otu] == ["0", "1"]:
                cat[cond[sid]]["Gram Negative Facultatives"] += abd
                sid_gramox_abd[sid]["Gram Negative Facultatives"] += abd
            elif gramox_data[otu] == ["0", "0"]:
                cat[cond[sid]]["Gram Negative Anaerobe"] += abd
                sid_gramox_abd[sid]["Gram Negative Anaerobe"] += abd
            elif "NA" in gramox_data[otu]:
                cat[cond[sid]]["Unknown"] += abd
                sid_gramox_abd[sid]["Unknown"] += abd

    # Calculate and collect percent abundance data.
    pct_cat = defaultdict(dict)
    pct_sid = defaultdict(dict)
    for condition, v in cat.iteritems():
        for gramox, abd in v.iteritems():
            tot_abd = cat[condition][gramox]
            total = cat[condition]["Total"]
            if gramox != "Total":
                pct_cat[condition][gramox] = tot_abd/total*100

    for sampleid, v1 in sid_gramox_abd.iteritems():
        for gramox1, abd in v1.iteritems():
            tot_abd1 = sid_gramox_abd[sampleid][gramox1]
            total1 = sid_gramox_abd[sampleid]["Total"]
            if sid_gramox_abd != "Total":
                pct_sid[sampleid][gramox1] = tot_abd1/total1*100

    # Write out the categorized gramox data
    with open(args.cat_gramox_out, "w") as kyu:
        kyu.write("condition\tgramox\tabd\ttotal_abd\tpct_abd\n")
        for cond, val in cat.iteritems():
            for gramox in sorted(val.keys()):
                if gramox != "Total":
                    kyu.write("{}\t{}\t{}\t{}\t{}\n".format(cond, gramox,
                              cat[cond][gramox], cat[cond]["Total"],
                              pct_cat[cond][gramox]))

    # Write out gramox data for each sample
    with open(args.sample_gramox_out, "w") as kyu:
        kyu.write("Sample\tGram Positive Facultatives\tGram Positive Anaerobe"
                  "\tGram Negative Facultatives\tGram Negative Anaerobe\t"
                  "Unknown\n")
        for sid, val in sid_gramox_abd.iteritems():
            kyu.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sid,
                      pct_sid[sid]["Gram Positive Facultatives"],
                      pct_sid[sid]["Gram Positive Anaerobe"],
                      pct_sid[sid]["Gram Negative Facultatives"],
                      pct_sid[sid]["Gram Negative Anaerobe"],
                      pct_sid[sid]["Unknown"],))
    if args.per_sid_abd:
        with open(args.per_sid_abd, "w") as qwe:
            qwe.write("sample\tgramox\tabd\ttotal_abd\tpct_abd\n")
            for sid, val in sid_gramox_abd.iteritems():
                for gramox in sorted(val.keys()):
                    if gramox != "Total":
                        qwe.write("{}\t{}\t{}\t{}\t{}\n".format(sid,
                                  gramox,
                                  sid_gramox_abd[sid][gramox],
                                  sid_gramox_abd[sid]["Total"],
                                  pct_sid[sid][gramox]))

    # Plot data
    plot_data = pd.read_csv(args.cat_gramox_out, sep="\t")
    plot_data = plot_data.sort_values("condition")
    palette_col = ["#253494", "#2c7fb8", "#41b6c4", "#a1dab4", "#ffffcc"]
    sns.set(style="whitegrid")
    fig = plt.figure(figsize=(20, 15))
    ax = sns.barplot(
        x="condition", y="pct_abd", hue="gramox", data=plot_data,
        hue_order=sorted(np.unique(plot_data.gramox)), palette=palette_col,
        saturation=1)
    if args.annotate_bars:
        for p in ax.patches:
            ax.annotate(np.round(p.get_height(), decimals=2),
                        (p.get_x()+p.get_width()/2., p.get_height()), ha='center',
                        va='bottom', xytext=(0, 5), textcoords='offset points',
                        size=80*p.get_width())
    ax.set_xlabel("")
    mpl.rc("font", family="Arial")  # define font for figure text
    ax.set_xticklabels(np.unique(plot_data.condition), size=18, alpha=1)
    ax.set_ylabel("Abundance Percentage (%)", labelpad=20, size=18, alpha=1)
    ax.set_ylim(top=int(max(plot_data.pct_abd) + 5))
    leg = ax.legend(loc=0, fontsize=18, frameon=True)
    leg.get_frame().set_linewidth(1)
    leg.get_frame().set_edgecolor('k')
    ax.grid(b=True, which='major', color='k', linestyle=':', alpha=0.5)
    if args.save_fig:
        fig.savefig(args.save_fig, dpi=300, format="svg", bbox_inches="tight",
                    pad_inches=0.5)
    else:
        plt.show()

if __name__ == "__main__":
    main()
