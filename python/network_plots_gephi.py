#!/usr/bin/env python
"""
Abstract: Generate network plots using NetworkX and Gephi.
Author: Akshay Paropkari
Date: 02/10/2016
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
    import pandas as pd
except ImportError as ie:
    importerrors.append(ie)
try:
    import networkx as nx
except ImportError:
    importerrors.append(ie)
try:
    from phylotoast import biom_calc as bc, otu_calc as oc
except ImportError:
    importerrors.append(ie)
if len(importerrors) > 0:
    for err in importerrors:
        print "Import Error: {}".format(err)
    sys.exit()


def get_relative_abundance(biomfile):
    """Return arcsine transformed relative abundance from a OTU table."""
    biomf = biom.load_table(biomfile)
    norm_biomf = biomf.norm(inplace=False)
    rel_abd = {}
    for sid in norm_biomf.ids():
        rel_abd[sid] = {}
        for otuid in norm_biomf.ids("observation"):
            otuname = oc.otu_name(norm_biomf.metadata(otuid, axis="observation")["taxonomy"])
            otuname = " ".join(otuname.split("_"))
            abd = norm_biomf.get_value_by_ids(otuid, sid)
            rel_abd[sid][otuname] = abd
    ast_rel_abd = bc.arcsine_sqrt_transform(rel_abd)
    return ast_rel_abd


def handle_program_options():
    """Parses the given options passed in at the command line."""
    parser = argparse.ArgumentParser(description="Create network plots based "
                                     "on correlation matrix.")
    parser.add_argument("biom_file", help="Biom file OTU table.")
    parser.add_argument("mapping_file", help="Mapping file for reading "
                        "sampleIDs and their groups.")
    parser.add_argument("condition_column", help="Column name in mapping file "
                        "denoting the categories.")
    parser.add_argument("in_corr_mat", help="Correlation matrix file. The "
                        "format for the tab-separated file should be: "
                        "Category -> Variable -> by Variable -> Correlation")
    parser.add_argument("cat_name", help="Category to be plotted.")
    parser.add_argument("-go", "--gexf_out",
                        help="Graph information written to this Graph Exchange"
                        " XML Format file. This file can be input to Gephi.")
    parser.add_argument("-fp", "--fil_pct", type=float, default=0.75,
                        help="Specify the minimum value of correlation "
                        "strength to display. By default, all correlations "
                        ">=0.75 will be shown.")
    parser.add_argument("-w", "--stats_out_fnh",
                        help="Write out graph statistics.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Error check input file
    try:
        with open(args.biom_file):
            pass
    except IOError as ioe:
        sys.exit("\nError in BIOM file path: {}\n".format(ioe))
    try:
        with open(args.mapping_file):
            pass
    except IOError as ioe:
        sys.exit("\nError in mapping file path: {}\n".format(ioe))
    try:
        with open(args.in_corr_mat):
            pass
    except IOError as ioe:
        sys.exit("\nError in correlation matrix file path: {}\n".format(ioe))

    relative_abundance = get_relative_abundance(args.biom_file)

    # get categorized nodes
    mapf = pd.read_csv(args.mapping_file, sep="\t")
    cat_sids = defaultdict(list)
    for rows in mapf.iterrows():
        row = rows[1]
        cat_sids[row["Condition"]].append(row["#SampleID"])

    # read the correlation data
    corr_data = pd.read_csv(args.in_corr_mat, sep="\t")

    # initialize graph
    G = nx.Graph()

    # get categorized graph edge data and edge colors
    for rows in corr_data.iterrows():
        row = rows[1]
        if row["Category"] == args.cat_name and abs(row["Correlation"]) > args.fil_pct:
            if row["Correlation"] > 0:
                G.add_weighted_edges_from(
                    [(" ".join(row["Variable"].split("_")),
                      " ".join(row["by Variable"].split("_")),
                      row["Correlation"])], color="#00FF00")
            else:
                G.add_weighted_edges_from(
                    [(" ".join(row["Variable"].split("_")),
                      " ".join(row["by Variable"].split("_")),
                      (row["Correlation"])*-1)], color="#FF0000")

    # add node attributes to graph node object
    for otu, attr in G.node.iteritems():
        total_abd = 0
        for sid in cat_sids[args.cat_name]:
            total_abd += relative_abundance[sid][otu]
        mean_otu_cat_abd = total_abd/len(cat_sids[args.cat_name])
        G.node[otu]["node_size"] = mean_otu_cat_abd

    # convert node labels to integers
    H = nx.convert_node_labels_to_integers(G, first_label=1,
                                           ordering="decreasing degree",
                                           label_attribute="id")

    # Write out GEXF file for Gephi
    if args.gexf_out:
        nx.write_gexf(H, args.gexf_out, version="1.2draft")

    # Write out betweenness centrality measure to file
    if args.stats_out_fnh:
        with open(args.stats_out_fnh, "w")as poi:
            poi.write("OTU Node\tDegree Centrality\tBetweenness Centrality\n")
            dc = nx.degree_centrality(G)
            bc = nx.betweenness_centrality(G, weight="weight")
            for key in sorted(bc.keys()):
                poi.write("{}\t{}\t{}\n".format(key, dc[key], bc[key]))

if __name__ == "__main__":
    sys.exit(main())
