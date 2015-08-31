#!/usr/bin/env python
'''
Abstract: Create network plots based on correlation matrix.

Date: 04/27/2015

Author: Akshay Paropkari
'''

import sys
import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict


def nodes_classify(gramox_inf, nodelist, lvl):
    '''
    Classify nodes of the network graph based on gram strain of the OTU's.

    :type gramox_inf: gramox data file
    :param gramox_inf: Master file of all gramox data for all OTU's in the
                       following tab-separated format:
                       OTU    Gram Status    Oxygen Requirement    Source

    :type nodelist: list
    :param nodelist: List of all nodes in the graph. Usually networkx
                     attribute such as G.nodes() would work.

    :type lvl: str
    :param lvl: Choose between genus(g) or species(s) phylogenetic level to
                use for classifying OTU's. Defaults to species(s) level.

    :type return: list
    :return: Returns 3 lists g_pos, g_neg, and unk, which have classified
             gram-positive, gram-negative and unknown/NA OTU's respectively.
    '''
    g_pos = []
    g_neg = []
    unk = []
    ctk = []
    with open(gramox_inf, 'rU') as gramoxf:
        if lvl == 'g':
            data = {line.strip().split('\t')[0]: line.strip().split('\t')[2]
                     for line in gramoxf.readlines()[1:]}
        elif lvl == 's':
            data = {line.strip().split('\t')[1]: line.strip().split('\t')[2]
                     for line in gramoxf.readlines()[1:]}
        for item in nodelist:
            if item in data.keys():
                if data[item] == '1':
                    g_pos.append(item)
                elif data[item] == '0':
                    g_neg.append(item)
                else:
                    unk.append(item)
            elif item[:2] == 'Hu':
                ctk.append(item)
            else:
                unk.append(item)
    return g_pos, g_neg, ctk, unk

def draw_network_graphs(pearson_inf, gramox_inf, lvl, filter_pct=None,
                        category=None):
    '''
    This function accepts statistically significant Pearson's correlation
    data to create graph nodes and edges and draws them out on a network
    graph.

    :type pearsoncorr: file path
    :param pearsoncorr: File with output of JMP correlation. The format
                        (first row) for the tab-separated file should be:
                        Category->Variable->by Variable->Correlation

    :type gramox_inf: gramox data file
    :param gramox_inf: Master file of all gramox data for all OTU's in the
                       following tab-separated format:
                       Genus->Species->Gram Status->Oxygen Requirement->Source

    :type lvl: str
    :param lvl: Choose between genus(g) or species(s) phylogenetic level to
                use for classifying OTU's. Defaults to species(s) level.

    :type filter_pct: float
    :param filter_pct: Specify the minimum value of correlation strength
                       to display. By default, all correlations will be
                       portrayed. Range is (0,1).

    :type category: str
    :param category: Provide for which category you want to create a
                     network graph, which should be one of the options
                     from the first column of pearsoncorr data file.

    :type return: network graph/figure
    :return: Returns a network graph with OTU or cytokines as nodes and
             their Pearson's correlation as edges. Green edges
             represent positive correlation and red edges denote negative
             correlation. Also, OTU nodes are colored based on gram strain,
             dark blue are gram-positive and light blue are gram-negative, '
             'yellow are cytokines.'
    '''
    # Read pearson correlation data into a dataframe.
    pdata = pd.read_csv(pearson_inf, sep='\t')

    # Creating a multigraph
    G = nx.MultiGraph()

    # Prep edges for graph
    pdata = pd.read_csv(pearson_inf, sep='\t')
    for rows in pdata.iterrows():
        row = rows[1]
        if category is None:
            G.add_edge(row['Variable'], row['by Variable'], weight=row['Correlation'])
        else:
            if row['Category'] == category:
                if filter_pct is None:
                    G.add_edge(row['Variable'], row['by Variable'], weight=row['Correlation'])
                else:
                    if row['Correlation'] >= filter_pct or row['Correlation'] <= -(filter_pct):
                        G.add_edge(row['Variable'], row['by Variable'], weight=row['Correlation'])

    print 'Length of nodes and edges:', len(G.nodes()), len(G.edges())

    # Classify positive or negative correlation edges
    pos_corr = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0]
    neg_corr = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] < 0]

    # Classify nodes based on gram strains
    g_pos, g_neg, ctk, unk = nodes_classify(gramox_inf, G.nodes(), lvl)
    print 'Length of gpos, gneg, ctk, unk:', len(g_pos), len(g_neg), len(ctk), len(unk)

    # Confirm if all nodes have been classified
    try:
        assert len(g_neg) + len(g_pos) + len(ctk) + len(unk) == len(G.nodes())
    except AssertionError:
        return 'Classified nodes do not add up to total number of nodes.'

    # Draw network graph
    plt.figure(figsize=(20, 20))
    pos = nx.spring_layout(G, iterations=200, k=0.2)
    nx.draw_networkx_nodes(G, pos, nodelist=g_pos, node_color='#3366ff',
                           node_size=500) # gpos: dark blue
    nx.draw_networkx_nodes(G, pos, nodelist=g_neg, node_color='#99ccff',
                           node_size=500) # gneg: light blue
    nx.draw_networkx_nodes(G, pos, nodelist=ctk, node_color='#FFDB19',
                           node_size=500)   # cytokines: yellow
    nx.draw_networkx_nodes(G, pos, nodelist=unk, node_color='#808080',
                           node_size=500)   # unknown: gray
    nx.draw_networkx_edges(G, pos, alpha=0.5, edgelist=pos_corr,
                           edge_color='#008000')  # poscorr: dark green
    nx.draw_networkx_edges(G, pos, alpha=0.5, edgelist=neg_corr,
                           edge_color='r')        # negcorr: red

    nx.draw_networkx_labels(G, pos)
    font = {'color': 'k', 'fontweight': 'bold', 'fontsize': 24}
    plt.axis('off')
    plt.show()


def prog_options():
    parser = argparse.ArgumentParser(description='Create network plots based '
                                    'on correlation matrix.')
    parser.add_argument('in_corr_mat',
                        help='Correlation matrix file. The format'
                             ' for the tab-separated file should be: '
                        'Category->Variable->by Variable->Correlation')
    parser.add_argument('in_gramox_fnh',
                        help='Master file of all gramox data for all OTU\'s '
                             'in the following tab-separated format: '
                       'Genus->Species->Gram Status->Oxygen Requirement->Source')
    parser.add_argument('phy_lvl',
                        help='Choose between genus(g) or species(s) '
                             'phylogenetic level to use for classifying '
                             'OTU\'s. Defaults to species(s) level')
    parser.add_argument('fil_pct', type=float,
                        help='Specify the minimum value of correlation '
                             'strength to display. By default, all '
                             'correlations will be portrayed. Range is (0,1)')
    parser.add_argument('cat_name',
                        help='Program will plot network graph for this '
                             'category only')
    return parser.parse_args()


def main():

    args = prog_options()

    draw_network_graphs(args.in_corr_mat, args.in_gramox_fnh, args.phy_lvl,
                        args.fil_pct, args.cat_name)

if __name__ == '__main__':
    main()
