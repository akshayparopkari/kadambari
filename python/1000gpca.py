#! /usr/bin/env python

"""
:Abstract: Principal Component Analysis on Chr21 of 1000 Genomes Project dataset.
:Author: Akshay Paropkari
:Date: 04/07/2018
"""

import sys
import gzip
import argparse
err = []
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc("font", family="Arial")
    mpl.rc("xtick", labelsize=8)  # set X axis ticksize
    mpl.rc("ytick", labelsize=8)  # set Y axis ticksize
except ImportError:
    err.append("matplotlib")
try:
    import numpy as np
except ImportError:
    err.append("numpy")
try:
    import pandas as pd
except ImportError:
    err.append("pandas")
try:
    from sklearn.decomposition import PCA
except ImportError:
    err.append("scikit-learn")
try:
    from palettable.colorbrewer.sequential import Blues_9    # Europe
    from palettable.colorbrewer.sequential import Greens_9   # East Asia
    from palettable.colorbrewer.sequential import Purples_9  # South Asia
    from palettable.colorbrewer.sequential import Reds_9     # Americas
    from palettable.colorbrewer.sequential import YlOrBr_9   # Sub-Saharan Africa
except ImportError:
    err.append("palettable")
try:
    assert len(err) == 0
except AssertionError:
    for error in err:
        print("Please install {}".format(error))
    sys.exit()


def tf_variants(input_list):
    """Get a list of variant calls and transform them to feature list"""
    feature_map = {"0|0": 0, "1|0": 1, "0|1": 1, "1|1": 2}
    return np.asarray([feature_map[item] if item in feature_map else 3
                       for item in input_list])


def standardize_features(input_list, mean_vec, std_vec):
    """From a list of features, center by mean and scale by standard variance"""
    try:
        assert mean_vec != 0
        assert std_vec != 0
    except AssertionError:
        return np.asarray(input_list)
    else:
        return np.asarray(([(input_list[i] - mean_vec) / std_vec
                            for i in range(len(input_list))]))


def get_valid_data(row_list, indices_to_keep_list):
    """From a row of variant call data, return specific entries for further processing"""
    return np.asarray([row_list[entry] for entry in indices_to_keep_list])


def handle_program_options():
    parser = argparse.ArgumentParser(description="Run PCA on Chr21 of 1000 genomes "
                                     " project.")
    parser.add_argument("-vcf", "--chr21_vcf_file", help="Path to input Chr21 VCF file")
    parser.add_argument("-md", "--map_fp",
                        help="Metadata mapping file corresponding to Chr21 VCF")
    parser.add_argument("-p", "--pca_in", help="Input PCA file.")
    parser.add_argument("-s", "--savefile",
                        help="Save the plot as sn SVG file.")
    parser.add_argument("-o", "--output_file",
                        help="Save consolidated data to this tab-separated file")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Get metadata for each genome
    if args.map_fp:
        with open(args.map_fp, "r") as md:
            md_data = dict()
            for line in md.readlines()[1:]:
                line = line.split()
                md_data[line[0]] = line[1]

    # Discard these population samples
    discarded_pop = tuple(["GIH", "CEU", "PEL", "MXL", "CLM", "PUR", "ASW", "ACB"])
    discarded_genomes = [genome for genome, pop in md_data.items()
                         if pop in discarded_pop]

    # Obtain feature list
    feature_data = []
    mean_row_vector = []
    std_row_vector = []
    feature_data_std = []
    if args.chr21_vcf_file:
        with gzip.open(args.chr21_vcf_file, "rb") as vcff:
            for a, line in enumerate(vcff):
                try:
                    line = line.decode()
                    assert line.startswith("#CHROM")
                except Exception:
                    try:
                        assert line.startswith("21")
                    except Exception:
                        continue
                    else:
                        line = line.strip().split("\t")
                        try:
                            assert line[3] in ["A", "T", "C", "G"]
                            assert line[4] in ["A", "T", "C", "G"]
                            common_var = (2 * (line[9:].count("1|1")) +
                                          line[9:].count("1|0") + line[9:].count("0|1")) \
                                          / 2 * len(line[9:])
                            assert common_var > 0.05
                        except AssertionError:
                            continue
                        else:
                            valid_data = get_valid_data(line[9:], genome_to_keep)
                            tf_data = np.asarray(tf_variants(valid_data))
                            feature_data.append(tf_data)
                            entry_mean = np.mean(tf_data, dtype=np.float64)
                            entry_std = np.std(tf_data, dtype=np.float64)
                            mean_row_vector.append([entry_mean])
                            std_row_vector.append([entry_std])
                            tf_features = np.asarray(standardize_features(tf_data,
                                                                          entry_mean,
                                                                          entry_std))
                            feature_data_std.append(tf_features)
                else:
                    line = line.strip().split("\t")
                    genome_to_keep = [i for i, genome in enumerate(line[9:])
                                      if genome not in discarded_genomes]
                    genome_order = get_valid_data(line[9:], genome_to_keep)

            # Run dimensionality reduction using PCA
            try:
                print("\nPCA starting...")
                pca = PCA(n_components=2)
                pca_ft = pca.fit_transform(np.asarray(feature_data_std).T)
            except Exception as ex:
                sys.exit(ex)
            else:
                if args.output_file:
                    with open(args.output_file, "w") as outf:
                        outf.write("Explained Variance:\t{}".
                                   format(pca.explained_variance_ratio_))
                        outf.write("Singular values:\t{}".format(pca.singular_values_))
                        for i, entry in enumerate(pca_ft):
                            outf.write("{0}\t{1}\t{2}\n".format(genome_order[i], entry[0],
                                                                entry[1]))
                print("Finished!\n")

    # Plot PCA
    if args.pca_in:
        colors = Blues_9.hex_colors[4:] + Greens_9.hex_colors[4:] +\
                 Purples_9.hex_colors[4:] + Reds_9.hex_colors[2:6] +\
                 YlOrBr_9.hex_colors[2:]
        x_order = tuple(["FIN", "GBR", "CEU", "IBS", "TSI", "CHS", "CDX", "CHB",
                         "JPT", "KHV", "GIH", "STU", "PJL", "ITU", "BEB", "PEL",
                         "MXL", "CLM", "PUR", "ASW", "ACB", "GWD", "YRI", "LWK",
                         "ESN", "MSL"])
        assigned_colors = {a: b for a, b in zip(x_order, colors)}
        for entry in discarded_pop:
            del assigned_colors[entry]
        pca_data = pd.read_csv(args.pca_in, sep="\t", header=None, skiprows=1,
                               index_col=False, names=["genome", "PC1", "PC2"])
        with mpl.style.context("ggplot"):
            plt.figure(figsize=(10, 7))
            for entry in pca_data.itertuples():
                try:
                    assert abs(entry[2]/np.std(pca_data["PC1"].values)) <= 3
                    assert abs(entry[3]/np.std(pca_data["PC2"].values)) <= 3
                except AssertionError:
                    continue
                else:
                    plt.scatter(entry[2], entry[3], marker="o", alpha=0.5,
                                c=assigned_colors[md_data[entry[1]]])
            plt.xlabel("Principal Component 1 (2.0 %)", fontsize=10)
            plt.ylabel("Principal Component 2 (0.61 %)", fontsize=10)
            legend_dict = {"Europe": Blues_9.hex_colors[4:][3],
                           "East Asia": Greens_9.hex_colors[4:][3],
                           "South Asia": Purples_9.hex_colors[4:][3],
                           "Sub-Saharan Africa": YlOrBr_9.hex_colors[2:][4]}
            l = [plt.scatter([], [], c=legend_dict[lab], s=50)
                 for lab in legend_dict.keys()]
            plt.legend(l, ["{}".format(lab) for lab in legend_dict.keys()],
                       fontsize=10, loc="best", scatterpoints=3)
            plt.tight_layout()
            if args.savefile:
                plt.savefig(args.savefile, dpi=300, format="svg",
                            bbox_inches="tight", pad_inches=0.25)
            else:
                plt.show()


if __name__ == "__main__":
    sys.exit(main())
