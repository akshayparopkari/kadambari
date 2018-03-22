#!/usr/bin/env python

"""
:Abstract: Run through multiple files and extract ancestry proportion for Chr21 from 1000
           Genomes project files.
:Author: Akshay Paropkari
:Date: 03/15/2018
"""


import os
import sys
import argparse
import subprocess as sp
from collections import defaultdict
err = []
try:
    import shlex
except ImportError:
    err.append("shlex")
try:
    import pandas as pd
except ImportError:
    err.append("pandas")
try:
    from palettable.cartocolors.qualitative import Safe_6    # 6 populations
except ImportError:
    err.append("palettable")
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc("font", family="Arial")
    mpl.rc("xtick", labelsize=10)  # set X axis ticksize
    mpl.rc("ytick", labelsize=10)  # set Y axis ticksize
except ImportError:
    err.append("matplotlib")
try:
    assert len(err) == 0
except AssertionError:
    for error in err:
        print("Please install {}".format(error))
    sys.exit()


def dir_iterator(sample_directory):
    """Iterate over a directory contents/files and yield files."""
    for root, dirs, files in os.walk(sample_directory):
        return root.strip("/"), files


def get_variant_counts(var_file):
    """Parse total variant counts file."""
    var_counts = defaultdict(dict)
    with open(var_file, "r") as vf:
        for line in vf.readlines()[1:]:
            line = line.strip().split("\t")
            var_counts[line[0]] = {line[1]: line[2]}
    return var_counts


def get_proportion_counts(dir_struct):
    """
    Run through a folder of .bed files and collect genome proportion on ancestry data
    for Chr21.

    :type dir: str
    :param dir: Directory containing .bed files. Input structure must conform to the
                following convention:
                ACB_bed/Sample1_A_final.bed, ACB_bed/Sample1_b_final.bed, ...

    :return: Data will be saved in a folder called "lai" in main directory as
             ACB_bed/lai/Sample1.txt, ACB_bed/lai/Sample2.txt, ...
    """
    for root, files in dir_struct.items():
        try:
            assert os.path.exists(os.path.join(root, "lai"))
        except AssertionError:
            os.makedirs(os.path.join(root, "lai"))
            output_dir = os.path.join(root, "lai")
        else:
            output_dir = os.path.join(root, "lai")
        for fnh1, fnh2 in zip(sorted(files)[0::2], sorted(files)[1::2]):
            file1 = os.path.join(root, fnh1)
            file2 = os.path.join(root, fnh2)
            try:
                assert fnh1.split("_")[0] == fnh2.split("_")[0]
            except AssertionError:
                print("ERROR! Files are not for same genome: {} {}".format(file1, file2))
                sys.exit()
            else:
                genome = fnh1.split("_")[0]
                in1 = "awk -F'\\t' 'BEGIN {OFS=FS} $1 == 21 {print $1, $2, $3, $4, $3 - "\
                      "$2}' %s %s" % (file1, file2)
                kwargs = shlex.split(in1)
                print(in1)
                out = sp.check_output(kwargs)
                with open(os.path.join(output_dir, genome)+".txt", "w") as outf:
                    outf.write(out.decode())


def get_plot_data(var_file, dir_struct):
    """
    Iterate through ancestry count files for each population, and save formatted local
    ancestry data to tab-separated file.
    """
    try:
        variant_counts = get_variant_counts(var_file)
    except Exception as ex:
        sys.exit("Error getting variant counts: {}".format(ex))
    else:
        for root, files in dir_struct.items():
            try:
                assert os.path.exists(os.path.join(root, "lai"))
            except AssertionError:
                sys.exit("Please run -gc before trying to. visualize the results.")
            else:
                lai_data = defaultdict(dict)
                input_dir = os.path.join(root, "lai")
                samples = [file.split(".")[0]
                           for root, dirs, files in os.walk(input_dir)
                           for file in files]
                population = root.split("_")[0]
                for sid in samples:
                    try:
                        lai_data[sid] = {"variants": variant_counts[sid][population],
                                         "AFR": 0, "NAT": 0, "EUR": 0,
                                         "source": root.split("_")[0]}
                    except Exception as ex:
                        continue
                for root, dirs, files in os.walk(input_dir):
                    for file in files:
                        fpath = os.path.join(input_dir, file)
                        sampleid = file.split(".")[0]
                        inf = pd.read_csv(fpath, sep="\t", header=0,
                                          names=["genome", "start", "end",
                                                 "population", "difference"])
                        inf["genome"] = inf["genome"].map({21: sampleid})
                        inf["prop"] = inf["difference"]/inf["difference"].sum()
                        inf = inf.groupby(["population"]).sum()
                        for row in inf.itertuples():
                            try:
                                lai_data[sampleid][row[0]] = row[-1]
                            except Exception:
                                continue
                lai_data_df = pd.DataFrame.from_dict(lai_data, orient="index")
                lai_data_df.to_csv("{}_plot.txt".format(root.split("_")[0]), sep="\t")


def prog_options():
    parser = argparse.ArgumentParser(
                description="Iterate through each folder in a file and process ancestry "
                            "information for Chr21 using 1000 Genomes Project bed files.")
    parser.add_argument("-sd", "--sample_dir", nargs='+',
                        help="Folder(s) containing input files.")
    parser.add_argument("-vf", "--var_file", help="Input tab-separated file containing "
                        "variant counts. Format of this file should follow this "
                        "convention: SampleID Population Variant Counts.")
    parser.add_argument("-gc", "--get_counts", action="store_true",
                        help="Supply this parameter to iterate through .bed files and "
                        "calculate ancestry proportion.")
    parser.add_argument("-gpd", "--get_plot_data", action="store_true",
                        help="Supply this option to calculate the ancestry proportions.")
    parser.add_argument("-pd", "--plot_data", nargs='+',
                        help="Supply this option to visualize the ancestry proportions.")
    parser.add_argument("-sp", "--save_plot", action="store_true",
                        help="Supply this option to save the ancestry proportions plot. "
                        "Output plot will be saved as ancestry.svg in the working "
                        "directory.")
    parser.add_argument("-dr", "--dry_run", action="store_true",
                        help="Dry run to see which directories and files will be used as "
                        "input. Use this flag to check that correct files are supplied. "
                        " Default is False.")
    return parser.parse_args()


def main():
    args = prog_options()

    # Check if input is single or multiple directory and get contents in a dict
    if args.sample_dir:
        dir_struct = {}
        try:
            assert len(args.sample_dir) == 1
        except AssertionError:
            if args.dry_run:  # dry run execution, else skip
                for input in args.sample_dir:
                    for root, dirs, files in os.walk(input):
                        print("\nroot:", root)
                        print("dir:", dirs)
                        print("files:", files)
                sys.exit()
            for input in args.sample_dir:
                name, files = dir_iterator(input)
                dir_struct[name] = files
        else:
            folder = args.sample_dir[0]
            if args.dry_run:  # dry run execution, else skip
                for root, dirs, files in os.walk(folder):
                    print("\nroot:", root)
                    print("dir:", dirs)
                    print("files:", files)
                    sys.exit()
            name, files = dir_iterator(folder)
            dir_struct[name] = files

    # Iterate through .bed files and save ancestry proportion
    if args.get_counts:
        get_proportion_counts(dir_struct)

    # Get variant counts for each genome
    if args.get_plot_data:
        get_plot_data(args.var_file, dir_struct)

    # Plot number of variants vs ancestry (NAT, AFR, EUR) per genome
    if args.plot_data:
        all_data = pd.DataFrame()
        colors = Safe_6.hex_colors
        pop_labels = tuple(["PEL", "MXL", "CLM", "PUR", "ASW", "ACB"])
        assigned_colors = {a:b for a, b in zip(pop_labels, colors)}
        for infile in args.plot_data:
            try:
                inf = pd.read_csv(infile, sep="\t")
            except Exception as ex:
                sys.exit("Error reading input file: {}".format(ex))
            else:
                all_data = pd.concat([all_data, inf], ignore_index=True,
                                     verify_integrity=True)
        fig, axarr = plt.subplots(3, sharex=True, figsize=(10,20))
        for entry in all_data.itertuples():
            try:
                axarr[0].scatter(entry[3], entry[8], marker=".",
                                 c=assigned_colors[entry[7]])
                axarr[1].scatter(entry[1], entry[8], marker=".",
                                 c=assigned_colors[entry[7]])
                axarr[2].scatter(entry[2], entry[8], marker=".",
                                 c=assigned_colors[entry[7]])
            except Exception:
                continue
        axarr[0].set_title("Native American Ancestry", fontdict={"fontsize": 18})
        axarr[1].set_title("African Ancestry", fontdict={"fontsize": 18})
        axarr[2].set_title("European Ancestry", fontdict={"fontsize": 18})
        for ax in axarr:
            ax.set_ylabel("Number of Variants per Genome in Chr21",
                          fontdict={"fontsize": 12})
            ax.grid(axis="x", linestyle=":")
            l = [plt.scatter([], [], c=assigned_colors[lab], s=150) for lab in pop_labels]
            ax.legend(l, ["{}".format(lab) for lab in pop_labels], fontsize=18,
                      bbox_to_anchor=(1,1), loc="upper left", scatterpoints=3)
        plt.savefig("ancestry.svg", dpi=300, format="svg", bbox_inches="tight",
                    pad_inches=0.1)

if __name__ == "__main__":
    sys.exit(main())
