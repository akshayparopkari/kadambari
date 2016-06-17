#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd
import subprocess as sp
import matplotlib.pyplot as plt
from os.path import join, relpath
try:
    import shlex
except ImportError as ie:
    sys.exit("Please install {} module before executing this script.".format(ie))


def prog_options():
    parser = argparse.ArgumentParser(
                description="Iterate through each FASTQ file in a directory and get "
                            "quality stats on input FASTQ files using FASTX Toolkit.")
    parser.add_argument("input_dir", help="Directory containing input files.")
    parser.add_argument("output_dir", help="Directory to save output files in.")
    parser.add_argument("-s", "--save_fig", default=None,
                        help="Path and filename to save the quality plots. Plots will be "
                        "saved in PNG format.")
    return parser.parse_args()


def main():
    args = prog_options()

    try:
        assert os.path.isdir(args.output_dir)
    except AssertionError:
        print "\nCreating output directory...\n"
        os.makedirs(args.output_dir)

    for root, dirs, files in os.walk(args.input_dir):
        for file in files:
            if file.endswith(".fastq") or file.endswith(".fq"):
                input_filename = relpath(join(args.input_dir, file))
                output_filename = relpath(join(args.output_dir, file.split(".")[0]))
                cmd = "fastx_quality_stats -i {} -o {}.txt".format(input_filename,
                                                                   output_filename)
                kwargs = shlex.split(cmd)
                print "\n", kwargs
                out = sp.check_output(kwargs)
                qual_data = pd.read_csv(output_filename+".txt", sep="\t", index_col=0)
                fig = plt.figure(figsize=(10, 7))
                plt.plot(qual_data.index, qual_data["mean"], color="#7570b3",
                         linewidth=2.0, label="Mean Quality Score")
                plt.axhline(y=20, linewidth=2, color="#ff0000", label="Q20")
                plt.ylim([0, 42])
                plt.legend(loc="best")
                plt.title("{}".format(file.split(".")[0]), size=12)
                plt.grid(True, linestyle=":", c="#808080")
                plt.xlabel("Nucleotide Base Number")
                plt.ylabel("Phred Quality Score")
                if args.save_fig:
                    plt.savefig(relpath(join(args.save_fig, file.split(".")[0])) + ".png",
                                dpi=200, facecolor="0.8", format="png",
                                bbox_inches="tight", pad_inches=0.2)
                else:
                    plt.show()
    return


if __name__ == "__main__":
    main()
