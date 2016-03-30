#!/usr/bin/env python

"""
Abstract: Iterate through all folders and execute FLASh command.

Date: 08/06/2015

Author: Akshay Paropkari
"""

import sys
import argparse
from os import chdir, walk
import subprocess as sp
from os.path import join, relpath
try:
    import shlex
except ImportError as ie:
    sys.exit("Please install {} module before executing this script."
             .format(ie))


def prog_options():
    parser = argparse.ArgumentParser(
                description="Iterate through each folder in a file and run "
                            "FLASh command on zipped Read 1 and Read 2 fastq "
                            "files.")
    parser.add_argument("sample_dir",
                        help="Directory containing sample folders.")
    return parser.parse_args()


def main():
    args = prog_options()

    for root, dirs, files in walk(args.sample_dir):
        if root != args.sample_dir:
            chdir(root)
            print root.split("/")[-1]

            for file in files:
                if file.endswith("R1_001.fastq.gz"):
                    R1 = file
                elif file.endswith("R2_001.fastq.gz"):
                    R2 = file
            in1 = "flash {} {} -M 300 --cap-mismatch-quals".format(R1, R2)
            kwargs = shlex.split(in1)
            print kwargs
            out = sp.check_output(kwargs)
            print out
            with open("flash_log.txt", "w") as outlog:
                outlog.write("{}".format(out))
    return


if __name__ == "__main__":
    main()
