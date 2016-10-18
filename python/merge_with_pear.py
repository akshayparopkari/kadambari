#!/usr/bin/env python

"""
Abstract: Iterate through all folders and execute PEAR command.
Date: 10/05/2016
Author: Akshay Paropkari
"""

import sys
import shlex
import argparse
import subprocess as sp
from os import chdir, walk


def prog_options():
    parser = argparse.ArgumentParser(
                description="Iterate through each folder in a file and merge adapter "
                            "trimmed Read 1 and Read 2 fastq files using PEAR. "
                            "Unassmebled or discarded reverse reads are not reverse "
                            "complemented while writing out. Consult PEAR documentation "
                            "at https://github.com/xflouris/PEAR for more info about "
                            "optional parameters.")
    parser.add_argument("sample_dir", help="Directory containing sample folders.")
    parser.add_argument("-n", "--min_length", default=200, type=int,
                        help="Filter out sequences less than '--min_length' length. "
                             "Default value is set to 200 base-pair length.")
    parser.add_argument("-y", "--memory", default="4G", type=str,
                        help="Specify  the  amount of memory to be used. Please check "
                             "PEAR documentation. Default value is 4G (4GB).")
    parser.add_argument("-j", "--threads", default=4, type=int,
                        help="Number of threads available for usage. Default is 4.")
    return parser.parse_args()


def main():
    args = prog_options()

    for root, dirs, files in walk(args.sample_dir):
        if root != args.sample_dir:
            chdir(root)
            print("{}".format(root))
            trimmed_list = [file for file in files if "-trimmed-" in file]
            for a, b in zip(trimmed_list[0::2], trimmed_list[1::2]):
                if a.endswith("trimmed-pair1.fastq"):
                    R1 = a
                    R2 = b
                elif a.endswith("trimmed-pair2.fastq"):
                    R1 = b
                    R2 = a
                try:
                    assert R1.split("-trimmed")[0] == R2.split("-trimmed")[0]
                except:
                    continue
                else:
                    output_prefix = R1.split("-trimmed")[0]
                    in1 = "pear -f {0} -r {1} -o {2} -n {3} -q 25 -y {4} -j {5} -k".\
                          format(R1, R2, output_prefix, args.min_length, args.memory,
                                 args.threads)
                    kwargs = shlex.split(in1)
                    print("{}\n".format(" ".join(kwargs)))
                out = sp.check_output(kwargs)
                print("{}".format(out))
                with open("{}_pear_merge_log.txt".format(output_prefix), "w") as outlog:
                    outlog.write("{}".format(out))
    return


if __name__ == "__main__":
    sys.exit(main())
