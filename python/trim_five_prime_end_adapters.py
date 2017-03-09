#!/usr/bin/env python

"""
Abstract: Iterate through all folders and perform adapter trimming on all
          reads. This script uses a sequencing mapping file to map wells to
          sampleids and trim 5'-end adapters of reads. Additionally, reads with
          lengths <200 nt are discarded. For more information on Skewer, visit
          http://www.biomedcentral.com/1471-2105/15/182
Date: 08/31/2015
Author: Akshay Paropkari
"""

import os
import sys
import argparse
import subprocess as sp
from csv import DictReader
try:
    import shlex
except ImportError as ie:
    sys.exit("Please install {} module before executing this script."
             .format(ie))


def create_cmd_to_run(threads, fwd_adapter, rev_adapter, primer, fwd_input_file,
                      rev_input_file, action="D"):
    """
    Creates the command which is given to subprocess call.
    3'-end of sequences are trimmed until desired quality is reached (Phred score = 25).
    Lowest mean quality value allowed before trimming is 25.
    Sequences < 200 base pairs are discarded.

    :type threads: int
    :param threads: Number of threads to use (default = 4)

    :type fwd_adapter: str
    :param fwd_adapter: Forward adapter sequence

    :type rev_adapter: str
    :param rev_adapter: Reverse adapter sequence (not reverse complimented)

    :type primer: str
    :param primer: Primer name

    :type fwd_input_file: filepath
    :param fwd_input_file: Forward reads file (usually R1.fastq)

    :type rev_input_file: filepath
    :param rev_input_file: Reverse reads file (usually R2.fastq)
    """
    if action == "D":
        cmd = "skewer -t {0} -x {1} -y {2} -m head -n -b -l 200 -q 25 -Q 25 -o {3} {4} \
              {5}".format(threads, fwd_adapter, rev_adapter, primer, fwd_input_file,
                          rev_input_file)
    else:
        cmd = "skewer -t {0} -x {1} -y {2} -m head -n -l 200 -q 25 -Q 25 -o {3} {4} {5}".\
              format(threads, fwd_adapter, rev_adapter, primer, fwd_input_file,
                     rev_input_file)
    kwargs = shlex.split(cmd)
    return kwargs


def prog_options():
    parser = argparse.ArgumentParser(
                description="Iterate through each folder in a file and trim 5'-end "
                            "adapters sample fastq files via Skewer. Using adapter-to-"
                            "primer mapping (adapters) and well-to-sampleid mapping "
                            "(ws_map), trim command and output files are generated. This "
                            "scripts discards reads with length <200 nt. 3'-end of "
                            "sequences are trimmed until desired quality (Q25) is reached"
                            ". Lowest mean quality value allowed before trimming is Q25. "
                            "Sequences < 200 base pairs are discarded.")
    parser.add_argument("sample_dir",
                        help="Directory containing sample folders.")
    parser.add_argument("ap_map",
                        help="Path to tab-delimited file containing adapter-primer pair "
                             "per line. Reverse primer adapters must be reverse "
                             "complemented entries.")
    parser.add_argument("ws_map",
                        help="Path to tab-separated sample ID to adapter mapping file. A "
                             "'NA' must be used for empty cells. Format of the file: well"
                             "->sampleid->adapters")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Optionally, users may specify number of threads to be used "
                             "for trimming. Default is 4.")
    parser.add_argument("-d", "--action", choices=["D", "T"], default="D",
                        help="Demultiplex (D) or Trim (T) sequences based on adapters. "
                             "Default action is demultiplexing sequencing.")
    parser.add_argument("-dr", "--dry_run", action="store_true",
                        help="Provide a dry run of the commands which will be executed "
                             "without actually executing the commands.")
    return parser.parse_args()


def main():
    args = prog_options()

    # Checking the necessary input arguments for validity.
    try:
        with open(args.ap_map):
            pass
    except IOError as ioe:
        sys.exit("\nError with adapter file: {}\n".format(ioe))

    try:
        with open(args.ws_map):
            pass
    except IOError as ioe:
        sys.exit("\nError with the well-to-sampleID mapping file: {}\n"
                 .format(ioe))

    # Read adapter sequences from file
    with open(args.ap_map, "rU") as apf:
        adapters = {line.strip().split("\t")[1]: line.strip().split("\t")[0]
                    for line in apf.readlines()[1:]}

    # Read mapping file data
    with open(args.ws_map, "rU") as csvf:
        map_data = [line for line in DictReader(csvf, delimiter="\t")]

    # Iterate through all directories and access their files
    for root, dirs, files in os.walk(args.sample_dir):
        if root != args.sample_dir:
            os.chdir(root)
            well = root
            print("\n{}".format(well))

    # Read the raw FASTQ file
            for file in files:
                if args.action == "D":
                    if file.endswith("R1_001.fastq.gz"):
                        R1_file = file
                    elif file.endswith("R2_001.fastq.gz"):
                        R2_file = file

    # Calculate skewer command for each sample
            for sample in map_data:
                if sample["well"] == well.split("/")[-1].split("-")[0]:
                    # Plate 1 - V1-V3
                    try:
                        if args.action == "T":
                            R1_file = "{}_V1-V3-assigned-A01-pair1.fastq".\
                                      format(sample["P1_sample"])
                            R2_file = "{}_V1-V3-assigned-A01-pair2.fastq".\
                                      format(sample["P1_sample"])
                        command = create_cmd_to_run(args.threads, sample["27F_P1"],
                                                    sample["519R_P1"],
                                                    sample["P1_sample"]+"_V1-V3",
                                                    R1_file, R2_file, args.action)
                        assert "NA" not in command
                    except:
                        pass
                    else:
                        print("Well: {} | Primer: {} | SampleID: {}".format(
                              sample["well"],
                              adapters[sample["27F_P1"]],
                              sample["P1_sample"]))
                        print("{}\n".format(" ".join(command)))
                        if not args.dry_run:
                            out = sp.check_output(command)
                            print("{}".format(out))

                    # Plate 1 - V4-V5
                    try:
                        if args.action == "T":
                            R1_file = "{}_V4-V5-assigned-A01-pair1.fastq".\
                                      format(sample["P1_sample"])
                            R2_file = "{}_V4-V5-assigned-A01-pair2.fastq".\
                                      format(sample["P1_sample"])
                        command = create_cmd_to_run(args.threads, sample["515F_P1"],
                                                    sample["806R_P1"],
                                                    sample["P1_sample"]+"_V4-V5",
                                                    R1_file, R2_file, args.action)
                        assert "NA" not in command
                    except:
                        pass
                    else:
                        print("Well: {} | Primer: {} | SampleID: {}".format(
                              sample["well"],
                              adapters[sample["515F_P1"]],
                              sample["P1_sample"]))
                        print("{}\n".format(" ".join(command)))
                        if not args.dry_run:
                            out = sp.check_output(command)
                            print("{}".format(out))

                    # Plate 2 - V1-V3
                    try:
                        if args.action == "T":
                            R1_file = "{}_V1-V3-assigned-A01-pair1.fastq".\
                                      format(sample["P2_sample"])
                            R2_file = "{}_V1-V3-assigned-A01-pair2.fastq".\
                                      format(sample["P2_sample"])
                        command = create_cmd_to_run(args.threads, sample["27F_P2"],
                                                    sample["519R_P2"],
                                                    sample["P2_sample"]+"_V1-V3",
                                                    R1_file, R2_file, args.action)
                        assert "NA" not in command
                    except:
                        pass
                    else:
                        print("Well: {} | Primer: {} | SampleID: {}".format(
                              sample["well"],
                              adapters[sample["27F_P2"]],
                              sample["P2_sample"]))
                        print("{}\n".format(" ".join(command)))
                        if not args.dry_run:
                            out = sp.check_output(command)
                            print("{}".format(out))

                    # Plate 2 - V4-V5
                    try:
                        if args.action == "T":
                            R1_file = "{}_V4-V5-assigned-A01-pair1.fastq".\
                                      format(sample["P2_sample"])
                            R2_file = "{}_V4-V5-assigned-A01-pair2.fastq".\
                                      format(sample["P2_sample"])
                        command = create_cmd_to_run(args.threads, sample["515F_P2"],
                                                    sample["806R_P2"],
                                                    sample["P2_sample"]+"_V4-V5",
                                                    R1_file, R2_file, args.action)
                        assert "NA" not in command
                    except:
                        pass
                    else:
                        print("Well: {} | Primer: {} | SampleID: {}".format(
                              sample["well"],
                              adapters[sample["515F_P2"]],
                              sample["P2_sample"]))
                        print("{}\n".format(" ".join(command)))
                        if not args.dry_run:
                            out = sp.check_output(command)
                            print("{}".format(out))
    return


if __name__ == "__main__":
    sys.exit(main())
