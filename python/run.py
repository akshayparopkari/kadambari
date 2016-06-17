#!/usr/bin/env python
import os
import sys
import argparse
import subprocess as sp
try:
    import shlex
except ImportError as ie:
    sys.exit("Please install {} module before executing this script.".format(ie))


def prog_options():
    parser = argparse.ArgumentParser(
                description="Iterate through each folder in a file and run some command.")
    parser.add_argument("sample_dir", help="Directory containing sample folders.")
    return parser.parse_args()


def main():
    args = prog_options()

    for root, dirs, files in os.walk(args.sample_dir):
        print "root:", root
        print "dir:", dirs
        print "files:", files
    return


if __name__ == "__main__":
    main()
