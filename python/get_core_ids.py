#!/usr/bin/env python
"""
:Abstract: Compute core IDs from a BIOM format file.
:Author: Akshay Paropkari
:Date: 03/09/2017
"""
import os
import sys
import argparse
importerrors = []
try:
    import biom
    from biom.util import biom_open
except ImportError:
    importerrors.append("biom")
try:
    from phylotoast.util import parse_map_file, gather_categories
except ImportError:
    importerrors.append("phylotoast")
if len(importerrors) > 0:
    for err in importerrors:
        print("Please install missing module: {}".format(err))
    sys.exit()


def prog_options():
    parser = argparse.ArgumentParser(description="Compute core IDs from a BIOM file.")
    parser.add_argument("in_biomf", help="Input BIOM format file.")
    parser.add_argument("out_fnh", help="Path and name of output directory, if single "
                        "category otherwise path to a directory. Multiple category data "
                        "output files will be saved under their category name.")
    parser.add_argument("map_fnh", help="Mapping file associated with input BIOM file.")
    parser.add_argument("-g", "--group_by",
                        help="Any mapping categories, such as treatment type, that will "
                        "be used to group the data in the  output iTol table. For example"
                        ", one category with three types will result in three data "
                        "columns in the final output. Two categories with three types "
                        "each will result in six data columns. Default is no categories "
                        "and all the data will be treated as a single group.")
    parser.add_argument("-cp", "--core_pct", default=0.8, type=float,
                        help="Value between 0 and 1 which indicates the proportion of "
                        "samples in BIOM format file that an observation should be "
                        "present to be considered as a part of core microbiome. Default "
                        "value os set to 0.8, corresponding to observations found in 80 "
                        "percent of the samples.")
    parser.add_argument("-b", "--biom_out", action="store_true",
                        help="Supply this parameter to output filtered BIOM format file "
                        "for the data. By default, BIOM files will not be created.")
    return parser.parse_args()


def main():
    args = prog_options()

    try:
        biomf = biom.load_table(args.in_biomf)
    except IOError as ioe:
        sys.exit("Error with input BIOM format file: {}".format(ioe))
    else:
        biomf_pa = biomf.pa(inplace=False)  # convert to presence/absence BIOM table
        obs_ids = biomf_pa.ids("observation")

    try:
        mheader, mdata = parse_map_file(args.map_fnh)
    except IOError as ioe:
        sys.exit("Error with input mapping file: {}".format(ioe))
    else:
        if args.group_by:
            sid_cat = gather_categories(mdata, mheader, [args.group_by])
        else:
            sid_cat = gather_categories(mdata, mheader)

    # calculate core
    core_calc = {k: set() for k in sid_cat.keys()}
    for idx in obs_ids:
        for cat, val in sid_cat.iteritems():
            obs_count = 0
            num_of_samples = len(val.sids)
            for sid in val.sids:
                try:
                    assert biomf_pa.get_value_by_ids(idx, sid) == 1
                except AssertionError:
                    continue
                else:
                    obs_count += 1
            try:
                assert obs_count > round(args.core_pct * num_of_samples)
            except AssertionError:
                continue
            else:
                core_calc[cat].add(idx)

    # Check if output directory exists, if not, create it
    try:
        assert os.path.exists(os.path.abspath(args.out_fnh)) is True
    except AssertionError:
        os.makedirs(os.path.abspath(args.out_fnh))
    finally:
        for k, v in core_calc.iteritems():
            print("{0} core IDs in {1}".format(len(v), k))
            idx_filename = os.path.join(os.path.abspath(args.out_fnh),
                                        k+"_80_pct_core_ids.txt")
            with open(idx_filename, "w") as of:
                of.write("{0}".format("\n".join(sorted(v))))
            filtered_biomf = biomf.filter(v, axis="observation", inplace=False)
            if args.biom_out:
                biom_filename = os.path.join(os.path.abspath(args.out_fnh),
                                             k+"_80_pct_core.biom")
                with biom_open(biom_filename, "w") as f:
                    filtered_biomf.to_hdf5(f, "CORE BIOM")


if __name__ == "__main__":
    main()
