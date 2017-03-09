#!/usr/bin/env python
"""
:Abstract: Populate the SQL database table with OTU metadata.
:Author: Akshay Paropkari
:Date: 12/28/2016
"""

import sys
import uuid
import argparse
try:
    import pandas as pd
except ImportError as ioe:
    sys.exit("Please install missing pandas module.\n{}".format(ioe))


def handle_program_options():
    parser = argparse.ArgumentParser(description="Populate a database wiht OTU metadata.")
    parser.add_argument("gramox_master_fnh",
                        help="Path to master gramox data file.")
    parser.add_argument("gramox_db_fnh",
                        help="Path and name of output file with SQL commands.")
    return parser.parse_args()


def main():

    args = handle_program_options()

    # Read in OTU metadata
    try:
        otu_md = pd.read_csv(args.gramox_master_fnh, sep="\t")
    except IOError as ioe:
        err_msg = "\nError opening master OTU metadata file: {}\n"
        sys.exit(err_msg.format(ioe))

    # Iterate over metadata dataframe and generate db entries for each column
    with open(args.gramox_db_fnh, "w") as outf:
        for rows in otu_md.iterrows():
            row = rows[1]
            genusID = uuid.uuid4()
            speciesID = uuid.uuid4()
            citationID = uuid.uuid4()
            outf.write("insert into genus (genusID, name) VALUES ({}, '{}');\n".
                       format(genusID, row["Genus"]))
            outf.write("insert into citation (citationID, url) VALUES ({}, '{}');\n".
                       format(citationID, row["Source"]))
            outf.write("insert into species (speciesID, name) VALUES ({}, '{}');\n".
                       format(speciesID, row["Species"]))
            outf.write("insert into oxygenrequirement (speciesID, citationID) VALUES ({},"
                       "'{}');\n".format(speciesID, citationID))
            outf.write("insert into gramstatus (speciesID, citationID) VALUES ({}, '{}');\n".
                       format(speciesID, citationID))


if __name__ == "__main__":
    sys.exit(main())
