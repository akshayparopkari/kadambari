#!/usr/bin/env python

"""
Abstract: Web Scraping for apartment listings on Craigslist.

Author: Akshay Paropkari
"""
importerrors = []
import sys
import argparse
from urlparse import urljoin
try:
    import requests
except ImportError as ie:
    importerrors.append("requests")
try:
    import pandas as pd
except ImportError as ie:
    importerrors.append("pandas")
try:
    from bs4 import BeautifulSoup
except ImportError as ie:
    importerrors.append("beautifulsoup4")
if len(importerrors) != 0:
    for package in importerrors:
        print "Please install %s package." % package
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description="Web scraping Craigslist \
                                                  apartment listing and \
                                                  getting results in a tab-\
                                                  separated file.")
    parser.add_argument("listings_url",
                        help="URL of the Craigslist search result page.")
    parser.add_argument("out_fnh", help="Output file name and path.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Obtain all the urls
    URL = args.listings_url
    BASE_URL = "https://columbus.craigslist.org/"

    response = requests.get(URL)
    soup = BeautifulSoup(response.content, "lxml")

    data_result = []

    for listing in soup.find_all("p", {"class": "row"}):
        timestamp = str(listing.find("span", {"class": "pl"}).contents[1])\
            .split('="')[2].split('"')[0]
        if listing.find("span", {"class": "price"}) is not None:
            price = listing.find("span", {"class": "price"}).contents[0]
        if listing.find("span", {"class": "pnr"}) is not None:
            address = str(listing.find("span", {"class": "pnr"}).contents[1])
            try:
                address = address.split("(")[1].split(")")[0]
            except IndexError:
                address = "Address information not available for this listing."
        link = urljoin(BASE_URL, listing.a["href"])
        if listing.find("span", {"class": "housing"}) is not None:
            housing = str(listing.find("span", {"class": "housing"})
                          .contents[0]).split("/")[1].strip()
        data_result.append({"TimeStamp": timestamp, "Address": address,
                            "Price": price, "Unit Info": housing,
                            "Link": link})
        result_df = pd.DataFrame(data_result)
        result_df.to_csv(args.out_fnh, sep='\t', index=False)

if __name__ == "__main__":
    main()
