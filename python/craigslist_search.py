#!/usr/bin/env python

"""
Abstract: Web Scraping for house and car listings on Craigslist.
Author: Akshay Paropkari
"""
import sys
import logging
import argparse
from urlparse import urljoin
importerrors = []
try:
    import requests
except ImportError:
    importerrors.append("requests")
try:
    from bs4 import BeautifulSoup
except ImportError:
    importerrors.append("beautifulsoup4")
if len(importerrors) != 0:
    for err in importerrors:
        print "Import Error: {}".format(err)
    sys.exit()


def get_listings(outfile, house_soup, base_url):
    """
    Get results for house/room listing. Parse through webpage and get timestamp, title,
    price, address and area info.
    """
    with open(outfile, "w") as outf:
        outf.write("TimeStamp\tTitle\tAddress\tPrice\tArea\tLink\n")
        for listing in house_soup.find_all("li", {"class": "result-row"}):
            timestamp = str(listing.find("time", {"class": "result-date"})).\
                        split(" ")[2].split('"')[1]
            try:
                price = str(listing.find("span", {"class": "result-price"})).\
                        split(">")[1].split("<")[0]
            except Exception:
                price = "NA"
            try:
                address = str(listing.find("span", {"class": "result-hood"})).\
                          split("(")[1].split(")")[0]
            except Exception:
                address = "NA"
            try:
                area = str(listing.find("span", {"class": "housing"})).split()[2].\
                       split("<")[0]
            except Exception:
                area = "NA"
            try:
                title = str(listing.find("a", {"class": "result-title hdrlnk"})).\
                        split(">")[1].split("<")[0]
                title = title.strip()
            except Exception:
                title = "NA"
            link = urljoin(base_url, listing.a["href"])
            outf.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".
                       format(timestamp, title, address, price, area, link))


def handle_program_options():
    parser = argparse.ArgumentParser(description="Web scraping Craigslist listing and "
                                     "getting results in a tsv file.")
    parser.add_argument("city_name", help="Name of the city you want to search."
                        "[REQUIRED]")
    parser.add_argument("-c", "--car_search", action="store_true",
                        help="Set to get listings of cars rather than housing. Default is"
                        " not reset.")
    parser.add_argument("-t", "--housing_type", default="apartment",
                        choices=["apartment", "room"],
                        help="The type of housing you are looking for. Default is set to "
                             "apartment")
    parser.add_argument("-d", "--distance", default=10, type=int,
                        help="Return results based on distance from postal zip code. "
                        "Default is set to 10 miles. This value must be an integer.")
    parser.add_argument("-z", "--zip_code", type=int,
                        help="Postal zip code to search around. This value must be an "
                             "integer.")
    parser.add_argument("-mp", "--max_price", type=int,
                        help="Maximum price for the apartment/room. This value must be an"
                             " integer.")
    parser.add_argument("-my", "--min_year", type=int,
                        help="Minimum make year for car. If provided, cars manufactured "
                        "before this year will not be listed in result. This value must "
                        "be an integer.")
    parser.add_argument("out_fnh", help="Output file name and path. [REQUIRED]")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Handle logging info
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s | %(levelname)s | %(funcName)s | %(message)s")

    # Obtain relevant URL
    BASE_URL = "https://{0}.craigslist.org/".format(args.city_name.lower())
    if args.car_search:
        URL = urljoin(BASE_URL, "search/cta?search_distance={0}&postal={1}&max_price={2}&"
                      "&min_auto_year={3}".
                      format(args.distance, args.zip_code, args.max_price, args.min_year))
    else:
        if args.housing_type == "apartment":
            URL = urljoin(BASE_URL, "search/apa?search_distance={0}&postal={1}&"
                          "max_price={2}&availabilityMode=0&private_room=1".
                          format(args.distance, args.zip_code, args.max_price))
        else:
            URL = urljoin(BASE_URL, "search/roo?search_distance={0}&postal={1}&"
                          "max_price={2}&availabilityMode=0&private_room=1".
                          format(args.distance, args.zip_code, args.max_price))

    # Download and read HTML page content
    try:
        response = requests.get(URL)
        response.raise_for_status()
    except Exception as ex:
        print "ERROR!"
        sys.exit(ex)
    else:
        soup = BeautifulSoup(response.content, "lxml")
        get_listings(args.out_fnh, soup, BASE_URL)


if __name__ == "__main__":
    sys.exit(main())
