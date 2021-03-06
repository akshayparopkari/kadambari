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


def get_vin_info(url):
    """GET VIN number for car search by parsing through result webpage."""
    try:
        c_response = requests.get(url)
        c_response.raise_for_status()
    except Exception as ex:
        print "ERROR!"
        sys.exit(ex)
    else:
        c_soup = BeautifulSoup(c_response.content, "lxml")
    for mdata in c_soup.find_all("p", {"class": "attrgroup"})[1:]:
        for minfo in [tag.text for tag in mdata.find_all("span")]:
            if "VIN:" in minfo:
                vin = minfo.split(": ")[1]
            else:
                vin = "NA"
    return vin

def get_search_results(house_soup, base_url):
    """
    Get results for Cragislist listing. Parse through webpage and get timestamp, title,
    price, address and area info.
    """

    search_results = set()
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
            title = title.strip().translate(None, "*~!")
        except Exception:
            title = "NA"
        link = urljoin(base_url, listing.a["href"])
        vin = get_vin_info(link)
        search_results.add((timestamp, price, address, area, title, link, vin))
    return search_results


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
        URL = urljoin(BASE_URL, "search/cta?sort=pricedsc&auto_title_status=1&max_auto_mi"
                      "les=90000&max_price={2}&min_auto_miles=1000&min_auto_year={3}&"
                      "min_price=3000&postal={1}&search_distance={0}"
                      .format(args.distance, args.zip_code,
                              args.max_price, args.min_year))
    else:
        if args.housing_type == "apartment":
            URL = urljoin(BASE_URL, "search/apa?search_distance={0}&postal={1}&"
                          "min_price=350&max_price={2}&availabilityMode=0&private_room=1".
                          format(args.distance, args.zip_code, args.max_price))
        else:
            URL = urljoin(BASE_URL, "search/roo?search_distance={0}&postal={1}&"
                          "min_price=350&max_price={2}&availabilityMode=0&private_room=1".
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
        search_res = get_search_results(soup, BASE_URL)

    # Write out results
    with open(args.out_fnh, "w") as outf:
        outf.write("Date\tPrice\tAddress\tArea\tTitle\tLink\tVIN\n")
        for res in search_res:
            outf.write("{0}\n".format("\t".join(list(res))))

if __name__ == "__main__":
    sys.exit(main())
