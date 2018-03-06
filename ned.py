import os
from sklearn.externals import joblib
import string
import urllib2
from bs4 import BeautifulSoup
import numpy as np
import time

data_dir = "/afs/ifh.de/user/s/steinrob/scratch/TDE_Output/TDE_NED"

names_directory = data_dir + "/names/"
coords_directory = data_dir + "/coords/"


def check_name(name):
    """Checks if a given name has already been used to query NED, and if not,
    then queries NED and saves the file. Otherwise loads the file from the
    previous query.

    :param name: Name of transient
    :return: Dictionary containing data from NED database
    """
    path = names_directory + name + ".html"
    if os.path.isfile(path):
        pass
    else:
        print "Making new entry for name", name
        download_name(name)
    return read_ned_query(path)


def download_name(name):
    """Downloads the file corresponding to the NED query with a given name,
    and writes it to file. Waits one second after querying to prevent
    overloading of NED server.

    :param name: Name of transient
    """
    rep_name = string.replace(name, "+", "%2B")
    sub_name = string.replace(rep_name, " ", "+")

    base_url = "http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname="
    end_url = "&extend=no&out_csys=Equatorial&out_equinox=J2000.0&of=xml_main"
    full_url = base_url + sub_name + end_url

    print full_url

    file_path = names_directory + name + ".html"
    with open(file_path, 'w') as f:
        data = urllib2.urlopen(full_url)
        htmlSource = data.read()
        f.write(htmlSource)
    time.sleep(1)


def check_coordinate(ra, dec):
    """Checks if a given ra/dec has already been used to query NED, and if not,
    then queries NED and saves the file. Otherwise loads the file from the
    previous query.

    :param name: Name of transient
    :return: Dictionary containing data from NED database
    """
    path = coords_directory + str(ra) + "_" + str(dec) + ".html"
    if os.path.isfile(path):
        pass
    else:
        print "Making new entry for RA:", ra, "DEC:", dec
        download_coordinate(ra, dec)
    return read_ned_query(path)


def download_coordinate(ra, dec):
    """Downloads the file corresponding to the NED query with a given ra/dec,
    and writes it to file. Waits one second after querying to prevent
    overloading of NED server. Returns objects within one arc minute of
    position.

    :param ra: Right Ascension
    :param dec: Declination
    """
    base_url = "http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=" \
        "Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon="
    mid_url = str(ra) + "d&lat=" + str(dec)
    end_url = "d&radius=1.0&out_csys=Equatorial&out_equinox=J2000.0&of=xml_main"

    full_url = base_url + mid_url + end_url

    print full_url

    file_path = coords_directory + str(ra) + "_" + str(dec) + ".html"
    with open(file_path, 'w') as f:
        data = urllib2.urlopen(full_url)
        htmlSource = data.read()
        f.write(htmlSource)
    time.sleep(1)


def read_ned_query(path):
    """Reads a query file from the NED database, and returns a dictionary
    containing the data

    :param path: Path of saved file
    :return: Dictionary containing data from NED database
    """
    with open(path) as f:
        parsed_html = BeautifulSoup(f, "xml")

        info = parsed_html.TABLEDATA

        if info != None:
            entry = dict()
            entry["TABLEDATA"] = info

            dt = np.dtype([
                ('No.', np.int),
                ('Object Name', "S50"),
                ("RA(deg)", np.float),
                ("DEC(deg)", np.float),
                ("Type", "S50"),
                ("Velocity", "S50"),
                ("Redshift", np.float),
                ("Redshift Flag", "S50"),
                ("Magnitude And Filter", "S50"),
                ("Distance (arcmin)", "S50"),
                ("References", np.int),
                ("Notes", "S50"),
                ("Photometry Points", np.int),
                ("Positions", np.int),
                ("Redshift Points", np.int),
                ("Diameter Points", np.int),
                ("Associations", np.int)
            ])

            i = 1

            for row in info:
                split = [str(x)[4:-5].strip() for x in row if x != u'\n']
                if len(split) == 0:
                    pass
                else:
                    if split[6] == "":
                        split[6] = np.nan
                    if i == 1:
                        entry["data_table"] = np.array([tuple(split)], dtype=dt)
                        i += 1
                    else:
                        entry["data_table"] = np.append(
                            entry["data_table"], np.array([tuple(split)], dtype=dt))

            # Creates a mask that returns only galaxies, if needed
            acceptable_classes = ["G", '*Cl']
            entry["mask"] = np.array([x in acceptable_classes
                                      for x in entry["data_table"]["Type"]])

            return entry
