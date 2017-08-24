# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 17:58:12 2017

@author: steinrob
"""


def return_candidates_with_peaks():
    """Returns those candidates which contain a full peak which includes rising
    as well as falling emission. There are three such candidates, selected
    manually through inspection of light curves. Also returns the bands in
    which the peak was well sampled

    :return: list of names, and list of bands.
    """
    names = ["PS1-11af", "PS1-10jh", "PTF09ge"]
    bands = ["r", "g", "i"]
    return names, bands


def xray_names():
    """Gives the names of the two X ray bands included in the data, as well
    as the associated variable for measuring emission in each band.

    :return: list containing band names, and list containing variable names
    """
    xray_bands = ["X-Ray (0.3-2.0) keV", "X-Ray (0.3-10) keV"]
    variable_names = ["luminosity", "countrate"]
    return xray_bands, variable_names


def dai_and_fang_list():
    """Returns the list of TDE events which were selected by Dai and Fang
    (https://arxiv.org/abs/1612.00011) as strong candidates for the source of
    Ice Cube neutrinos on the grounds that they were nearby, bright,
    and coincided with Ice Cube operation time.

    :return: List of TDE names
    """
    names = ["UGC 03317",
             "PGC 1185375",
             "PGC 1190358",
             "PGC 015259",
             "iPTF16fnl",
             "XMMSL1 J0740-85",
             "ASASSN-15oi",
             "ASASSN-14li",
             "ASASSN-14ae",
             "Swift J1644+57"
             ]

    return names
