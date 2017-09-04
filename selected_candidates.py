# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 17:58:12 2017

@author: steinrob
"""
import copy

def xray_variables():
    """Gives the names of the two X ray bands included in the data, as well
    as the associated variable for measuring emission in each band.

    :return: list containing band names, and list containing variable names
    """
    xray_bands = ["X-Ray (0.3-2.0) keV", "X-Ray (0.3-10) keV"]
    variable_names = ["luminosity", "countrate"]
    return xray_bands, variable_names


def peak_dict():
    """Returns those candidates which contain a full peak which includes rising
    as well as falling emission. There are three such candidates, selected
    manually through inspection of light curves. Also returns the bands in
    which the peak was well sampled

    :return: list of names, and list of bands.
    """
    p_dict = dict()
    p_dict["names"] = ["PS1-11af", "PS1-10jh", "PTF09ge"]
    p_dict["bands"] = ["r", "g", "i"]
    return p_dict


def dai_and_fang_dict():
    """Returns the list of TDE events which were selected by Dai and Fang
    (https://arxiv.org/abs/1612.00011) as strong candidates for the source of
    Ice Cube neutrinos on the grounds that they were nearby, bright,
    and coincided with Ice Cube operation time.

    :return: List of TDE names
    """
    df_dict = dict()
    df_dict["names"] = ["UGC 03317",
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
    df_dict["bands"] = None

    return df_dict


def agrr_xray_dict():
    agrr_dict = dict()
    agrr_dict["names"] = ["ASASSN-14li",
                          "Swift J1644+57",
                          "Swift J2058+05",
                          "XMMSL1 J0740-85"
                          ]
    agrr_dict["bands"] = None
    return agrr_dict

def agrr_likely_xray_dict():
    agrr_dict = dict()
    agrr_dict["names"] = ["2MASX J0249",
                          "3XMM J152130.7+074916",
                          "IGR J17361-4441",
                          "NGC 247",
                          "OGLE16aaa",
                          "PTF10iya",
                          "SDSSJ1201",
                          "SDSSJ1311",
                          "SDSSJ1323"
                          ]
    agrr_dict["bands"] = None
    return agrr_dict


def agrr_possible_xray_dict():
    agrr_dict = dict()
    agrr_dict["names"] = ["ASASSN-15oi",
                          "D3-13",
                          "LEDA 095953",
                          "NGC 3599",
                          "NGC 5905",
                          "RBS 1032",
                          "RX J1242-11A",
                          "RX J1420+53",
                          "RX J1624+75",
                          "SDSSJ0159",
                          "Swift J1112-82",
                          "WINGS J1348"
                          ]
    agrr_dict["bands"] = None
    return agrr_dict


def agrr_veiled_dict():
    agrr_dict = dict()
    agrr_dict["names"] = ["ASASSN-14ae",
                          "ASASSN-15lh",
                          "D1-9",
                          "D23H-1",
                          "DES14C1kia",
                          "iPTF16fnl",
                          "PS1-10jh",
                          "PS1-11af",
                          "PS1-12yp",
                          "PTF09axc",
                          "PTF09djl",
                          "PTF09ge",
                          "SDSSJ0748",
                          "SDSSJ0952",
                          "SDSSJ1342",
                          "SDSSJ1350",
                          "TDE2"
                          ]
    agrr_dict["bands"] = None
    return agrr_dict


def agrr_unknown_dict():
    agrr_dict = dict()
    agrr_dict["names"] = ["Dougie",
                          "PGC 1185375",
                          "PGC 1190358",
                          "PTF10nuj",
                          "PTF11glr"
                          ]
    agrr_dict["bands"] = None
    return agrr_dict


def agrr_xray_likely_and_veiled_dict():
    agrr_dict = copy.deepcopy(agrr_xray_dict())
    likely = agrr_likely_xray_dict()
    agrr_dict["names"].extend(likely["names"])
    veiled = agrr_veiled_dict()
    agrr_dict["names"].extend(veiled["names"])
    return agrr_dict


def full_dict():
    fulldict = dict()
    fulldict["names"] = None
    fulldict["bands"] = None
    return fulldict


def xray_dict():
    xraydict = dict()
    xraydict["names"] = None
    xraydict["bands"] = xray_variables()
    return xraydict


def return_all_catalogues():
    all_catalogues = dict()
    all_catalogues["full"] = full_dict()
    all_catalogues["daifang"] = dai_and_fang_dict()
    all_catalogues["peaked"] = peak_dict()
    all_catalogues["agrr_xray"] = agrr_xray_dict()
    all_catalogues["agrr_likely"] = agrr_likely_xray_dict()
    all_catalogues["agrr_possible"] = agrr_possible_xray_dict()
    all_catalogues["agrr_veiled"] = agrr_veiled_dict()
    all_catalogues["agrr_unknown"] = agrr_unknown_dict()
    return all_catalogues


def catalogues_to_export():
    export_catalogues = dict()
    export_catalogues["full_TDE"] = full_dict()
    export_catalogues["Dai_Fang_TDE"] = dai_and_fang_dict()
    export_catalogues["AGRR_TDE"] = agrr_xray_likely_and_veiled_dict()
    export_catalogues["AGRR_Xray_TDE"] = agrr_xray_dict()
    return export_catalogues
