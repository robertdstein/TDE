
import numpy as np
from astropy import units as u
from classes import FullSet, config_path
import selected_candidates as sc
import os
import ConfigParser
from scipy.interpolate import interp1d

fs_core = "/afs/ifh.de/user/s/steinrob/scratch/flarestack__data/input" \
          "/catalogues/TDEs/"
fs_single_source_dir = fs_core + "individual_TDEs/"


def run(data):
    """Creates a catalogue based on the TDE dataset, for use in stacking
    analysis.

    Produces:
    -   One complete catalogue with every entry containing sufficient data
    -   One catalogue based on the Dai&Fang TDE list
        (list defined in selected_candidates.py)

    :param data: A FullSet object containing TDE data.
    """
    print data

    print "Data contains", len(data.TDEs.__dict__.keys()), "entries"

    export_catalogues = sc.catalogues_to_export()

    for dir_path in [fs_core, fs_single_source_dir]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    fs_custom_dtype = [
        ("ra", np.float),
        ("dec", np.float),
        ("Relative Injection Weight", np.float),
        ("Ref Time (MJD)", np.float),
        ("Start Time (MJD)", np.float),
        ("End Time (MJD)", np.float),
        ('Distance (Mpc)', np.float),
        ('Name', 'a30')
    ]

    for category in set(data.data_table["Category"]):
        print "Category =", category

        mask = data.data_table["Category"] == category

        tab = data.data_table[mask]

        sources = np.empty(np.sum(mask), dtype=fs_custom_dtype)

        sources['ra'] = np.deg2rad(tab["RA"])
        sources['dec'] = np.deg2rad(tab["Dec"])
        sources["Relative Injection Weight"] = np.ones_like(sources['ra'])

        sources['Distance (Mpc)'] = np.array(tab["Luminosity Distance (Mpc)"])
        sources['Name'] = tab["Name"]
        sources["Ref Time (MJD)"] = np.array(tab["Max Date"])
        sources["Start Time (MJD)"] = np.array(
            tab["Window Start Date"])
        sources["End Time (MJD)"] = np.array(
            tab["Window End Date"])
        path = fs_core + "TDE_" + category + "_catalogue.npy"
        np.save(path, sources)

        print "Exporting catalogue with", np.sum(mask), "entries, to", path

    for category in set(data.data_table["Category"]):
        print "Category =", category

        mask = data.data_table["Category"] == category

        tab = data.data_table[mask]

        sources = np.empty(np.sum(mask), dtype=fs_custom_dtype)

        sources['ra'] = np.deg2rad(tab["RA"])
        sources['dec'] = np.deg2rad(tab["Dec"])
        sources["Relative Injection Weight"] = np.ones_like(tab["Weight"])

        sources['Distance (Mpc)'] = np.array(tab["Luminosity Distance (Mpc)"])
        sources['Name'] = tab["Name"]
        sources["Ref Time (MJD)"] = np.array(tab["Max Date"])
        sources["Start Time (MJD)"] = np.array(
            tab["Window Start Date"])
        sources["End Time (MJD)"] = np.array(
            tab["Window End Date"])
        path = fs_core + "TDE_" + category + "_catalogue_weighted.npy"
        np.save(path, sources)

        print "Exporting catalogue with", np.sum(mask), "entries, to", path

    # Loops over catalogues to be created
    for catalogue_name, dictionary in export_catalogues.iteritems():
        ra = []
        dec = []
        distance = []
        discovery_date = []
        max_date = []
        names = []
        starts = []
        ends = []

        variable_names = ["ra_deg", "dec_deg", "lumdist", "mjddisc", "name"]

        if dictionary["names"] is None:
            dictionary["names"] = vars(data.TDEs)

        for tde_name in dictionary["names"]:
            try:
                tde = getattr(data.TDEs, tde_name)
            except AttributeError:
                print sorted(data.TDEs.__dict__.keys(), key=str.lower)
                raise Exception("Name not in catalogue.")
            include = True
            for v in variable_names:
                if not hasattr(tde, v):
                    include = False
                    print tde_name, "is missing", v

            if include:
                ra.append(tde.ra_deg)
                dec.append(tde.dec_deg)

                if str(tde.lumdist.u_value) == "Mpc":
                    distance.append(tde.lumdist.value)
                else:
                    print vars(tde.lumdist)
                    raise Exception("Units of distance are not correct!")
                discovery_date.append(tde.mjddisc)
                max_date.append(tde.mjdmax)
                names.append(tde_name)
                starts.append(data.data_dict[tde_name]["start"])
                ends.append(data.data_dict[tde_name]["end"])

        n_sources = len(ra)

        sources = np.empty(n_sources, dtype=fs_custom_dtype)

        sources['ra'] = np.deg2rad(ra)
        sources['dec'] = np.deg2rad(dec)
        sources["Relative Injection Weight"] = np.ones_like(sources['ra'])

        sources['Distance (Mpc)'] = np.array(distance)
        sources['Name'] = names
        sources["Ref Time (MJD)"] = np.array(max_date)
        sources["Start Time (MJD)"] = np.array(starts)
        sources["End Time (MJD)"] = np.array(ends)
        path = fs_core + catalogue_name + "_catalogue.npy"
        np.save(path, sources)

        print "Exporting catalogue with", n_sources, "entries, to", path

    print "Creating individual source catalogues..."

    for name, tde_dict in data.data_dict.iteritems():

        maxtime = float(tde_dict["mjdmax"])

        sources = np.empty(1, dtype=fs_custom_dtype)

        sources['ra'] = np.deg2rad(tde_dict["ra_deg"])
        sources['dec'] = np.deg2rad(tde_dict["dec_deg"])
        sources["Relative Injection Weight"] = np.array([1.0])
        sources['Distance (Mpc)'] = np.array(tde_dict["dist"])
        sources['Name'] = name
        sources["Ref Time (MJD)"] = np.array(maxtime)
        sources["Start Time (MJD)"] = np.array(tde_dict["start"])
        sources["End Time (MJD)"] = np.array(tde_dict["end"])
        path = fs_single_source_dir + name + "_catalogue.npy"
        np.save(path, sources)
