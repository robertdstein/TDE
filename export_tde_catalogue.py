
import numpy as np
from astropy import units as u
from classes import FullSet, config_path
import selected_candidates as sc
import os
import ConfigParser
from scipy.interpolate import interp1d

core = "/afs/ifh.de/user/s/steinrob/scratch/The-Flux-Evaluator__Data"
root = core + "/Input/Catalogues/"
single_source_dir = root + "Individual_TDEs/"

sourcepath = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"

fs_core = "/afs/ifh.de/user/s/steinrob/scratch/flarestack__data/input" \
          "/catalogues/TDEs/"
fs_single_source_dir = fs_core + "individual_TDEs/"

# pre_window = 365
# post_window = 100
# length = pre_window + post_window

sim_length = 10

spectral_indices = [1.8, 2.0, 2.5, 3.0]

spectral_indices = ["{0:.1f}".format(x) for x in np.linspace(1.8, 3.0, 13)]


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

    for dir_path in [fs_core, root, single_source_dir, fs_single_source_dir]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    data_dir = sourcepath + "data_configs/"

    tfe_custom_dtype = [
        ("ra", np.float), ("dec", np.float), ("distance", np.float),
        ("flux", np.float), ("weight", np.float),
        ("discoverydate_mjd", np.float), ("Start Time (MJD)", np.float),
        ("End Time (MJD)", np.float), ('name', 'a30'),
    ]

    fs_custom_dtype = [
        ("ra", np.float), ("dec", np.float),
        ("Relative Injection Weight", np.float),
        ("Ref Time (MJD)", np.float),
        ("Start Time (MJD)", np.float),
        ("End Time (MJD)", np.float),
        ('Distance (Mpc)', np.float), ('Name', 'a30'),
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

        # sources = np.empty(n_sources, dtype=tfe_custom_dtype)
        #
        # sources['ra'] = np.deg2rad(ra)
        # sources['dec'] = np.deg2rad(dec)
        # sources['distance'] = np.ones_like(sources['ra'])
        # sources['flux'] = np.ones_like(sources['ra'])
        # normalisation = n_sources * 1.e-9 / np.sum(sources['flux'])
        # sources['flux'] = sources['flux'] * normalisation
        # sources['weight'] = np.ones_like(sources['ra'])
        #
        # sources['discoverydate_mjd'] = np.array(discovery_date)
        # sources['name'] = names
        # sources["Start Time (MJD)"] = np.array(starts)
        # sources["End Time (MJD)"] = np.array(ends)
        # path = root + catalogue_name + "_catalogue.npy"
        # np.save(path, sources)
        #
        # print "Exporting catalogue with", n_sources, "entries, to", path

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

    combo_datasets = []

    data_config = ConfigParser.ConfigParser()

    data_config.read(config_path)

    datastart = float(data_config.get("IC40", "start_mjd"))
    dataend = datastart

    for season in data_config.sections():
        new_maxdate = float(data_config.get(season, "end_mjd"))
        if new_maxdate > dataend:
            dataend = new_maxdate

    config = ConfigParser.ConfigParser()

    if not os.path.isdir(single_source_dir):
        os.makedirs(single_source_dir)

    for name, tde_dict in data.data_dict.iteritems():

        # if name in sc.jetted_dict()["names"]:
        if True:

            # veto = ["Earlier", "Later"]
            veto = []

            hits = [x for x in tde_dict["season"] if x in veto]
            seasons = tde_dict["season"]

            if len(hits) > 0:
                pass
            elif len(seasons) == 0:
                pass
            else:
                seasons = tde_dict["season"]
                print name, tde_dict

                combo_name = "+".join(seasons) + ".ini"
                if len(seasons) > 1:
                    combo_datasets.append(combo_name)

                maxtime = float(tde_dict["mjdmax"])

                # new_catalogue = np.empty(1, dtype=tfe_custom_dtype)
                # new_catalogue['ra'] = np.deg2rad(tde_dict["ra_deg"])
                # new_catalogue['dec'] = np.deg2rad(tde_dict["dec_deg"])
                # new_catalogue['distance'] = np.array([1.0])
                # new_catalogue['flux'] = np.array([1.e-9])
                # new_catalogue['weight'] = np.array([1.0])
                # new_catalogue['discoverydate_mjd'] = np.array(maxtime)
                # new_catalogue['name'] = name
                # new_catalogue["Start Time (MJD)"] = np.array(tde_dict["start"])
                # new_catalogue["End Time (MJD)"] = np.array(tde_dict["end"])
                # cat_path = single_source_dir + name + ".npy"
                #
                # np.save(cat_path, new_catalogue)

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

                # t0 = float(tde_dict["start"]) - maxtime
                # length = float(tde_dict["end"]) - float(tde_dict["start"])
                #
                # for gamma in spectral_indices:
                #     section_name = "Individual_TDEs/" + \
                #                    "".join([x for x in name if x != " "]) + \
                #                    "_gamma=" + str(gamma) + "/"
                #
                #     config.add_section(section_name)
                #     config.set(section_name, "UseEnergy", True)
                #     config.set(section_name, "FitGamma", True)
                #     config.set(section_name, "FixedGamma", gamma)
                #     config.set(section_name, "InjectionGamma", gamma)
                #     config.set(section_name, "UseTime", True)
                #     config.set(section_name, "SimTimeModel", "Box")
                #     config.set(section_name, "SimTimeParameters",
                #                {"t0": t0, "length": sim_length})
                #     config.set(section_name, "ReconTimeModel", "Box")
                #     config.set(section_name, "ReconTimeParameters",
                #                {"t0": t0, "length": length})
                #     config.set(section_name, "FitWeights", False)
                #     config.set(section_name, "UseBox", False)
                #     config.set(section_name, "CatName", cat_path)
                #     config.set(section_name, "DataConfig", combo_name)
                #
                #     sens_save_path = core + "/Output/PS_k_Sensitivities/" + \
                #         seasons[0] + ".npy"
                #
                #     data = np.load(sens_save_path)[0].T
                #
                #     f = interp1d(data[0], data[1])
                #     sens = 150 * f(np.sin(np.deg2rad(tde_dict["dec_deg"])))
                #
                #     if float(tde_dict["start"]) < datastart:
                #         delta = datastart - float(tde_dict["start"])
                #         frac = (length - delta)/length
                #
                #         sens *= 1/frac
                #
                #     elif float(tde_dict["end"]) > dataend:
                #         delta = float(tde_dict["end"]) - dataend
                #         frac = (length - delta)/length
                #
                #         sens *= 1/frac
                #
                #     sens *= 12 ** ((float(gamma) - 2.0)/0.25)
                #
                #     sens *= length / 450
                #
                #     config.set(section_name, "MaxK", sens)

    # analysis_path = sourcepath + "analysis_config/individual_TDEs.ini"
    #
    # with open(analysis_path, "wb") as f:
    #     config.write(f)
    #
    # combo_datasets = list(set(combo_datasets))
    #
    # for dataset in combo_datasets:
    #     seasons = dataset[:-4].split("+")
    #
    #     new_config = ConfigParser.ConfigParser()
    #
    #     for season in seasons:
    #         old = ConfigParser.ConfigParser()
    #         old_path = data_dir + season + ".ini"
    #         old.read(old_path)
    #
    #         new_config.add_section(season)
    #
    #         for (var, val) in old.items(season):
    #             new_config.set(season, var, val)
    #
    #     new_path = data_dir + dataset
    #
    #     with open(new_path, "wb") as file:
    #         new_config.write(file)
