
import numpy as np
from classes import FullSet, config_path
import selected_candidates as sc
import os
import ConfigParser
from scipy.interpolate import interp1d

core = "/afs/ifh.de/user/s/steinrob/scratch/The-Flux-Evaluator__Data"
root = core + "/Input/Catalogues/"

sourcepath = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"

pre_window = 365
post_window = 100
length = pre_window + post_window

sim_length = 10

spectral_indices = [1.8, 2.0, 2.5, 3.0]
spectral_indices = [1.8, 2.0]


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

    data_dir = sourcepath + "data_configs/"

    custom_dtype = [
            ("ra", np.float), ("dec", np.float), ("flux", np.float),
            ("n_exp", np.float), ("weight", np.float),
            ("weight_acceptance", np.float), ("weight_time", np.float),
            ("weight_distance", np.float), ("norm_time", np.float),
            ("global_weight_norm_time", np.float),
            ("discoverydate_mjd", np.float), ("distance", np.float),
            ('name', 'a30'),
            ]

    # Loops over catalogues to be created
    for catalogue_name, dictionary in export_catalogues.iteritems():
        ra = []
        dec = []
        distance = []
        discovery_date = []
        names = []

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
                distance.append(tde.lumdist)
                discovery_date.append(tde.mjddisc)
                names.append(tde_name)

        n_sources = len(ra)

        sources = np.empty(n_sources, dtype=custom_dtype)

        sources['ra'] = np.deg2rad(ra)
        sources['dec'] = np.deg2rad(dec)
        sources['distance'] = np.ones_like(sources['ra'])
        sources['flux'] = np.ones_like(sources['ra'])
        normalisation = n_sources * 1.e-9 / np.sum(sources['flux'])
        sources['flux'] = sources['flux'] * normalisation
        sources['weight'] = np.ones_like(sources['ra'])

        sources['discoverydate_mjd'] = np.array(discovery_date)
        sources['name'] = names
        path = root + catalogue_name + "_catalogue.npy"
        np.save(path, sources)

        print "Exporting catalogue with", n_sources, "entries, to", path

    print "Creating individual source catalogues..."

    combo_datasets = []

    single_source_dir = root + "Individual_TDEs/"

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

        veto = ["Earlier", "Later"]

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

            new_catalogue = np.empty(1, dtype=custom_dtype)
            new_catalogue['ra'] = np.deg2rad(tde_dict["ra_deg"])
            new_catalogue['dec'] = np.deg2rad(tde_dict["dec_deg"])
            new_catalogue['distance'] = np.array([1.0])
            new_catalogue['flux'] = np.array([1.e-9])
            new_catalogue['weight'] = np.array([1.0])
            new_catalogue['discoverydate_mjd'] = np.array(maxtime)
            new_catalogue['name'] = name
            cat_path = single_source_dir + name + ".npy"

            np.save(cat_path, new_catalogue)

            for gamma in spectral_indices:

                section_name = "Individual_TDEs/" + \
                               "".join([x for x in name if x != " "]) + \
                               "_gamma=" + str(gamma) + "/"

                config.add_section(section_name)
                config.set(section_name, "UseEnergy", True)
                config.set(section_name, "FitGamma", True)
                config.set(section_name, "FixedGamma", gamma)
                config.set(section_name, "UseTime", True)
                config.set(section_name, "SimTimeModel", "Box")
                config.set(section_name, "SimTimeParameters",
                           {"t0": -pre_window, "length": sim_length})
                config.set(section_name, "ReconTimeModel", "Box")
                config.set(section_name, "ReconTimeParameters",
                           {"t0": -pre_window, "length": length})
                config.set(section_name, "FitWeights", False)
                config.set(section_name, "UseBox", False)
                config.set(section_name, "CatName", cat_path)
                config.set(section_name, "DataConfig", combo_name)

                sens_save_path = core + "/Output/PS_k_Sensitivities/" + \
                    seasons[0] + ".npy"

                data = np.load(sens_save_path)[0].T

                f = interp1d(data[0], data[1])
                sens = 50 * f(np.sin(np.deg2rad(tde_dict["dec_deg"])))

                start_window = maxtime - pre_window

                end_window = maxtime + post_window

                if start_window < datastart:
                    delta = datastart - start_window
                    frac = (length - delta)/length

                    sens *= 1/frac

                elif end_window > dataend:
                    delta = end_window - dataend
                    frac = (length - delta)/length

                    sens *= 1/frac

                config.set(section_name, "MaxK", sens)

    analysis_path = sourcepath + "analysis_config/individual_TDEs.ini"

    with open(analysis_path, "wb") as f:
        config.write(f)

    combo_datasets = list(set(combo_datasets))

    for dataset in combo_datasets:
        seasons = dataset[:-4].split("+")

        new_config = ConfigParser.ConfigParser()

        for season in seasons:
            old = ConfigParser.ConfigParser()
            old_path = data_dir + season + ".ini"
            old.read(old_path)

            new_config.add_section(season)

            for (var, val) in old.items(season):
                new_config.set(season, var, val)

        new_path = data_dir + dataset

        with open(new_path, "wb") as file:
            new_config.write(file)
