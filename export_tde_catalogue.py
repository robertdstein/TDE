
import numpy as np
from classes import FullSet
import selected_candidates as sc

root = "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/Catalogue/"


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

        sources = np.empty(n_sources, dtype=[
            ("ra", np.float), ("dec", np.float), ("flux", np.float),
            ("n_exp", np.float), ("weight", np.float),
            ("weight_acceptance", np.float), ("weight_time", np.float),
            ("weight_distance", np.float), ("norm_time", np.float),
            ("global_weight_norm_time", np.float),
            ("discoverydate_mjd", np.float), ("distance", np.float),
            ('name', 'a30'),
            ])

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
