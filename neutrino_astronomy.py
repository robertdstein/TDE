import cPickle as pickle
import math
import numpy as np
from classes import FullSet
import astropy
from astropy import units as u
from astropy import cosmology
from astropy.cosmology import WMAP5, WMAP7
# cosmology.set_current(WMAP7)

from export_tde_catalogue import sim_length, core, spectral_indices


def run(data):
    tde_pickle = core + "/Output/pickle_results/Individual_TDEs.pkl"

    with open(tde_pickle) as f:
        tde_dict = pickle.load(f)

    fields = ["Redshift", "RA", "Dec", "Max Date"]

    # gamma = 2.0
    Emax = 10 * u.PeV
    Emin = 100 * u.GeV

    for row in data.data_table:
        name = "".join([x for x in row["Name"] if x != " "])
        if name in tde_dict.keys():
            for gammakey in tde_dict[name].keys():

                gamma = float(gammakey.split("=")[1])

                entry = tde_dict[name][gammakey]

                for field in fields:
                    entry[field] = row[field]

                z = row["Redshift"]

                entry["sens"] = sens = (
                    (entry["analytic"] * (u.TeV ** (gamma-1))) /
                    (u.second * u.cm ** 2)) * (10 ** -12)

                print name, sens,

                model_length = 10 ** 5 * u.second
                window_length = sim_length * (60 * 60 * 24) * u.second

                entry["Time-integrated Sensitivity"] = window_sens = (
                    sens * window_length)
                lumdist = astropy.coordinates.Distance(z=z)
                entry["Luminosity Distance"] = lumdist

                entry["Neutrino Fluence"] = fluence = window_sens * (
                    ((1/(10 * u.TeV)) ** (gamma-1)) -
                    ((1/(10 * u.PeV)) ** (gamma-1)))

                if gamma == 2:
                    e_integral = np.log10(Emax / Emin)
                else:
                    power = 2 - gamma
                    e_integral = ((Emax ** power) - (Emin ** power))/power

                # print e_integral,

                entry["Neutrino Energy"] = etot = ((
                    window_sens * (4 * math.pi * (lumdist.to(u.cm)) ** 2)) *
                                                            e_integral).to(u.erg)

                print etot

            # print z, lumdist
            # raw_input("prompt")

    # for (tde, vals) in tde_dict.items():
    #     print tde, vals["Neutrino Fluence"], np.log10(vals["Neutrino Energy"])
