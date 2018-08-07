import cPickle as pickle
import math
import numpy as np
from classes import FullSet
from tabulate import tabulate
import astropy
from astropy import units as u
import selected_candidates as sc
from astropy import cosmology
from astropy.cosmology import WMAP5, WMAP7
# cosmology.set_current(WMAP7)

from export_tde_catalogue import sim_length, core, spectral_indices
import export_table_to_tex


def run(data):
    tde_pickle = core + "/Output/pickle_results/Individual_TDEs.pkl"

    with open(tde_pickle) as f:
        sens_dict = pickle.load(f)

    fields = ["Redshift", "RA", "Dec", "Max Date"]

    # gamma = 2.0
    Emax = (10 * u.PeV).to(u.GeV)
    Emin = (100 * u.GeV).to(u.GeV)
    E0 = 1 * u.GeV

    f_pi = 0.1
    waxmann_bachall = (3./8.) * f_pi

    baryon_loading = 100

    shared = [
        ('Name', "S50"),
        ("Redshift", np.float),
        ("RA", np.float),
        ("Dec", np.float),
        ("Max Date", np.float)
    ]

    specific = [
        "Sensitivity", "Time Integrated Sensitivity",
        "Integrated Sensitivity", "Total N Neutrino", "pretty Neutrino Energy",
        "Neutrino Energy", "CR Energy", "Rad Energy",
        "Average Neutrino Luminosity"
        ]

    extra = [(x + " (" + str(y) + ")", np.float) for y
             in sorted(sens_dict.itervalues().next().keys()) for x in specific]

    full = shared + extra

    dt = np.dtype(full)

    table = np.zeros(len(sens_dict), dtype=dt)

    i = 0

    for row in data.data_table:
        name = "".join([x for x in row["Name"] if x != " "])

        # TDE_entry = getattr(data.TDEs, row["Name"])

        if name in sens_dict.keys():

            base = ["Name", "Redshift", "RA", "Dec", "Max Date"]

            new_row = [row[x] for x in base]

            if row["Name"] in sc.jetted_dict().keys():
                new_name = r"\b{" + base["Name"] + "}"
                new_row[0] = new_name

            for gammakey in sorted(sens_dict[name].keys()):
                gamma = float(gammakey.split("=")[1])

                print gamma

                entry = sens_dict[name][gammakey]

                for field in fields:
                    entry[field] = row[field]

                z = row["Redshift"]

                entry["Sensitivity"] = sens = (
                    (entry["analytic"] * (10 ** -9)) * (
                    1. / (u. GeV * u.second * u.cm ** 2)))

                # model_length = 10 ** 5 * u.second
                window_length = sim_length * (60 * 60 * 24) * u.second

                entry["Time Integrated Sensitivity"] = window_sens = (
                    sens * window_length)

                # entry["pretty TIS"] = window_sens.value
                lumdist = astropy.coordinates.Distance(z=z)
                entry["Luminosity Distance"] = lumdist

                entry["Illuminated Area"] = area = (
                    4 * math.pi * (lumdist.to(u.cm)) ** 2)

                phi_power = 1 - gamma
                # phi_integral = (1 / phi_power) * (
                #     (Emax ** phi_power) - (Emin ** phi_power))
                #
                # phi_integral = (1 / phi_power) * (
                #     (Emax ** phi_power) - (Emin ** phi_power)) * u.GeV

                phi_integral = ((1./phi_power) * (E0 ** gamma) * (
                    (Emax ** phi_power) - (Emin ** phi_power))).value * u.GeV

                entry["Integrated Sensitivity"] = dNdA = \
                    (window_sens * phi_integral).to(u.cm**-2)

                entry["Total N Neutrino"] = dNdA * area

                if gamma == 2:
                    e_integral = np.log10(Emax / Emin)
                    e_integral = np.log10(Emax / Emin) * (u.GeV ** 2)
                else:
                    power = 2 - gamma

                    print (1./power) * (E0 ** gamma)
                    print (Emax ** power) - (Emin ** power)

                    print ((1./power) * (E0 ** gamma) * (
                            (Emax ** power) - (Emin ** power)))

                    # Get around astropy power rounding error (does not give
                    # EXACTLY 2)

                    e_integral = ((1./power) * (E0 ** gamma) * (
                            (Emax ** power) - (Emin ** power))
                                  ).value * u.GeV**2


                entry["Neutrino Energy"] = etot = (
                    window_sens * area * e_integral).to(u.erg)

                print gamma, sens, window_sens, dNdA, etot,

                # print window_sens, e_integral, etot, phi_power
                # raw_input("prompt")

                entry["pretty Neutrino Energy"] = etot/(u.erg)
                print entry["pretty Neutrino Energy"]

                raw_input("prompt")

                entry["CR Energy"] = etot / waxmann_bachall

                entry["Rad Energy"] = entry["CR Energy"] / baryon_loading

                entry["Average Neutrino Luminosity"] = etot / window_length

                new_row += [entry[x].value for x in specific]

            table[i] = np.array([tuple(new_row)], dtype=dt)
            i += 1

    # table = table[table['Time Integrated Sensitivity (gamma=2.0)'] > 0]

    for y in sorted(sens_dict.itervalues().next().keys()):
        key = "(" + str(y) + ")"

        to_print = ["Name", "Dec", "Time Integrated Sensitivity " + key,
                    "pretty Neutrino Energy " + key]
        headers = [r"Name", "Dec",
                   r"Time-Integrated Sensitivity $ \frac{1}{Tev times cm^{2}}$",
                   r"Required Neutrino Energy (erg)"]

        table = np.sort(
            table, order=['Time Integrated Sensitivity ' + key],
            axis=0).view()

        print tabulate(table[to_print], headers)

        title = r"Assumed Neutrino Spectrum $ \frac {d \phi}{dE} \propto " + \
                "E ^ {-" + str(y) + "} $"

        fields_1 = ["Name", "Time Integrated Sensitivity " + key]
        headers_1 = [
            r"Name", r"Time-Integrated Sensitivity [$ (\frac{TeV}{cm^{2}} $]"]

        export_table_to_tex.run(
            table[fields_1], headers_1, title, "sens_" + key)

        title = r"Assumed Neutrino Spectrum $ \frac {d \phi}{dE} \propto " + \
                "E ^ {-" + str(y) + "} $"

        fields_1 = ["Name", "Integrated Sensitivity " + key]
        headers_1 = [
            r"Name", r"Integrated Sensitivity $\frac{1}{cm^{2}} $"]

        export_table_to_tex.run(
            table[fields_1], headers_1, title, "integrated_sens_" + key)

        fields_2 = ["Name", 'pretty Neutrino Energy ' + key]

        headers_2 = [r"Name", r"Required Neutrino Energy (erg)"]

        table = np.sort(
            table, order=['Neutrino Energy ' + key], axis=0).view()

        export_table_to_tex.run(
            table[fields_2], headers_2, title, "Energy_" + key)
