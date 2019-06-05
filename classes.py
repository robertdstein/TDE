import zipfile
import json
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import os.path
import string
from ned import check_name, check_coordinate
from tabulate import tabulate
import lightcurves as lc
import plotting as p
import readswift
import readsupernova as rs
import selected_candidates as sc
import loglikelihoodminimisation as llh
import ConfigParser
from astropy import units as u
# from astropy import cosmology
# from astropy.cosmology import WMAP5, WMAP7
from astropy.coordinates import Distance
# cosmology.set_current(WMAP7)
import math

# Sets the minimum number of points required for a curve fit
# Number of measurements must at least equal number of model parameters

minimum_points = len(lc.return_histogram_labels())

# Sets the root directory containing scripts and config files/subdirectories
root = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"

conf = ConfigParser.ConfigParser()
config_path = root + "data_configs/full.ini"

def wrap_around_180(ra_deg):
    ra = np.deg2rad(ra_deg)
    ra[ra > np.pi] -= 2*np.pi
    return ra

def magnitude_to_luminosity(absolute_magnitude):
    """Convert an absolute magnitude to a distance.

    :param absolute_magnitude: Absolute Magnitude
    :return: luminosity
    """
    solar_luminosity = 3.826 * 10 ** 33
    solar_absolute_magnitude = 4.75
    exponent = (solar_absolute_magnitude - absolute_magnitude) / 2.5
    luminosity = solar_luminosity * 10 ** exponent
    return luminosity

class FullSet:
    """The full set of all TDE candidates
    """
    def __init__(self, path):
        self.path = path
        self.TDEs = container()
        self.total_entries = 0
        self.candidates_with_photometry = 0
        self.prime_candidates = 0
        self.all_available_bands = []
        self.best_list = []
        self.extract()
        self.update_time = datetime.datetime.now()
        self.creation_time = datetime.datetime.now()
        self.data_table = np.nan
        self.data_dict = np.nan
        self.neutrino_table = np.nan
        self.data_start = 54561.4746759
        self.list_catalogue()

        for cat in ["gold", "silver", "jetted", "obscured"]:
            self.export_catalogue_to_wikitext(cat)

    def extract(self):
        """Extract all data from tde.zip.
        The zipfile contains individual json files for each candidate.
        """
        zf = zipfile.ZipFile(self.path, 'r')
        self.update_time = datetime.datetime.now()
        self.creation_time = datetime.datetime(*zf.infolist()[0].date_time)

        file_list = zf.namelist()
        remove = ["tde-1980-2025-master/", 'tde-1980-2025-master/.gitignore',
                  'tde-1980-2025-master/---.json']
        file_list = [x for x in file_list if x not in remove]

        bn = set([])

        for name in file_list:
            with zf.open(name) as f:
                data = f.read()
                d = json.loads(data)
                info = d.values()

                self.total_entries += 1

                tde = TDE_Candidate(info)
                setattr(self.TDEs, tde.name, tde)
                if tde.has_photometry:
                    self.candidates_with_photometry += 1
                    for band_name in tde.bandlist:
                        bn.add(band_name)

        if hasattr(self.TDEs, "NGC247") and hasattr(self.TDEs, "NGC 247"):
            print "Removing duplicate entry of NGC 247"
            self.total_entries -= 1
            if self.TDEs.NGC247.has_photometry:
                self.candidates_with_photometry -= 1
            del self.TDEs.NGC247

        self.all_available_bands = list(bn)
        print "Dataset fully extracted. Of", self.total_entries, "we have", \
            self.candidates_with_photometry, "with at least", minimum_points, \
            "four photometric observations in one channel."
        # self.list_best()

    def creationtime(self):
        """Returns the time that the online catalogue was last updated
        """
        print "Catalogue version created:", self.creation_time, "\n"
        print "Dataset last expanded or updated on:", self.update_time, "\n"

    def info_candidate(self, name):
        """Plots all available light curves for a given candidate.

        :param name: Name of candidate
        """
        list_of_names = sorted(self.TDEs.__dict__.keys(), key=str.lower)

        print "\n"

        # Checks if name is valid. If not, prints possible names and requests
        # new input of name.

        while not hasattr(self.TDEs, name):
            original_name = str(name)
            testname = name.lower()
            testname = testname.replace(" ", "")
            testname = testname.replace("-", "")
            if "#" in name:
                no = int(name.replace("#", ""))
                name = list_of_names[no]
                print "Successfully matched candidate ID #"+str(no), \
                    "to candidate name", name
            else:
                matches = []
                for i, entry in enumerate(list_of_names):
                    matcher = entry.lower()
                    matcher = matcher.replace(" ", "")
                    matcher = matcher.replace("-", "")
                    if testname in matcher:
                        matches.append((i, entry))

                if len(matches) == 1:
                    name = matches[0][1]
                    print "Successfully matched", original_name, \
                        "to candidate name", name

                else:
                    if len(matches) == 0:
                        print "No TDE with this name was found!"
                        matches = list(enumerate(list_of_names))
                    else:
                        print "Multiple matches to this name were found. " \
                              "Please refine search."
                    print "Available candidates are:"
                    print matches
                    print "\n"
                    print "Please reenter candidate name. Candidate ID " \
                          "numbers can be entered in the form #XYZ"
                    name = raw_input("Candidate name or ID: \t")

        tde = getattr(self.TDEs, name)
        tde.fit_xray()
        tde.plot_candidate(fit=True)
        self.update_time = datetime.datetime.now()

    def plot_all_candidates(self, fit=False):
        """Plots all available light curves for each candidate.

        :param fit: Toggle True/False for whether to fit a model to the
        lightcurves, or simply plot them.
        """
        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            tde.fit_xray()
            tde.plot_candidate(fit)
            if fit:
                self.update_time = datetime.datetime.now()

    def plot_all_bands(self):
        """Plots all available light curves for each band, so that each
        candidate contributes one light curve."""
        for band in self.all_available_bands:
            title = band

            folder = "bands/"
            combined_folder = "combinedbands/"

            all_results = []
            combined_results = []

            for name in vars(self.TDEs):
                tde = getattr(self.TDEs, name)
                if tde.has_photometry:
                    if band in tde.bandlist:
                        i = tde.bandlist.index(band)

                        var = tde.varlist[i]
                        res = tde.return_photometry_dataset(band, var)

                        res["variable"] = var
                        res["label"] = name

                        all_results.append(res)
                        if (var == "magnitude") & (len(res["time"]) > 0):
                            combined_results.append(res)

            p.scatter_photometry(folder, title, all_results)
            p.scatter_photometry(combined_folder, title, combined_results,
                                 combined=True, tshift=True, sn=True)

    def list_catalogue(self):

        dt = np.dtype([
            ('Name', "S50"),
            ("Lowercase Name", "S50"),
            ("Alias", "S50"),
            ("Host", "S50"),
            ("Redshift", np.float),
            ("RA", np.float),
            ("Dec", np.float),
            ("Disc Date", np.float),
            ("Max Date", np.float),
            ("Last Upper Limit", np.float),
            ("Window Start Date", np.float),
            ("Window End Date", np.float),
            ("Window Length", np.float),
            ("NED RA", np.float),
            ("NED Dec", np.float),
            ("NED Redshift", np.float),
            ("Max Apparent Magnitude", np.float),
            ("Max Absolute Magnitude", np.float),
            ("Max Luminosity", np.float),
            ("Luminosity Distance (Mpc)", np.float),
            ("Category", "S50"),
            ("Weight", np.float)
        ])

        table = np.zeros(self.total_entries, dtype=dt)

        self.data_dict = dict()

        pre_emission_window = 365
        post_emission_window = 100

        min_pre_window = 30.

        for i, name in enumerate(vars(self.TDEs)):
            tde = getattr(self.TDEs, name)


            if not np.isnan(tde.lastul):
                window_start = float(tde.lastul)

                if (tde.mjdmax - window_start < min_pre_window) and (
                        tde.tde_category != "jetted"):
                    window_start = tde.mjdmax - min_pre_window

            else:
                window_start = tde.mjdmax - pre_emission_window

            if True:
                window_end = tde.mjdmax + post_emission_window

            print tde.redshift_float
            if np.isnan(float(tde.redshift_float)):
                lum_d = np.nan
            else:
                lum_d = Distance(z=float(tde.redshift_float)).to("Mpc").value

            try:
                weight = tde.weight.value
            except AttributeError:
                weight = np.nan

            # if str(tde.lumdist.u_value) == "Mpc":
            #     distance.append(tde.lumdist.value)
            # else:
            #     print vars(tde.lumdist)
            #     raise Exception("Units of distance are not correct!")

            table[i] = np.array([(
                str(tde.name), str(tde.name).lower(),
                str(tde.alias), str(tde.host),
                tde.redshift_float,
                tde.ra_deg, tde.dec_deg, tde.mjddisc, tde.mjdmax, tde.lastul,
                window_start, window_end, window_end - window_start,
                tde.NED_ra, tde.NED_dec, tde.NED_redshift,
                tde.app_mag_max, tde.abs_mag_max,  tde.lum, lum_d,
                tde.tde_category, weight
            )], dtype=dt)

            tde_dict = dict()
            tde_dict["ra_deg"] = tde.ra_deg
            tde_dict["dec_deg"] = tde.dec_deg
            tde_dict["mjdmax"] = tde.mjdmax
            tde_dict["mjddisc"] = tde.mjddisc
            tde_dict["lastul"] = tde.lastul
            tde_dict["start"] = window_start
            tde_dict["end"] = window_end
            tde_dict["z"] = tde.redshift_float
            tde_dict["dist"] = lum_d
            tde_dict["category"] = tde.tde_category
            tde_dict["season"] = []
            self.data_dict[tde.name] = tde_dict

        self.data_table = np.sort(
            table, order=['Lowercase Name'], axis=0).view()

        emission_window = pre_emission_window + post_emission_window

        mask = self.data_table["Window End Date"] < self.data_start

        headers = ["Name", "Category", "Redshift", "Luminosity Distance (Mpc)",
                   "RA", "Dec",
                   "Max Date", "Last Upper Limit",
                   "Window Start Date", "Window End Date", "Window Length",
                   "Weight"
                   ]

        print "Of", len(self.data_table["Max Date"]), "entries, we remove",
        print len(self.data_table['Window Start Date'][mask]),
        print "events that finished before IC40 (",
        print self.data_start, "), leaving",
        print len(self.data_table['Window Start Date'][~mask]), "events."

        print

        print "Rejected: \n"

        print tabulate(self.data_table[mask][headers], headers)

        print

        print "Accepted: \n"

        print tabulate(self.data_table[~mask][headers], headers)

        print

        self.data_table = np.sort(
            table, order=['Window Start Date'], axis=0).view()

        if os.path.isfile(config_path):
            conf.read(config_path)

            data_start = float(conf.get("IC40", "start_mjd"))

            pre_mask = self.data_table["Window Start Date"] < self.data_start

            for tde in self.data_table[pre_mask]["Name"]:
                self.data_dict[tde]["season"].append("Earlier")

            print "Early Entries \n"

            print "Before", self.data_start, "onwards, or unknown date."

            print "(", len(self.data_table[pre_mask]), ") entries \n"

            print tabulate(self.data_table[pre_mask][headers], headers)

            for season in conf.sections():
                print season, "\n"
                start_time = float(conf.get(season, "start_mjd"))
                end_time = float(conf.get(season, "end_mjd"))

                print start_time, "to", end_time,

                pre_mask = self.data_table["Window Start Date"] < end_time

                post_mask = self.data_table["Window End Date"] > start_time

                mix_mask = np.logical_and(pre_mask, post_mask)

                print "(", len(self.data_table[mix_mask]), ") entries \n"

                print tabulate(self.data_table[mix_mask][headers], headers)

                print

                for tde in self.data_table[mix_mask]["Name"]:
                    self.data_dict[tde]["season"].append(season)

            # Does lightcurve peak close enough to the end of data, or after
            # the end of data, meaning that emission extends beyond 7 year?

            other_mask = (self.data_table["Window End Date"]) < end_time

            print "Remaining Entries \n"

            print "From", end_time, "onwards, or unknown date."

            print tabulate(self.data_table[~other_mask][headers], headers)

            print

            for tde in self.data_table[~other_mask]["Name"]:
                self.data_dict[tde]["season"].append("Later")

        # print self.data_table["Category"]

        for category in set(self.data_table["Category"]):
            print "Category", category

            mask = self.data_table["Category"] == category

            print tabulate(self.data_table[mask][headers], headers)

    def export_catalogue_to_wikitext(self, category):
        path = "wikitext_" + category + ".txt"

        mask = self.data_table["Window Start Date"] < self.data_start

        dt = self.data_table[~mask]

        ct = dt[dt["Category"] == category]

        with open(path, "w") as f:

            variables = ["Name",  "Luminosity Distance (Mpc)",
                         "RA", "Dec", "Window Start Date",
                         "Max Date", "Window End Date",
                         "Window Length", "Category"]

            f.write("{| class ='wikitable sortable' border='1'|thumb| \n")
            headers = "!"
            for var in variables:
                headers += " '''" + var + "'''||"

            f.write(headers[:-2] + " \n")

            for row in ct[variables]:
                f.write("|- \n")
                new = "|"
                for entry in row:
                    new += "  " + str(entry) + "  ||"

                # print row, new
                f.write(new[:-2] + " \n")

            f.write("|}")

    def list_best(self):
        """Lists the TDE candidates which are most promising. The criteria
        for this are set in the IsBest() function for the TDE class.
        """
        best = []
        ratio = []
        print "\n The best candidates, in terms of available data, are: \n"
        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            if tde.best:
                best.append(name)
                print name
            if tde.ratiotest:
                ratio.append(name)
        self.best_list = best
        print "\n \n"
        print "Those passing the ratio test are: \n"
        print ratio

    def list_fits(self):
        """Prints each candidate that has been fitted"""
        print "The following candidates have been fitted: \n"
        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            if len(tde.fits) > 0:
                print tde.name, ":"
                for d in tde.fits:
                    print "\t", d, "\t", list(tde.fits[d])

    def plot_general_properties(self, names=None, suffix=""):
        """ Plots binned distributions of candidates for each variable in
        'variables' list.
        """
        variables = ["redshift", "comovingdist", "lumdist", "hostoffsetdist",
                     "mjddisc", "nbands", "maxabsmag", "ra_float",
                     "dec_float", "ebv"]

        titles = ["Redshift", "Comoving Distance", "Luminosity Distance",
                  "Host Offset Distance", "Date of Discovery",
                  "Number of observational bands", "Max Absolute Magnitude",
                  "Right ascension", "Declination", "EBV"]

        x_labels = ["Redshift (z)", "Distance (MPc)", "Distance (MPc)",
                    "Angular Diameter Distance", "Time(MJD)", "n", "", "hours",
                    "degrees", "Magnitudes"]

        if names is None:
            names = vars(self.TDEs)

        values = []
        for i in range(0, len(variables)):
            values.append([])

        for name in names:
            tde = getattr(self.TDEs, name)
            for i in range(0, len(variables)):
                var = variables[i]
                if hasattr(tde, var):
                    val = getattr(tde, var)
                    if not isinstance(val, float):
                        val = getattr(tde, var).value

                    if not np.isnan(float(val)):
                        values[i].append(float(val))

        print "Plotting histograms for the following:", titles
        p.histograms_plot(titles, x_labels, values, suffix=suffix)

    def plot_catalogue_general_properties(self):
        all_catalogues = sc.return_all_catalogues()
        for key, dictionary in all_catalogues.iteritems():
            suffix = "_" + key
            self.plot_general_properties(dictionary["names"], suffix)

    def scatter_data_2d(self, title, variable_config, dictionary=sc.full_dict(),
                        loose=False):
        """Plots scatter distributions of different variables

        :param title: Title of the Scatter Graph
        :param variables: List containing the attribute-name of variable to plot
        :param variable_names: List containing the names of the variable
        :param pairIDs: List of variable ID pairs in the form [x,y],
        each corresponding to a subplot scatter graph in which variable
        number x in variables is plotted against variable number y
        :param dictionary: Contains a list of the names of candidates to be
        included. (If None, the default full list will be used.) Also contains
        a list of bands to be plotted. (If None, all will be included).
        :param loose: Toggle True/False to swap between loose and restricted
        fitting bounds
        """
        folder = "misc/"

        names = dictionary["names"]
        bands = dictionary["bands"]

        if names is None:
            names = vars(self.TDEs)

        variables = variable_config["variables"]
        variable_names = variable_config["variable_names"]

        suffixes = [""]

        if loose:
            suffixes.append("_loose")

        for suffix in suffixes:

            allvals = []
            for pair in variable_config["pairIDs"]:
                res = dict()

                letters = ["x", "y", "z"]
                allvars = []

                for i, ID in enumerate(pair):
                    letter = letters[i]
                    res[letter + "var"] = variables[ID]
                    res[letter + "_label"] = variable_names[ID]
                    allvars.append(variables[ID])

                res["vars"] = allvars

                allvals.append(res)

            for res in allvals:
                x = []
                y = []
                z = []
                for name in names:
                    tde = getattr(self.TDEs, name)
                    if tde.has_photometry:
                        vals = []
                        for var in res["vars"]:
                            path = (var+suffix).split(".")
                            val = tde
                            while len(path) > 0:
                                if (path[0] == "*") & (isinstance(val, dict)):
                                    data = []
                                    for key in val.keys():
                                        d = val[key]
                                        if (bands is not None) and \
                                                (key not in bands):
                                            pass
                                        else:
                                            if isinstance(d, list):
                                                data.extend(d[path[1]])
                                            else:
                                                data.append(d[path[1]])
                                    val = data
                                    break

                                try:
                                    val = val[path[0]]
                                except AttributeError:
                                    if hasattr(val, path[0]):
                                        val = getattr(val, path[0])
                                    else:
                                        break
                                except KeyError:
                                    break

                                if isinstance(val, unicode):
                                    val = float(val)

                                path = path[1:]

                            if isinstance(val, float) or isinstance(val, list):
                                vals.append(val)

                        if len(vals) == len(res["vars"]):
                            for k, v in enumerate(vals):
                                if isinstance(v, list):
                                    [x, y, z][k].extend(v)
                                elif isinstance(v, float):
                                    [x, y, z][k].append(v)
                res["x"] = x
                res["y"] = y
                if len(z) > 0:
                    res["z"] = z

            p.scatter_distribution(folder, title + suffix, allvals)

    def plot_2d_distributions(self):
        """Plots a 2d distribution for general candidate parameters"""
        title = "2D_distributions"

        variable_config = p.config_for_general_scatter()

        all_dictionaries = sc.return_all_catalogues()

        for name, dictionary in all_dictionaries.iteritems():
            self.scatter_data_2d(title + "_" + name, variable_config,
                                 dictionary)


    def plot_2d_fit_distributions(self, dictionary=sc.full_dict(),
                                  title="2D_fit_distributions"):
        """Plots 2D distributions for parameters of the fit.

        :param dictionary: Contains a list of the names of candidates to be
        included. (If None, the default full list will be used.) Also contains
        a list of bands to be plotted. (If None, all will be included).
        :param title: Title of 2D scatter plots
        """
        variable_config = p.config_for_fit_scatter()

        all_dictionaries = sc.return_all_catalogues()

        for name, dictionary in all_dictionaries.iteritems():
            self.scatter_data_2d(title + "_" + name, variable_config,
                                 dictionary, loose=True)

    def plot_2d_fit_xray(self):
        bands, varnames = sc.xray_variables()
        title = "2D_fit_distribution_xray"
        self.plot_2d_fit_distributions(names=None, bands=bands, title=title)

    def plot_fit_parameters(self, names=None, bands=None,
                            savepath="graphs/optical_fits/",
                            root="combined/"):
        """Plots binned distributions of candidates for each variable
        in 'variables'.

        :param names: Names of candidates to be included. If default of None,
        uses all available TDEs.
        :param bands: Selects the bands to be included
        :param savepath: Path for saving the graphs
        :param root: Root of folder for combined histograms
        """

        all_dictionaries = sc.return_all_catalogues()

        for key, dictionary in all_dictionaries.iteritems():

            full = [[], []]
            titles, xlabels = lc.return_all_parameter_names()

            names = dictionary["names"]
            if names is None:
                names = vars(self.TDEs)

            for name in names:
                tde = getattr(self.TDEs, name)
                if len(tde.fits) > 0:
                    path = savepath + tde.name + "/"
                    if not os.path.exists(path):
                        os.mkdir(path)
                    data = tde.fit_parameter_histogram(path, bands)

                    if len(data[0][0]) == 0:
                        os.rmdir(path)

                    for i, subset in enumerate(data):
                        for j in range(len(subset)):
                            if len(full[i]) < len(subset):
                                full[i].append([])
                            full[i][j].extend(subset[j])

            for i, name in enumerate(["", "_loose"]):
                print "\n"
                print name
                print "\n"
                print "Parameters \t Max \t \t Min \n"

                for j, param in enumerate(xlabels):
                    print param, \
                        "\t \t", "{0:.3}".format(float(min(full[i][j]))),\
                        "\t \t", "{0:.3}".format(float(max(full[i][j])))
                print "\n"
                p.histograms_plot(
                    xlabels, xlabels, full[i],
                    savepath + root + key + name + "_")

    def plot_fit_parameters_xrays(self):
        xrbs, xvars = sc.xray_variables()
        savepath = "graphs/xrays/"
        root = "all"
        self.plot_fit_parameters(bands=xrbs, savepath=savepath, root=root)

    def plot_spectra(self):
        """Plots a spectrum for each candidate with spectral data. """
        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            if tde.has_spectra:
                tde.plot_spectrum()

    def plot_skymap(self):
        """Plots a map of the distribution of TDE candidates on the sky.
        Uses the redshift as a colour scale.
        """
        ra = []
        dec = []
        rs = []
        markers = []

        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            if (hasattr(tde, "ra_deg")) and (hasattr(tde, "redshift")):
                ra.append(float(tde.ra_deg))
                dec.append(float(tde.dec_deg))
                rs.append(float(tde.redshift.value))
                if tde.best:
                    markers.append("*")
                else:
                    markers.append("o")

        plt.figure()
        plt.subplot(111, projection="aitoff")

        cm = plt.cm.get_cmap('RdYlGn_r')
        sc = plt.scatter(
            wrap_around_180(ra), np.deg2rad(dec), c=rs, s=35, cmap=cm)
        cbar = plt.colorbar(sc)
        cbar.set_label('Redshift(z)')

        path = "graphs/misc/skymap.pdf"
        print "Saving to", path
        plt.tight_layout()
        plt.savefig(path)
        plt.close()

    def add_supernova_fits(self):
        self.sn_fits = rs.plot()

    def plot_supernova(self):
        if not hasattr(self, "sn_fits"):
            self.add_supernova_fits()
        rs.fit_scatter(self.sn_fits)
        rs.fit_histogram(self.sn_fits)

    def add_swift(self):
        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            tde.swiftdata()


class container:
    """A container for new data
    """
    def __init__(self):
        pass


class TDE_Candidate:
    """A TDE Candidate from the catalogue
    """
    def __init__(self, jsonfile):
        self.raw = jsonfile[0]
        self.name = jsonfile[0]["name"]
        self.has_photometry = False
        self.has_spectra = False
        self.ratiotest = False
        self.has_xray = False
        self.type = "TDE"
        self.best = False
        self.fits = dict()

        self.alias = []

        self.mjddisc = np.nan
        self.mjdmax = np.nan
        self.lastul = np.nan
        self.ra_deg = np.nan
        self.dec_deg = np.nan
        self.redshift_float = np.nan
        self.app_mag_max = np.nan
        self.abs_mag_max = np.nan
        self.lum = np.nan

        self.NED_ra = np.nan
        self.NED_dec = np.nan
        self.NED_redshift = np.nan
        self.NED_name = np.nan

        self.tde_category = np.nan

        for key, val in jsonfile[0].items():

            if key == "alias":
                for entry in val:
                    self.alias.append(str(entry["value"]))

            elif key == "host":
                pass

            elif key != "photometry":
                if isinstance(val, list):

                    # if len(val) > 1:
                    #     print self.name, key

                    con = container()
                    for key2, val2 in val[0].items():
                        setattr(con, key2, val2)
                    setattr(self, key, con)
                else:
                    setattr(self, key, val)
            else:
                self.group_photometry(val)

        # In the Open TDE Catalog, ra/dec are only provided if different from
        # host. Thus, sets ra/dec to equal host ra/dec if these are not given.

        if not (hasattr(self, "ra")) and not (hasattr(self, "dec")):
            if hasattr(self, "hostra") and hasattr(self, "hostdec"):
                self.ra = self.hostra
                self.dec = self.hostdec

        self.coordinates()

        accept = False

        to_check = list(self.alias)

        if "host" in jsonfile[0].keys():
            val = jsonfile[0]["host"]
            for host in val:
                to_check.append(host["value"])

        hosts = []
        NED_hosts = []
        for name in to_check:

            print "Checking", name

            entry = check_name(name)
            if entry is not None:
                if len(entry["data_table"][entry["mask"]]) > 0:
                    NED_hosts.append(entry)
                    hosts.append(entry["data_table"][entry["mask"]][
                                     "Object Name"])
                # else:
                #     print entry["data_table"]
                    # raw_input("prompt")

        if len(NED_hosts) == 1:

            NED_galaxy = NED_hosts[0]["data_table"][0]
            hosts = str(list(hosts[0])[0])

        elif len(NED_hosts) == 0:

            entry = check_coordinate(self.ra_deg, self.dec_deg)

            if entry is not None:

                val = jsonfile[0]["alias"]

                candidate = entry["data_table"][entry["mask"]][0]

                new_name = string.replace(candidate["Object Name"],
                                          " ", "")


                try:

                    potential_hosts = [
                        string.replace(x["value"], " ", "") for x in val]

                    alt = string.replace(new_name, "GALEXASC", "GALEX")
                    # print self.name, self.alias, potential_hosts, alt
                    #
                    # print self.name, self.alias, potential_hosts, new_name

                    if new_name in potential_hosts:
                        hosts = str(val[potential_hosts.index(new_name)]["value"])
                        accept = True

                    elif alt in potential_hosts:
                        accept = True
                        hosts = str(val[potential_hosts.index(alt)]["value"])

                    if accept:
                        print "Match!"
                        NED_galaxy = candidate

                except KeyError:
                    pass
        else:

            identical = True

            for host in NED_hosts:
                name = str(host["data_table"]["Object Name"])
                for alt_host in NED_hosts:
                    alt_name = str(alt_host["data_table"]["Object Name"])
                    if name != alt_name:
                        identical = False


            if identical:
                NED_galaxy = NED_hosts[0]["data_table"][0]
                print hosts,
                hosts = str(list(hosts[0])[0])
                print hosts

            else:
                raise Exception("Found too many hosts!")

        # else:
        #     print "No host or", self.name
        #     for name in self.alias:
        #         entry = check_name(name)
        #         if entry is not None:
        #             print entry["data_table"]
        #
        #     print self.name, self.alias
        try:
            self.NED_name = str(NED_galaxy["Object Name"])
            self.NED_ra = float(NED_galaxy["RA(deg)"])
            self.NED_dec = float(NED_galaxy["DEC(deg)"])
            self.NED_redshift = float(NED_galaxy["Redshift"])

            self.host = hosts

        except UnboundLocalError:
            self.host = np.nan

        if hasattr(self, "redshift"):
            self.redshift_float = self.redshift.value

        if hasattr(self, "maxappmag"):
            self.app_mag_max = self.maxappmag.value

        if hasattr(self, "maxabsmag"):
            self.abs_mag_max = self.maxabsmag.value

        if hasattr(self, "category"):
            self.tde_category = self.category.value

        self.swiftdata()
        self.spectrum()

        self.mjd_time()
        self.is_best()

    def add_ned(self):
        best = None

        if hasattr(self, "host"):
            print self.name, self.host.value, vars(self.host)

            # entry = check_name(str(self.host.value))
            # if entry is not None:
            #     best = entry["data_table"]
            #     # if best is None:
            #     #
            #     # else:
            #     #     if best["Object Name"] == \
            #     #             entry["data_table"]["Object Name"]:
            #     #         pass
            #     #     else:
            #     #         pass
            #             # print self.alias
            #             # print best
            #             # print entry["data_table"]
            #             # raise Exception("Conflict with aliases. Each matches "
            #             #                 "to a different object!")
            #
            #     if len(best) > 1:
            #         raise Exception("Too many entries")



        # entry = check_coordinate(self.ra.deg, self.dec.deg)
        # if entry is not None:
        #     try:
        #         best = entry["data_table"][entry["mask"]][0]
        #         # if best is not None:
        #         #     if best["Object Name"] != new["Object Name"]:
        #         #         name_mask = entry["data_table"]["Object Name"] == \
        #         #                     best["Object Name"]
        #         #
        #         #         alt = entry["data_table"][name_mask]
        #         #
        #         #         if len(alt) == 1:
        #         #             best = alt
        #         #         else:
        #         #             check =[
        #         #             x for x in entry["data_table"]["Object Name"]
        #         #             if any(y.replace(" ", "") in x.replace(" ", "")
        #         #             for y in self.alias)]
        #         #
        #         #             if len(check) > 0:
        #         #
        #         #                 best = entry["data_table"][
        #         #                     entry["data_table"]["Object Name"] ==
        #         #                     check[0]]
        #         #             else:
        #         #                 best = new
        #         #     else:
        #         #         best = new
        #         # else:
        #         #     best = new
        #     except IndexError:
        #         pass

        if best is not None:
            self.alias.append(best["Object Name"])
            self.NED_name = str(best["Object Name"])
            self.NED_ra = float(best["RA(deg)"])
            self.NED_dec = float(best["DEC(deg)"])
            self.NED_redshift = float(best["Redshift"])
            self.NED_coords = SkyCoord(
                str(self.NED_ra) + " " + str(self.NED_dec),
                unit=(u.deg, u.deg))
            self.NED_offset = self.NED_coords.separation(self.coords).arcsec

    def swiftdata(self):
        """Reads datasets downloaded from SWIFT in the 0.3-10KeV range,
        saved in location swiftdata/NAME.qdp .
        """
        datapath = "swiftdata/" + self.name + ".qdp"
        if os.path.isfile(datapath):

            print "Found file:", datapath,

            band, metadata, entries = readswift.run(datapath, self.name)
            newentries = readswift.condense(entries)

            variable = "countrate"

            self.swiftmetadata = metadata
            self.swift_raw = entries
            self.bandlist.append(band)
            self.varlist.append(variable)
            self.nbands += 1
            setattr(self.photometry, band + "_" + variable, newentries)

    def coordinates(self):
        """Converts the coordinates given in the json files to floats.
        """
        if hasattr(self, "ra") and hasattr(self, "dec"):
            self.ra_float = get_degrees(self.ra.value)
            self.dec_float = get_degrees(self.dec.value)
            c = SkyCoord(
                self.ra.value, self.dec.value, unit=(u.hourangle, u.deg))
            self.ra_deg = c.ra.degree
            self.dec_deg = c.dec.degree

    def is_best(self):
        """Assesses whether a given candidate is promising.
        """
        # self.fit_xray()
        if self.has_photometry:
            if hasattr(self, "mjdmax") and hasattr(self, "mjdvismax"):
                if self.nbands > 2.0:
                    self.best = True

    def mjd_time(self):
        """Converts the discovery date to MJD units.
        Does the same for the maximum date, and maximum visual date.
        """

        print self.name

        if hasattr(self, "discoverdate"):
            val = self.discoverdate.value.split("/")
            if len(val) == 3:
                [y, m, d] = val
                d = d.split(".")[0]
                discdate = y + "-" + m + "-" + d + "T00:00:00"
                t = Time(discdate)
                setattr(self, "mjddisc", float(t.mjd))

        if hasattr(self, "maxdate"):
            [y, m, d] = self.maxdate.value.split("/")
            maxdate = y + "-" + m + "-" + d + "T00:00:00"
            t = Time(maxdate)
            setattr(self, "mjdmax", float(t.mjd))

        if hasattr(self, "maxvisualdate"):
            [y, m, d] = self.maxvisualdate.value.split("/")
            maxdate = y + "-" + m + "-" + d + "T00:00:00"
            t = Time(maxdate)
            setattr(self, "mjdvismax", float(t.mjd))

        if hasattr(self, "lastupperlim"):
            [y, m, d] = self.lastupperlim.value.split("/")
            uldate = y + "-" + m + "-" + d + "T00:00:00"
            t = Time(uldate)
            setattr(self, "lastul", float(t.mjd))

    def group_photometry(self, val):
        """Extracts the photometry data from the json file, and groups it
        into separate classes for each available band.
        """
        allbands = container()
        bandlist = []
        varlist = []
        for entry in val:
            if "band" in entry.keys():
                band = entry["band"]
                variable = "magnitude"

            elif "energy" in entry.keys():
                band = "X-Ray (" + \
                       str(entry["energy"][0]) + "-" + str(entry["energy"][1]) \
                       + ") " + str(entry["u_energy"])
                variable = "luminosity"
            elif "frequency" in entry.keys():
                band = "Gamma Ray"
                variable = "fluxdensity"
            else:
                print entry
                raise Exception("No band info!")

            if variable in entry.keys():
                pass
            elif "countrate" in entry.keys():
                variable = "countrate"
            elif "flux" in entry.keys():
                variable = "flux"
            else:
                print entry
                raise Exception("No suitable variable found!")

            if hasattr(allbands, band + "_" + variable):
                pass
            else:
                setattr(allbands, band + "_" + variable, [])
                bandlist.append(str(band))
                varlist.append(str(variable))

            getattr(allbands, band + "_" + variable).append(entry)

        self.bandlist = bandlist
        self.varlist = varlist
        self.photometry = allbands
        self.nbands = float(len(bandlist))
        for band in vars(allbands):
            if len(getattr(self.photometry, band)) > minimum_points:
                self.has_photometry = True

    def fit_xray(self):
        """Checks to see if there's a light curve visible in X Ray measurements.
        Fits this distribution to see if it is close to the characteristic
        t^-5/3.
        """
        xrbands, varnames = sc.xray_variables()
        path = "graphs/xrays/" + self.name + "/"
        if hasattr(self, "bandlist"):
            for (band, var) in zip(xrbands, varnames):
                if hasattr(self.photometry, band + "_" + var):
                    fitdata = self.dataset_to_fit(band, var)
                    if "t" in fitdata.keys():
                        npoints = len(fitdata["t"])
                        if npoints > 2:
                            if "upt" in fitdata.keys():
                                npoints += len(fitdata["upt"])
                            if npoints > minimum_points:
                                fig = plt.figure()

                                print "Minimising for candidate", self.name, \
                                    "in band", band

                                info = llh.run(fitdata, self.name)
                                self.fits[band] = info

                                fig.set_size_inches(10, 10)

                                if not os.path.exists(path):
                                    os.mkdir(path)
                                file_name = path + band + ".pdf"
                                plt.savefig(file_name)
                                plt.close()
                                print "Saving to", file_name

    def fit_optical(self):
        """Checks to see if there is a light curve visible in X Ray
        measurements. Fits this distribution to see if it is close to the
        characteristic t^-5/3.
        """
        var = "magnitude"
        path = "graphs/optical_fits/" + self.name + "/"
        if hasattr(self, "bandlist"):
            for band in self.bandlist:
                if hasattr(self.photometry, band + "_" + var):
                    fitdata = self.dataset_to_fit(band, var)
                    if "t" in fitdata.keys():
                        npoints = len(fitdata["t"])
                        if npoints > 2:
                            if "upt" in fitdata.keys():
                                npoints += len(fitdata["upt"])
                            if npoints > minimum_points:
                                fig = plt.figure()

                                print "Minimising for candidate", self.name, \
                                    "in band", band

                                info = llh.run(fitdata, self.name)
                                self.fits[band] = info

                                fig.set_size_inches(7, 10)

                                if not os.path.exists(path):
                                    os.mkdir(path)
                                file_name = path + band + ".pdf"
                                plt.savefig(file_name)
                                plt.close()
                                print "Saving to", file_name

        # filename = path + "histogram.pdf"
        if len(self.fits) > 0:
            self.fit_parameter_histogram(path)

    def fit_parameter_histogram(self, path, bands=None):
        """Produces histograms to show the distribution of fit parameters
        """
        alllists = [[], []]

        for i, suffix in enumerate(["", "_loose"]):
            data = self.fits

            if bands is None:
                bands = data.keys()
                xrbs, xrvars = sc.xray_variables()
                bands = [x for x in bands if x not in xrbs]

            titles, xlabels = lc.return_all_parameter_names()

            for j in range(len(xlabels)):
                alllists[i].append([])

            for band in bands:
                if band in data.keys():
                    entry = data[band]
                    for j, var in enumerate(xlabels):
                        alllists[i][j].append(entry[var+suffix])

            if len(alllists[i][0]) > 0:
                p.histograms_plot(titles, xlabels, alllists[i], path,
                                  suffix=suffix)
        return alllists

    def dataset_to_fit(self, band, var):
        """Checks to see if there is a light curve visible in swift X Ray
        measurements. Fits this distribution to see if it is close to the
        characteristic t^-5/3.
        """
        fitdata = dict()
        fitdata["var"] = var
        fitdata["band"] = band

        if hasattr(self, "bandlist"):
            if band in self.bandlist:

                results = self.return_photometry_dataset(band, var)

                if len(results["y"]) > 0:
                    fitdata["t"] = []
                    fitdata["y"] = []
                    fitdata["up"] = []
                    fitdata["down"] = []
                    fitdata["ul_t"] = results["upper_limit_time"]
                    fitdata["ul_y"] = results["upper_limit_y"]

                    for i in range(len(results["time"])):
                        t = results["time"][i]
                        fitdata["t"].append(float(t))
                        fitdata["y"].append(float(results["y"][i]))
                        if results["sigma_upper"][i] > 0:
                            fitdata["up"].append(results["sigma_upper"][i])
                            fitdata["down"].append(results["sigma_lower"][i])
                        else:
                            fitdata["up"].append(float(results["y"][i])*0.5)
                            fitdata["down"].append(float(results["y"][i])*0.5)

                    if hasattr(self, "mjdvismax"):
                        fitdata["maxtime"] = self.mjdvismax
                    else:
                        my = max(fitdata["y"])
                        index = fitdata["y"].index(my)
                        fitdata["maxtime"] = fitdata["t"][index]

        return fitdata

    def spectrum(self):
        """Checks to see if the candidate has spectral data.
        """
        if hasattr(self, "spectra"):
            if self.spectra.u_fluxes != "Uncalibrated":
                self.has_spectra = True

    def plot_spectrum(self):
        """Plots the spectral data for a given candidate.
        """
        if self.has_spectra:
            plt.figure()

            wl = []
            f = []

            spec = self.spectra

            for entry in spec.data:
                f.append(entry[0])
                wl.append(entry[1])

            plt.plot(f, wl)

            plt.xlabel(spec.u_wavelengths)
            plt.ylabel(spec.u_fluxes)

            title = self.name
            if hasattr(spec, "redshift"):
                title += " [z =" + str(spec.redshift) + "]"
            if hasattr(self, "mjdmax") & hasattr(spec, "time"):
                time_after_max = int(float(spec.time) - float(self.mjdmax))
                if time_after_max < 0:
                    title += " [" + str(-1 * time_after_max) + " days before " \
                                                              "max]"
                else:
                    title += " [" + str(time_after_max) + " days after max]"

            plt.title(title)

            path = "graphs/spectra/" + self.name + ".pdf"
            print "Saving to", path
            plt.savefig(path)
            plt.close()

    def return_photometry_dataset(self, band, var):
        if hasattr(self.photometry, band + "_" + var):
            data = getattr(self.photometry, band + "_" + var)

            results = dict()
            results["time"] = []
            results["y"] = []
            results["upper_limit_time"] = []
            results["upper_limit_y"] = []
            results["sigma_upper"] = []
            results["sigma_lower"] = []

            if hasattr(self, "mjdmax"):
                results["maxdate"] = [self.mjdmax]
            else:
                results["maxdate"] = [None]

            if hasattr(self, "mjdvismax"):
                results["maxdate"].append(self.mjdvismax)
            else:
                results["maxdate"].append(None)

            for entry in data:
                if "time" in entry.keys():
                    t = entry["time"]
                    if isinstance(t, list):
                        t = t[0]
                    if "upperlimit" in entry.keys():
                        results["upper_limit_y"].append(float(entry[var]))
                        results["upper_limit_time"].append(float(t))
                    else:
                        results["y"].append(float(entry[var]))
                        results["time"].append(float(t))

                        if "e_upper_" + var in entry.keys():
                            sigup = float(entry["e_upper_" + var])
                            sigdown = float(entry["e_lower_" + var])
                        elif "e_" + var in entry.keys():
                            sigup = float(entry["e_" + var])
                            sigdown = sigup
                        else:
                            print "Can't find value for error in " + var
                            sigup = 0.0
                            sigdown = sigup
                        results["sigma_upper"].append(sigup)
                        results["sigma_lower"].append(sigdown)

            return results
        else:
            return None

    def plot_candidate(self, fit=False):
        if fit:
            self.fit_optical()
        elif self.has_photometry:
            folder = "candidates/"
            opticalfolder = "optical/"
            title = self.name

            allresults = []
            opticalresults = []

            for band in sorted(self.bandlist, key=str.lower):
                i = self.bandlist.index(band)

                var = self.varlist[i]
                res = self.return_photometry_dataset(band, var)

                res["variable"] = var
                res["label"] = band + "_" + var
                res["band"] = band

                allresults.append(res)
                if var == "magnitude":
                    opticalresults.append(res)

            p.scatter_photometry(folder, title, allresults)
            if len(opticalresults) > 0:
                p.scatter_photometry(opticalfolder, title, opticalresults,
                                     combined=True)
        else:
            print "No photometry was found for", self.name


def get_degrees(time_str):
    """Converts a HH:MM:SS time string to a float (number of hours).
    """
    if time_str.count(':') > 1.:
        s, arcm, arcs = time_str.split(':')
        if float(s) < 0.0:
            return -(abs(float(s)) + (float(arcm) / 60.) + (float(arcs) /
                                                            3600.))
        else:
            return float(s) + (float(arcm) / 60.)+(float(arcs) / 3600.)
    else:
        arcm, arcs = time_str.split(':')
        if float(arcm) < 0.0:
            return -(abs(float(arcm) / 60.)+(float(arcs) / 3600.))
        else:
            return (float(arcm) / 60.) + (float(arcs) / 3600.)
