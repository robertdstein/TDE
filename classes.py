import zipfile
import math
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

import lightcurves as lc
import plotting as p
import readswift
import readsupernova as rs
import selected_candidates as sc
import loglikelihoodminimisation as llh

minpoints = len(lc.return_histogram_labels())


def wrapAround180(ra_deg):
    ra = np.deg2rad(ra_deg)
    ra[ra > np.pi] -= 2*np.pi
    return ra


def magtolum(absmag):
    """Convert an absolute magnitude to a distance.
    """
    solarlum = 3.826 * 10**33
    solarabsmag = 4.75
    exponent = (solarabsmag - absmag)/2.5
    lum = solarlum * 10 ** exponent
    return lum


class FullSet:
    """The full set of all TDE candidates
    """
    def __init__(self, path):
        self.path = path
        self.TDEs = container()
        self.totalentries = 0
        self.candidates_with_photometry = 0
        self.prime_candidates = 0
        self.all_available_bands = []
        self.best_list = []
        self.extract()

    def extract(self):
        """Extract all data from tde.zip.
        The zipfile contains individual json files for each candidate.
        """
        zf = zipfile.ZipFile(self.path, 'r')
        self.update_time = datetime.datetime.now()
        self.creation_time = datetime.datetime(*zf.infolist()[0].date_time)

        filelist = zf.namelist()[2:]

        bn = set([])
        for name in filelist:
            with zf.open(name) as f:
                data = f.read()
                d = json.loads(data)
                info = d.values()

                self.totalentries += 1

                tde = TDE_Candidate(info)
                setattr(self.TDEs, tde.name, tde)
                if tde.has_photometry:
                    self.candidates_with_photometry += 1
                    for bandname in tde.bandlist:
                        bn.add(bandname)

        self.all_available_bands = list(bn)
        print "Dataset fully extracted. Of", self.totalentries, "we have", \
            self.candidates_with_photometry, "with at least", minpoints, \
            "four photometric observations in one channel."
        self.list_best()

    def creationtime(self):
        """Returns the time that the online catalogue was last updated
        """
        print "Catalogue version created:", self.creation_time, "\n"
        print "Dataset last expanded or updated on:", self.update_time, "\n"

    def info_candidate(self, name):
        """Plots all available light curves for a given candidate.
        """
        list_of_names = sorted(self.TDEs.__dict__.keys(), key=str.lower)

        print "\n"

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
        """
        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            tde.fit_xray()
            tde.plot_candidate(fit)
            if fit:
                self.update_time = datetime.datetime.now()

    def plot_all_bands(self):
        """Plots all available light curves for each band.
        """
        for band in self.all_available_bands:
            title = band

            folder = "bands/"
            combifolder = "combinedbands/"

            allresults = []
            combiresults = []

            for name in vars(self.TDEs):
                tde = getattr(self.TDEs, name)
                if tde.has_photometry:
                    if band in tde.bandlist:
                        i = tde.bandlist.index(band)

                        var = tde.varlist[i]
                        res = tde.return_photometry_dataset(band, var)

                        res["variable"] = var
                        res["label"] = name

                        allresults.append(res)
                        if (var == "magnitude") & (len(res["time"]) > 0):
                            combiresults.append(res)

            p.scatter_photometry(folder, title, allresults)
            p.scatter_photometry(combifolder, title, combiresults,
                                 combined=True, tshift=True, sn=True)

    def list_best(self):
        """Lists the TDE candidates which are most promising.
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
        print "The following candidates have been fitted: \n"
        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            if len(tde.fits) > 0:
                print tde.name, ":"
                for d in tde.fits:
                    print "\t", d, "\t", list(tde.fits[d])

    def plot_distibutions(self):
        """ Plots binned distributions of candidates for each variable in 'variables'.
        """
        variables = ["redshift", "comovingdist", "lumdist", "hostoffsetdist",
                     "mjdmax", "nbands","maxabsmag", "ra_float",
                     "dec_float", "ebv"]

        titles = ["Redshift", "Comoving Distance", "Luminosity Distance",
                  "Host Offset Distance", "Date of Max Luminosity",
                  "Number of observational bands", "Max Absolute Magnitude",
                  "Right ascension", "Declination", "EBV"]

        xlabels = ["Redshift (z)", "Distance (MPc)", "Distance (MPc)",
                   "Angular Diameter Distance", "Time(MJD)", "n", "", "hours",
                   "degrees", "Magnitudes"]

        values = []
        for i in range (0, len(variables)):
            values.append([])

        for name in vars(self.TDEs):
            tde = getattr(self.TDEs, name)
            for i in range(0, len(variables)):
                var = variables[i]
                if hasattr(tde, var):
                    val = getattr(tde, var)
                    if not isinstance(val, float):
                        val = getattr(tde, var).value
                    values[i].append(float(val))

        print "Plotting histograms for the following:", titles
        p.histograms_plot(titles, xlabels, values)

    def scatter_data_2D(self, title, variables, variable_names, pairIDs,
                        names=None, bands=None, loose=False):
        """Plots scatter distributions of different variables
        """
        folder = "misc/"

        if names is None:
            names = vars(self.TDEs)

        suffixes = [""]

        if loose:
            suffixes.append("_loose")

        for suffix in suffixes:

            allvals = []
            for pair in pairIDs:
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
                                        d=val[key]
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

                            if  isinstance(val, float) or isinstance(val, list):
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

            p.scatter_distibution(folder, title+suffix, allvals)

    def plot_2d_distributions(self):
        title = "2D_distributions"
        variables = ["redshift.value", "maxabsmag.value", "ebv.value",
                     "nbands", "hostoffsetdist.value"]
        variablenames = ["Redshift", "Max Absolute Magnitude", "EBV",
                         "Number of observational bands",
                         "Host Offset Distance"]

        pairIDs=[[0, 1], [0, 2], [0, 3], [0, 4]]
        self.scatter_data_2D(title, variables, variablenames, pairIDs)

    def plot_2d_fit_distributions(self, names=None, bands=None,
                                  title="2D_fit_distributions"):
        variables, variablenames, pairIDs = p.config_sd()
        self.scatter_data_2D(title, variables, variablenames, pairIDs, names,
                             bands, loose=True)

    def plot_2d_fit_with_peaks(self):
        names, bands = sc.return_candidates_with_peaks()
        title = "2D_fit_distribution_with_peaks"
        self.plot_2d_fit_distributions(names, bands, title)

    def plot_2d_fit_xray(self):
        bands, varnames = sc.xray_names()
        title = "2D_fit_distribution_xray"
        self.plot_2d_fit_distributions(names=None, bands=bands, title=title)

    def plot_fit_parameters(self, names=None, bands=None,
                            savepath="graphs/optical_fits/", root="combined/all"):
        """ Plots binned distributions of candidates for each variable in 'variables'.
        """
        full = [[],[]]
        titles, xlabels = lc.return_all_parameter_names()
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
                        if len(full[i]) <len(subset):
                            full[i].append([])
                        full[i][j].extend(subset[j])

        for i, name in enumerate(["", "_loose"]):
            print "\n"
            print name
            print "\n"
            print "Parameters \t Max \t \t Min \n"

            for j, param in enumerate(xlabels):
                print param, "\t \t","{0:.3}".format(float(min(full[i][j]))),\
                    "\t \t","{0:.3}".format(float(max(full[i][j])))
            print "\n"
            p.histograms_plot(titles, xlabels, full[i],
                              savepath + root + name + "_")

    def plot_fit_parameters_with_peaks(self):
        names, bands = sc.return_candidates_with_peaks()
        savepath="graphs/optical_fits/"
        root="combined/with_peak"
        self.plot_fit_parameters(names, bands, savepath, root)

    def plot_fit_parameters_xrays(self):
        xrbs, xvars = sc.xray_names()
        savepath = "graphs/xrays/"
        root="all"
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
            wrapAround180(ra), np.deg2rad(dec), c=rs, s=35, cmap=cm)
        cbar = plt.colorbar(sc)
        cbar.set_label('Redshift(z)')
        plt.show()

        path = "graphs/misc/skymap.pdf"
        print "Saving to", path
        plt.savefig(path)
        plt.close()

        print rs, max(rs), len(rs)

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
        self.name = jsonfile[0]["name"]
        self.has_photometry = False
        self.has_spectra = False
        self.ratiotest = False
        self.has_xray = False
        self.type = "TDE"
        self.best = False
        self.fits = dict()
        for key, val in jsonfile[0].items():
            if key != "photometry":
                if isinstance(val, list):
                    con = container()
                    for key2, val2 in val[0].items():
                        setattr(con, key2, val2)
                    setattr(self, key, con)
                else:
                    setattr(self, key, val)
            else:
                self.group_photometry(val)
        self.swiftdata()
        self.spectrum()
        self.coordinates()
        self.mjd_time()
        self.isbest()

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
            c = SkyCoord(self.ra.value, self.dec.value, unit=(u.hourangle, u.deg))
            self.ra_deg = c.ra.degree
            self.dec_deg = c.dec.degree

    def isbest(self):
        """Assesses whether a given candiate is promising.
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
            if len(getattr(self.photometry, band)) > minpoints:
                self.has_photometry = True

    def fit_xray(self):
        """Checks to see if there's a light curve visible in X Ray measurements.
        Fits this distribution to see if it is close to the characteristic
        t^-5/3.
        """
        xrbands, varnames = sc.xray_names()
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
                            if npoints > minpoints:
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
                            if npoints > minpoints:
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
                xrbs, xrvars = sc.xray_names()
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
