# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 13:55:08 2017

@author: steinrob
"""
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.markers import CARETDOWN
import lightcurves as lc
import readsupernova as rs


def config_for_fit_scatter():
    info = dict()
    info["variables"] = ["fits.*.chi2_per_dof", "fits.*.tt",  "fits.*.A0",
                         "fits.*.a", "fits.*.offset", "fits.*.c"]
    info["variable_names"] = [
        r"$\chi^{2}$ per degree of freedom",
                             "Transition Times",
                             r"$A_{0}$ in $A_{0} \exp^{-a x^{2}}$",
                             r"a in $A_{0} \exp^{-a x^{2}}$",
                     "Displacement from peak to Brightest Observation Date",
                     r"$\gamma$"]
    info["pairIDs"]=[
        [1, 2, 0], [1, 3, 0], [4, 2, 0], [5, 2, 0], [5, 1, 0], [5, 3, 0]]
    return info


def config_for_general_scatter():
    info = dict()
    info["variables"] = ["redshift.value", "maxabsmag.value", "ebv.value",
                         "nbands", "hostoffsetdist.value"]

    info["variable_names"] = ["Redshift", "Max Absolute Magnitude", "EBV",
                              "Number of observational bands",
                              "Host Offset Distance"]

    info["pairIDs"] = [[0, 1], [0, 2], [0, 3], [0, 4]]
    return info


def scatter_distribution(folder, title, allresults):
    """Plots scatter distributions
    """
    fig = plt.figure()
    plt.suptitle(title)
    npoints = len(allresults)

    if npoints > 1:
        ncolumns=2
    else:
        ncolumns=1

    nrows = int(float(npoints)/2.) + (npoints % 2)
    fig.set_size_inches(ncolumns*7, nrows*5)
    fig.subplots_adjust(hspace=.5)

    if npoints > 0:
        for i in range(0, npoints):
            res = allresults[i]

            ax = plt.subplot(nrows, ncolumns,i+1)

            plt.xlabel(res["x_label"])
            plt.ylabel(res["y_label"])

            if "z" in res.keys():
                cm = plt.cm.get_cmap('RdYlGn_r')
                if res["zvar"].split(".")[-1] == "chi2_per_dof":
                    ul = 2.0
                    ll = 0.0
                else:
                    ul = None
                    ll = None

                scattergraph = plt.scatter(
                    res["x"], res["y"], c=res["z"], vmin=ll, vmax=ul, cmap=cm)
                cbar = plt.colorbar(scattergraph)
                cbar.set_label(res["z_label"])

            else:
                ax.scatter(res["x"], res["y"])
            ax.legend()

        path = "graphs/" + folder + title + ".pdf"

        plt.savefig(path)
        plt.close()
        print "Saving to", path


def histograms_plot(titles, xlabels, values, path=None, bnds=False, suffix=""):
    """Plots histograms for a given set of variables
    """
    fig = plt.figure()
    npoints = len(titles)
    nrows = int(float(npoints)/2.) + (npoints % 2)

    if npoints > 0:
        for i in range(0, npoints):
            ax = plt.subplot(2, nrows, i+1)

            nbins = np.minimum(
                np.maximum(int(float(len(values[i]))/4.) + 1, 5), 20)

            median = np.median(values[i])

            if not math.isnan(median):
                n, bins, _ = plt.hist(values[i], bins=nbins,
                                      histtype='stepfilled', label=xlabels[i])

                plt.title(
                    titles[i] + " (Median = " + "{0:.3}".format(median) + ")")
                plt.xlabel(xlabels[i])

    fig.set_size_inches(nrows*4, 7)
    fig.subplots_adjust(hspace=.5)

    if path is None:
        path = "graphs/misc/histogram" + suffix + ".pdf"
    else:
        split = path.split("/")
        name = split[-1].split(".")
        title = name[0]+suffix
        plt.suptitle(title)
        path += "histogram" + suffix + ".pdf"

    plt.tight_layout()
    plt.savefig(path)
    print "Saving to", path
    plt.close()


def scatter_photometry(folder, title, allresults, combined=False, tshift=False,
                       sn=False):
    """Produces scatter plots for a given set of data.
    """
    fig = plt.figure()
    plt.suptitle(title)
    n_points = len(allresults)

    if n_points > 1:
        n_columns = 2
    else:
        n_columns = 1

    if not combined:
        n_rows = int(float(n_points)/2.) + (n_points % 2)
        fig.set_size_inches(n_columns*5, n_rows*3)
        fig.subplots_adjust(hspace=.5)
    else:
        n_rows = 1
        ax = plt.subplot(1, 1, 1)
        plt.gca().invert_yaxis()
        fig.set_size_inches(10, 7)
        fig.subplots_adjust(hspace=.5)

    if n_points > 0:
        for i in range(0, n_points):
            res = allresults[i]

            if tshift:
                maxdate = res["maxdate"][1]
                sourcet = np.array(res["time"])
                check = sourcet > (maxdate-150)
                t = sourcet[check] - maxdate
                y = np.array(res["y"])[check]
                sigup = np.array(res["sigma_upper"])[check]
                sigdown=np.array(res["sigma_lower"])[check]
                plt.xlabel("Time since maximum visual luminosity (days)")
            else:
                t = res["time"]
                y = res["y"]
                upt = res["upper_limit_time"]
                upy = res["upper_limit_y"]
                sigup = res["sigma_upper"]
                sigdown = res["sigma_lower"]
                plt.xlabel("Time (MJD)")

            label = res["label"]
            var = res["variable"]

            md = res["maxdate"]

            if not combined:
                ax = plt.subplot(n_rows, n_columns,i+1)
                plt.title(label)
                if md[0] is not None:
                    plt.axvline(md[0], label="max date", color="orange")
                    x_formatter = matplotlib.ticker.ScalarFormatter(
                        useOffset=md[0])
                    ax.xaxis.set_major_formatter(x_formatter)

                if md[1] is not None:
                    plt.axvline(md[1], color="purple", label="max vis date")
                if len(upt) > 0:
                    ax.scatter(
                        upt, upy, marker=CARETDOWN, s=150, label="upper limits")
                if len(t) > 0:
                    ax.errorbar(
                        t, y, yerr=[sigup, sigdown], label="measurements",
                        color="green", ecolor="red", fmt="o")

                if var == "magnitude":
                    plt.gca().invert_yaxis()
                else:
                    ax.set_yscale('log')

            elif len(t) > 0:
                ax.errorbar(t, y, yerr=[sigup, sigdown], label=label, fmt="o")

            plt.ylabel(var)
            ax.legend()

        if sn and (len(title) == 1):
            sn_data = rs.run()
            lims = ax.get_ylim()
            for data in sn_data:
                if title in data.keys():
                    ydata = data[title]
                    ymax = max(ydata)
                    ymin = min(ydata)

                    diff = ymin - lims[1] + 3.0
                    plt.plot(
                        data["t"], ydata - diff, linestyle="--", marker="",
                        label=data["name"])

                    print lims, ymax, ymin, diff

            plt.gca().set_ylim(bottom=lims[0])
            ax.legend()

        path = "graphs/" + folder + title + ".pdf"

        plt.savefig(path)
        plt.close()
        print "Saving to", path

    else:
        print "Not Saved!"
