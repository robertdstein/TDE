# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:53:56 2017

@author: steinrob
"""
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

labels = ["A0", "a", "tt", "c", "offset"]
titles = ["Amplitude of peak", "Width of Peak", "Transition Time",
          "Power Law Index",
          "Time displacement of peak from maximum luminosity"]


def return_parameter_names():
    """Returns a list containing the full descriptive name of each light curve
    parameter, as well as its associated label in parametrisation

    :return: Title list, label list
    """
    return list(titles), list(labels)


def return_all_parameter_names():
    """Returns a list containing the full descriptive name of each light curve
    parameter, as well as its associated label in parametrisation. Also
    includes the Chi squared per degree of freedom as a parameter

    :return: Title list, label list
    """
    a = list(titles)
    a.append(r"$\chi^{2}$ per degree of freedom")
    b = list(labels)
    b.append("chi2_per_dof")
    return a, b


def return_histogram_labels():
    """Returns a list of the associated labels for each parameter"""
    return list(labels)


def return_parameter_bounds(maximum_luminosity=20):
    """Return the stricter parameter bounds for the fit, based on physical
    motivations and the general trend observed in the three peaked TDEs.

    :param maximum_luminosity: Maximum OBSERVED luminosity
    :return: Array containing bounds for parameters
    """
    return [(maximum_luminosity, maximum_luminosity + 3),
            (3 * 10 ** -4, 8 * 10 ** -3), (2., 350), (-8., -0.2),
            (-400, 400)]


def return_loose_bounds(maxlum=None):
    """Return the looser parameter bounds for the fit, based on purely
    on physical motivations.

    :return: Array containing bounds for parameters
    """
    return[(None,None), (10**-6, None), (2., 350),
           (None, -10**-6), (None, None)]


def default(maximum_luminosity=30):
    """Returns an array with typical values for each fit parameter,
    lying well within the strict parameter bounds

    :param maximum_luminosity: Maximum OBSERVED luminosity
    :return: Array containing starting values
    """
    return [maximum_luminosity, 5 * 10 ** -3, 20, -2., 0.0]


def logpeak(x, p=default()):
    """Returns a value for the peak part of the lightcurve, which is a
    Gaussian that in logspace becomes a quadratic function

    :param x: Shifted Time (x - offset)
    :param p: Fit parameter array
    :return: Value of lightcurve at time x for peak law
    """
    model = p[0] - p[1]*(x**2)
    return model


def logcontinuity(p=default()):
    """Function to find the point of intersection between the two models,
    and calculate the consequent gradient

    :param p: Fit Parameter Array
    :return: Values of Time, Light Curve and Gradient at transition
    """
    ytr = logpeak(p[2], p)
    xtr = p[2]
    gradtr = -2 * p[1] * p[2]
    return xtr, ytr, gradtr


def logpowerlaw(x, p=default()):
    """Returns a value for the Power Law part of the lightcurve,
    which requires calculation of the transition ponit, in order to ensure
    continuity of the two components

    :param x: Shifted Time (x - offset)
    :param p: Fit parameter array
    :return: Value of Lightcurve at time x for power law
    """
    xtr, ytr, gradtr = logcontinuity(p)
    power = p[3]
    x0 = xtr - power/gradtr
    b = ytr - power*np.log(xtr-x0)
    return b + power*np.log(x-x0)


def fitfunc(x_unshifted, p=default()):
    """Returns a value for the model at a time X_unshifted, using the
    parameter array p. Checks whether the point lies in the peak or power
    component, and returns the corresponding value

    :param x_unshifted: Unshifted time
    :param p: Fit parameter array
    :return: Value of Light Curve at x
    """
    x = x_unshifted+p[4]
    xtr, ytr, gradtr = logcontinuity(p)
    if x < xtr:
        return logpeak(x, p)
    else:
        return logpowerlaw(x, p)


def plot():
    """Plots the light curve model with various parameters, to show the range
    of possible forms.
    """
    xvals = np.arange(-50, 250, step=0.1)

    fig = plt.figure()
    plt.suptitle("Gaussian with smooth transition to power law")

    A0vals = [10, 11]
    avals = [5*10**-3, 10**-3, 5*10**-4]
    ttvals = [10., 50., 100.]
    cvals = [-0.1, -0.9, -5./3., -4.]
    offset = [-30, 0.0, 30]

    paramvals = [A0vals, avals, ttvals,cvals, offset]
    titles, labels = return_parameter_names()

    nplots = len(paramvals)

    for i in range(nplots):
        plt.subplot(nplots, 1, i+1)
        vals = paramvals[i]
        for j in range(len(vals)):
            pset = list(default())
            pset[i] = vals[j]
            yvals=[]
            ypower=[]
            ypeak=[]
            for x in xvals:
                yvals.append(fitfunc(x, pset))
                ypeak.append(logpeak(x,pset))
                if x > 0:
                    ypower.append(logpowerlaw(x,pset))
            label = labels[i] + "="+str(vals[j])
            plt.plot(xvals, yvals, label = label)

        plt.title(titles[i])
        plt.legend()

    fig.set_size_inches(15, 30)
    plt.savefig("graphs/misc/lightcurve_models.pdf")
    plt.close()


#plot()
