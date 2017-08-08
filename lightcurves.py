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
titles =["Amplitude of peak", "Width of Peak", "Transition Time", "Power Law Index", "Time displacement of peak from maximum luminosity"]

def return_parameter_names():
	return list(titles), list(labels)
	
def return_all_parameter_names():
	a = list(titles)
	a.append(r"$\chi^{2}$ per degree of freedom")
	b = list(labels)
	b.append("chi2_per_dof")
	return a, b
	
def return_histogram_labels():
	return list(labels)
	
def return_parameter_bounds(maxlum=20):
	return [(maxlum, maxlum+3), (3*10**-4, 8*10**-3), (2., 350), (-8., -0.2), (-400, 400)]

def return_loose_parameter_bounds():
	return[(None, None), (10**-6, None), (1., None), (-15., 0.01), (None, None)]
	
def default(maxlum= 30):
	return [maxlum, 5*10**-3, 20, -2., 0.0]

def logpeak(x, p=default()):
	model = p[0] - p[1]*(x**2)
	return model
	
def logcontinuity(p=default()):
	ytr = logpeak(p[2], p)
	xtr = p[2]
	gradtr = -2 * p[1] * p[2]
	return xtr, ytr, gradtr
		
def logpowerlaw(x, p=default()):
	xtr, ytr, gradtr = logcontinuity(p)
	power = p[3]
	x0 = xtr - power/gradtr
	b = ytr - power*np.log(xtr-x0)
	return b + power*np.log(x-x0)
		
def fitfunc(x_unshifted, p=default()):
	x = x_unshifted+p[4]
	xtr, ytr, gradtr = logcontinuity(p)
	if x < xtr:
		return logpeak(x, p)
	else:
		return logpowerlaw(x, p)
		
def plot():
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
	
