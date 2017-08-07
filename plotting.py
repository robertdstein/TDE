# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 13:55:08 2017

@author: steinrob
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.markers import CARETDOWN
import lightcurves as lc

def scatter_plot(folder, title, allresults, combined=False, tshift=False):
	"""Produces scatter plots for a given set of data.
	""" 	
	fig = plt.figure()
	plt.suptitle(title)
	npoints = len(allresults)

	if npoints > 1:			
		ncolumns=2
	else:
		ncolumns=1
		
	if not combined:
		nrows = int(float(npoints)/2.) + (npoints % 2)
		fig.set_size_inches(ncolumns*7, nrows*5)
		fig.subplots_adjust(hspace=.5)
	else:
		nrows = 1
		ax=plt.subplot(1,1,1)
		plt.gca().invert_yaxis()
		fig.set_size_inches(10, 7)
		fig.subplots_adjust(hspace=.5)
	
	if npoints > 0:
		for i in range(0, npoints):
			res=allresults[i]
			
			if tshift:
				maxdate = res["maxdate"][1]
				sourcet = np.array(res["time"])
				check = sourcet>(maxdate-150)
				t = sourcet[check] - maxdate
				y = np.array(res["y"])[check]
				sigup = np.array(res["sigma_upper"])[check]
				sigdown=np.array(res["sigma_lower"])[check]
				plt.xlabel("Time since maximum visual luminosity (days)")
			else:
				t = res["time"]
				y = res["y"]
				upt=res["upper_limit_time"]
				upy = res["upper_limit_y"]
				sigup = res["sigma_upper"]
				sigdown=res["sigma_lower"]
				plt.xlabel("Time (MJD)")

			label = res["label"]
			var = res["variable"]
			
			md = res["maxdate"]
			
			if not combined:
				ax = plt.subplot(nrows, ncolumns,i+1)
				plt.title(label)
				if md[0] != None:
					plt.axvline(md[0], label="max date", color="orange")
					x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=md[0])
					ax.xaxis.set_major_formatter(x_formatter)
			
				if md[1] != None:
					plt.axvline(md[1], color = "purple", label="max vis date")
				if len(upt) > 0:
					ax.scatter(upt, upy, marker=CARETDOWN, s=150, label="upper limits")
				if len(t) > 0:
					ax.errorbar(t, y, yerr=[sigup, sigdown], label="measurements", color="green", ecolor="red", fmt="o")

				if var == "magnitude":
					plt.gca().invert_yaxis()
				else:
					ax.set_yscale('log')
			
			elif len(t) > 0:
				ax.errorbar(t, y, yerr=[sigup, sigdown], label=label, fmt="o")							
			
			plt.ylabel(var)
			ax.legend()
		
		path = "graphs/" + folder + title + ".pdf"
		
		plt.savefig(path)
		plt.close()
		print "Saving to", path
		
	else:
		print "Not Saved!"

def histograms_plot(titles, xlabels, values, path=None, bnds=False):
	"""Plots histograms for a given set of variables
	"""
	fig = plt.figure()
	npoints = len(titles)
	nrows = int(float(npoints)/2.) + (npoints % 2)
	
	if bnds:
		bnds = lc.return_parameter_bounds(max(values[0]))
	
	if npoints > 0:
		for i in range(0, npoints):
			ax = plt.subplot(nrows,2, i+1)
			
			nbins = np.maximum(int(float(len(values[i]))/2.) + 1, 20)
			
			if bnds:
				xmin = bnds[i][0]
				xmax = bnds[i][1]
				if xmin == None:
					xmin = min(values[i])*0.9
				if xmax == None:
					xmax  = max(values[i])*1.1
				histrange=[xmin, xmax]
			
			else:
				histrange=None
				
			n, bins, _ = plt.hist(values[i], bins=nbins, range=histrange, histtype='stepfilled', label=xlabels[i])	
	
			plt.title(titles[i])
			plt.xlabel(xlabels[i])
			
	fig.set_size_inches(25, nrows*5)
	fig.subplots_adjust(hspace=.5)
	
	if path == None:	
		path = "graphs/histogram.pdf"
		
	else:
		split = path.split("/")
		name  = split[-1].split(".")
		title = name[0]
		plt.suptitle(title)
	
	plt.savefig(path)
	print "Saving to", path
	plt.close()