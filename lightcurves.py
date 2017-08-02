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

#[Maxlum, a, transitiontime, c, f]
pinit = [10**43, 10**-3, 100, -5./3., 10**-8, 0.0]

def peak(x, p=pinit):
	model = p[0]*(np.exp(-p[1]*(x)**2))
	return np.maximum(model, 1.0)
	
def continuity(p=pinit):
	ytr = peak(p[2], p)
	xtr = math.sqrt(-np.log(ytr/p[0])/p[1])
	gradtr = -2 * xtr * ytr*p[1]
	return xtr, ytr, gradtr
		
def powerlaw(x, p=pinit):
	xtr, ytr, gradtr = continuity(p)
	power = p[3]
	if power != 0.0: 
		b = ((gradtr/power)**power)/(ytr ** (power-1))
	else:
		b=ytr	
	if power != 0.0:
		x0 = xtr - (ytr/b)**(1./(power))
	else:
		x0 = 0.0
	return b*((x-x0)**power)
	
				
def fitfunc(x_unshifted, p=pinit):
	x = x_unshifted+p[5]
	xtr, ytr, gradtr = continuity(p)
	if x < xtr:
		return peak(x, p) + p[4]*p[0]
	else:
		return powerlaw(x, p) + p[4]*p[0]
		
def plot():
	xvals = np.arange(-100, 500, step=0.1)
	
	fig = plt.figure()
	plt.suptitle("Gaussian with smooth transition to power law")
	
	A0vals = [10**43, 10**44]
	avals = [10**-2, 10**-3, 10**-4]
	ttvals = [10., 50., 100.]
	cvals = [-0.1, -0.9, -5./3., -4.]
	fvals = [0.0, 10**-5, 10**-3]
	offset = [-100, 0.0, 100, 200]
	
	labels = ["A0", "a", "tt", "c", "f", "offset"]
	paramvals = [A0vals, avals, ttvals,cvals, fvals, offset]
	titles = ["Amplitude of peak", "Width of Peak", "Transition Time", "Power Law Index", "Underlying Luminosity as fraction of maximum", "Time displacement of peak from maximum luminosity"]
	
	nplots = len(paramvals)
	
	for i in range(nplots):
		plt.subplot(nplots, 1, i+1)
		vals = paramvals[i]
		for j in range(len(vals)):
			pset = list(pinit)
			pset[i] = vals[j]
			yvals=[]
			ypower=[]
			ypeak=[]
			for x in xvals:
				yvals.append(fitfunc(x, pset))
				ypeak.append(peak(x,pset))
				if x > 0:
					ypower.append(powerlaw(x,pset))
			label = labels[i] + "="+str(vals[j])
			plt.plot(xvals, yvals, label = label)
		
		plt.yscale("log")
		plt.title(titles[i])
		plt.gca().set_ylim(bottom=(pinit[0]*10**-10))
		plt.legend()
	
	fig.set_size_inches(15, 30)						
	plt.savefig("graphs/lightcurve_models.pdf")
	plt.close()
	
plot()