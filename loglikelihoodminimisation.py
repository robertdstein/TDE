# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 09:37:42 2017

@author: steinrob
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.markers import CARETDOWN
from scipy import optimize
import lightcurves as lc
from tqdm import tqdm

minpoints = len(lc.return_histogram_labels())

def run(fitdata, name):
	fitt = np.array(fitdata["t"])
	fity = np.array(fitdata["y"])
	
	print fity
	
	uncertaintyfrac = 0.1
	
	fitup = np.array(fitdata["up"])
	fitdown=np.array(fitdata["down"])
	fit_ul_t=np.array(fitdata["ul_t"])
	fit_ul_y=np.array(fitdata["ul_y"])

	if fitdata["var"] == "magnitude":
		fitup += uncertaintyfrac
		fitdown += uncertaintyfrac
		fity = -fity
		fit_ul_y = -fit_ul_y
		plt.ylabel(fitdata["var"])
		maxlum = max(fity)
		maxtime = fitdata["maxtime"]
	else:
		fitup = np.log10((1.+uncertaintyfrac)*fity + fitup) - np.log10(fity)
		fitdown = np.log10(fity) - np.log10(np.maximum(10**-6, ((1.-uncertaintyfrac)*fity - fitdown)))
		fity = np.log10(fity)
		fit_ul_y = np.log10(fit_ul_y)
		plt.ylabel("Log(" +  fitdata["var"]+")")
		maxlum = max(fity)
		mask = fity == maxlum
		maxtime = fitt[mask]
		
	fitsig = np.array(zip(fitup, fitdown))
	
	offset = np.linspace(-30.0, 120., num=6)
	A0vals = [maxlum]
	avals = [10**-3]
	ttvals = [10., 100.]
	cvals = [-0.2, -1.0, -2.0, -4.0]
	
	N = len(avals)*len(ttvals)*len(cvals)*len(offset)*len(A0vals)
	
	def llh(p, to_print=False):
		ll=0.0

		for k, t in enumerate(fitt):
			time = t - maxtime + p[4]
			y = fity[k]
			model = lc.fitfunc(time, p)
			if y > model:
				err=fitsig[k][0]
			else:
				err=fitsig[k][1]
			ll += ((y - model)/err)**2
			if to_print:
				print time, y, model, err, ((y - model)/err)**2, ll
					
		for k, ul_t in enumerate((fit_ul_t)):
			time = ul_t - maxtime + p[4]
			y = fit_ul_y[k]
			model = lc.fitfunc(time, p)
			if y < model:
				ll += -np.log(0.05)				
		return ll
	
	info = dict()
	suffix=["", "_loose"]
	
	for j, f in enumerate([lc.return_parameter_bounds, lc.return_loose_bounds]):
		
		ax = plt.subplot(2,1,j+1)

		ax.errorbar(fitt, fity, yerr=[fitdown, fitup], label="measurements", color="b", ecolor="r", fmt='o')
		ax.axvline(maxtime, label="Peak Visual Magnitude",color="green", linestyle="--")
		
		plt.xlabel("t")
		
		params=[]
		covars = []
		vals = []
	
		with tqdm(total=N) as pbar:
			for i in offset:
				pbest = lc.default(maxlum)
				
				fittedpeaktime = maxtime-i
	
				bnds = lc.return_parameter_bounds(maxlum)
				bestout = optimize.minimize(llh, pbest, bounds=bnds)		
				
				llbest = np.sum(bestout.fun)
				
				for A0 in A0vals:
					for a in avals:
						for tt in ttvals:
							for c in cvals:
								pinit = [A0, a, tt, c, i]
								bnds = f(A0)
								out = optimize.minimize(llh, pinit, method='L-BFGS-B', bounds=bnds)
								ll = np.sum(out.fun)
								pbar.update(1)
								if ll < llbest:
									llbest= ll
									bestout=out
									pbest=pinit
				
				pfinal = bestout.x
				covar = bestout.success					
				covars.append(covar)
				params.append(pfinal)
				vals.append(bestout.fun)
				
				best = min(vals)
				index = vals.index(best)
				bestparams = params[index]
		
		best = min(vals)
		index = vals.index(best)
	
		bestparams = params[index]
		
		for i, var in enumerate(lc.return_histogram_labels()):
			info[var+suffix[j]] = bestparams[i]
			
		info["llh"+suffix[j]] = best
		
		npoints = len(fitt) + len(fit_ul_t)
		if len(fitt) > minpoints:
			info["chi2_per_dof"+suffix[j]] = float(best)/float(npoints - minpoints)
		else:
			raise Exception("Too few datapoints (" + str(len(fitt)) + ") for fitting a model with (" + str(minpoints) + ") parameters")
		
		covar = covars[index]
		
		t_max_to_peak=bestparams[4]
		
		fittedpeaktime = maxtime-t_max_to_peak
		
		plottimes = np.linspace(fitt[0]-100, fitt[-1]+100, 5000)		
		shiftedtimes = plottimes-fittedpeaktime
				
		def fit(times):
			models=[]
			for t in times:
				models.append(lc.fitfunc(t, bestparams))
			return models
			
		title= name +suffix[j]+ " (["
		for a in bestparams:
			title += "{0:.6}".format(str(a))
			title += ", "
			
		title += "] " + "{0:.4}".format(str(float(info["chi2_per_dof"+suffix[j]]))) +")" 
		
		plt.title(title)
	
		ax.plot(plottimes, fit(shiftedtimes), label="Fit")
		
		ax.axvline(fittedpeaktime -t_max_to_peak + bestparams[2], label="Model Transition",color="purple", linestyle="--")
		
		if fitdata["var"] == "magnitude":
			plt.ylabel("Magnitude")
		
		else:
			plt.ylabel("Log("+ fitdata["var"]+ ")")
		
		if len(fit_ul_t) > 0:
			ax.scatter(fit_ul_t, fit_ul_y, marker=CARETDOWN, s=150, label="Upper Limits")	

		lowerlim = min(fity)-1
		upperlim = max(fity) + 1
	#		
		ax.set_ylim(bottom=lowerlim, top=upperlim)
		ax.set_xlim(fitt[0]-(100+np.abs(t_max_to_peak)))
		ax.legend()		
		
		plt.legend()
		
		
	
	return info
