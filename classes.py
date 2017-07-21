import zipfile
import math
import json
import datetime
import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import optimize
from scipy.optimize import curve_fit
from numpy import vectorize
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

minpoints = 3

def wrapAround180(ra_deg):
    ra = np.deg2rad(ra_deg)
    ra[ra > np.pi] -= 2*np.pi
    return ra
				
def magtolum(absmag):
	"""Convert an absolute magnitude to a distance.
	"""
	solarlum = 3.826 * 10**33
	solarabsmag = 4.75
	exponent = (solarabsmag - float(absmag))/2.5
	lum = solarlum * 10 **(exponent)
	return lum

class Full_set:
	"""The full set of all TDE candidates
	"""
	def __init__(self, path):
		self.path=path
		self.TDEs= container()
		self.totalentries=0
		self.candidates_with_photometry=0
		self.prime_candidates=0
		self.all_available_bands=[]
		self.best_list=[]
		self.extract()
		
	def extract(self):
		"""Extract all data from tde.zip.
		The zipfile contains individual json files for each candidate.
		"""
		zf = zipfile.ZipFile(self.path, 'r')
		self.creation_time=datetime.datetime(*zf.infolist()[0].date_time)
		
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
		print "Dataset fully extracted. Of", self.totalentries, "we have", self.candidates_with_photometry, "with at least four photometric observations in one channel."
		print "The following channels were covered \n", str(self.all_available_bands)
		self.list_best()
		
	def creationtime(self):
		"""Returns the time that the online catalogue was last updated
		"""
		print "Catalogue version created:", self.creation_time, "\n"
		
	def info_candidate(self, name):
		"""Plots all available light curves for a given candidate.
		"""
		if hasattr(self.TDEs, name):
			tde = getattr(self.TDEs, name) 
			print "TDE found, with photometry =", tde.has_photometry
			tde.plot_candidate()
		else:
			raise Exception("No TDE with this name was found!")
				
	def plot_all_candidates(self):
		"""Plots all available light curves for each candidate.
		"""
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			tde.plot_candidate()
			
	def plot_all_bands(self):
		"""Plots all available light curves for each band.
		"""
		for band in self.all_available_bands:
			title = band
			
			folder = "bands/"
			
			allt = []
			ally = []
			allupt = []
			allupy = []
			labels =[]
			variables = []
			maxdate = []
			
			for name in vars(self.TDEs):
				tde = getattr(self.TDEs, name)
				if tde.has_photometry:
					if band in tde.bandlist:
						i = tde.bandlist.index(band)
						
						var = tde.varlist[i]
						t, y, upt, upy, md = tde.return_photometry_dataset(band, var)
						
						labels.append(tde.name)
						allt.append(t)
						ally.append(y)
						allupt.append(upt)
						allupy.append(upy)
						variables.append(var)
						maxdate.append(md)
				
			scatter_plot(folder, title, allt, ally, allupt, allupy, labels, variables, maxdate, sharex=True)
		
	def list_best(self):
		"""Lists the TDE candidates which are most promising.
		"""
		best=[]
		ratio=[]
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
		
	
	def plot_distibutions(self):
		""" Plots binned distributions of candidates for each variable in 'variables'.
		"""
		variables = ["redshift", "comovingdist", "lumdist", "hostoffsetdist", "mjdmax", "nbands","maxabsmag",  "nhhost", "ra_float", "dec_float", "ebv"]
		
		titles = ["Redshift", "Comoving Distance", "Luminosity Distance", "Host Offset Distance", "Date of Max Luminosity", "Number of observational bands", "Max Absolute Magnitude", "Something cool, probably (nnhost)", "Right ascension", "Declination", "EBV"]
		xlabels = ["Redshift (z)", "Distance (MPc)", "Distance (MPc)", "Distance (?)", "Time(MJD)", "n", "", "???", "hours", "degrees","Magnitudes"]
		
		values = []
		for i in range (0, len(variables)):
			values.append([])
		
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			if tde.has_photometry:
				for i in range (0, len(variables)):
					var = variables[i]
					if hasattr(tde, var):
						val = getattr(tde, var)
						if not isinstance(val, float):
							val = getattr(tde, var).value
						values[i].append(float(val))
		
		print "Plotting histograms for the following:", titles
		histograms_plot(titles, xlabels, values)
		
	def plot_spectra(self):
		"""Plots the spectrum for each candidate  which has one
		"""
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			if tde.has_spectra:
				tde.plot_spectrum()
				
	def plot_skymap(self):
		"""Plots a map of the distribution of TDE candidates on the sky.
		Uses the redshift as a colour scale.
		"""
		ra=[]
		dec = []
		rs=[]
		markers=[]
		
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			if hasattr(tde, "ra_deg") and hasattr(tde, "redshift"):
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
		sc = plt.scatter(wrapAround180(ra),np.deg2rad(dec),c=rs, s=35, cmap=cm)
		cbar = plt.colorbar(sc)
		cbar.set_label('Redshift(z)')
		plt.show()
		
		path = "graphs/skymap.pdf"
		print "Saving to", path
		plt.savefig(path)
		plt.close()
		
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
		self.best=False
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
		self.spectrum()	
		self.coordinates()
		self.mjd_time()
		self.isbest()
		
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
		self.check_xray()
		if self.has_photometry:
			if hasattr(self, "mjdmax") and hasattr(self, "mjdvismax"):
				if self.nbands > 2.0:
					self.best=True
	
	def mjd_time(self):
		"""Converts the discovery date to MJD units.
		Does the same for the maximum date, and maximum visual date.
		"""
		if hasattr(self, "discoverdate"):
			val = self.discoverdate.value.split("/")
			if len(val)==3:
				[y, m, d] = val
				d = d.split(".")[0]
				discdate = y+"-"+m+ "-"+ d+"T00:00:00"
				t = Time(discdate)
				setattr(self, "mjddisc", float(t.mjd))
		
		if hasattr(self, "maxdate"):
			[y, m, d] = self.maxdate.value.split("/") 
			maxdate = y+"-"+m+ "-"+ d+"T00:00:00"
			t = Time(maxdate)
			setattr(self, "mjdmax", float(t.mjd))

		if hasattr(self, "maxvisualdate"):
			[y, m, d] = self.maxvisualdate.value.split("/") 
			maxdate = y+"-"+m+ "-"+ d+"T00:00:00"
			t = Time(maxdate)
			setattr(self, "mjdvismax", float(t.mjd))
					
	def group_photometry(self, val):
		"""Extracts the photometry data from the json file, and groups it
		into seperate claases for each available band.
		"""
		allbands = container()
		bandlist = []
		varlist=[]
		for entry in val:
			if "band" in entry.keys():
				band = entry["band"]
				variable = "magnitude"
				
			elif "energy" in entry.keys():
				band = "X-Ray ("+str(entry["energy"][0]) + "-" + str(entry["energy"][1]) + ") " + str(entry["u_energy"])
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
			
			if hasattr(allbands, band+"_"+variable):
				pass
			else:
				setattr(allbands, band+"_"+variable, [])
				bandlist.append(str(band))
				varlist.append(str(variable))
				
			getattr(allbands, band+"_"+variable).append(entry)
		
		self.bandlist = bandlist
		self.varlist = varlist
		self.photometry = allbands
		self.nbands = float(len(bandlist))
		for band in vars(allbands):
			if len(getattr(self.photometry, band)) > minpoints:
				self.has_photometry = True
				
	def check_xray(self):
		"""Checks to see if there is a light curve visible in X Ray measurements.
		Fits this distribution to see if it is close to the charecteristic t^-5/3.
		"""
		xrband = "X-Ray (0.3-2.0) keV"
		if hasattr(self, "bandlist"):
			if xrband in self.bandlist:
				data = getattr(self.photometry, xrband+"_"+"luminosity")					
				
				time=[]
				y=[]
				sigma=[]
				upt = []
				upy =[]
				upsig=[]
				lowersig = []
				
				for entry in data:
					if "time" in entry.keys():
						t = entry["time"]
						if isinstance(t, list):
							t = t[0]
						if "upperlimit" in entry.keys():
							upy.append(float(entry["luminosity"]))
							upt.append(float(t))
						else:
							lum = float(entry["luminosity"])
							y.append(lum)
							time.append(float(t))
							meansigma = 0.5*(float(entry["e_upper_luminosity"]) + float(entry["e_lower_luminosity"]))
							if meansigma==0.0:
								meansigma = 0.5*lum
								upsig.append(0.5*lum)
								lowersig.append(0.5*lum)
								
							sigma.append(meansigma)
							upsig.append(float(entry["e_upper_luminosity"]))
							lowersig.append(float(entry["e_lower_luminosity"]))
					
				if len(y) > 0:
					maxlum = float(max(y))
					maxtime = time[y.index(maxlum)]
					
					print self.name, maxlum, y, maxtime, time
					
					allt = time+upt
					
					lastt=0.0			
					
					for j in range(len(allt)):
						t = allt[j]
						if lastt < t < maxtime:
							lastt=t
							
					preceding = []
					for i in range(len(upt)):
						t = upt[i]
						if float(t) < (maxtime-50):
							preceding.append(float(upy[i]))
							
#					print preceding, "Max lum:", maxlum
					
					if len(preceding) > 0:
						rmax = float(maxlum)/float(max(preceding))					
						rmean = float(maxlum)/float(np.mean(preceding))
						rmin= float(maxlum)/float(min(preceding))
						
						if  rmax> 1.5:
							if rmean > 4:
								if rmin > 20:
									self.ratiotest=True
								
#						print self.name, rmax, rmean, rmin, self.ratiotest
#					else:
#						print self.name, "has insufficient data"
					
#					if self.ratiotest:
#						print self.name, "has passed the ratio test!"
#					else:
#						print self.name, "has FAILED the ratio test!"
						
					fitt=[]
					fity=[]
					fitsig=[]
					fitup=[]
					fitdown=[]
					for i in range(len(time)):
						t = time[i]
						if float(t) >= maxtime:
							fitt.append(float(t))
							fity.append(float(y[i]))
							fitsig.append(sigma[i])
							fitup.append(upsig[i])
							fitdown.append(lowersig[i])
							
							
					if len(fitt) > 2:
						
						fig = plt.figure()	
						ax1 = plt.subplot(312)
						
						if hasattr(self, "mjdmax"):
							plt.axvline(self.mjdmax, color="k", label="max date")
							
						if hasattr(self, "mjdmax"):
							plt.axvline(self.mjdvismax, color = "green", label="max vis date")
							
						plt.suptitle(self.name)
						ax1.set_yscale('log')
						
						fitt = np.array(fitt)
						fity = np.array(fity)
						fitsig = np.array(zip(fitup, fitdown))
						
#						print fity, fitt, fitup, fitdown
#						raw_input("prompt")
						
						ax1.errorbar(fitt, fity, yerr=[fitup, fitdown], label="measurements", color="b", ecolor="r", fmt='o')
						plt.xlabel("t")
						plt.ylabel("Log(Luminosity)")
						plt.legend()
						
						params=[]
						covars = []
						vals = []
						
						scale = np.log10(maxtime  - lastt)			
						
						offset = np.logspace(0, scale, num=200)
#						offset = np.array(range())
#						print "Range", offset
						
						for i in offset:
							starttime = maxtime - (i + 1e-5)
#							if hasattr(self, "mjdmax"):
#								starttime = self.mjdmax
							
							newt = fitt-starttime
							logfitt = np.log10(newt)
							
							def indices(y,model):
								if y > model:
									return 0
								else:
									return 1
									
							vfunc = vectorize(indices)
							
							
							def fitfunc(p, x):
								return (p[0]*(x**p[1]))
								
							def err(indices, errs):
								return errs[indices]
								
							verr = vectorize(err)
								
							def llh(p, x, y, err2):
								model = fitfunc(p, x)
								err_indices = vfunc(y, model)
								err_range = np.arange(len(err_indices))
								print err_indices, err_range
								err = verr(err_indices, err2)
								print err2, err
								raw_input("prompt")
								llh = ((y - model)/err)**2
								return llh
							
							pinit = [maxlum, -1.0]
							out = optimize.leastsq(llh, pinit,
							                       args=(newt, fity, fitsig), full_output=1)
							
							pfinal = out[0]
							covar = out[1]
							covars.append(covar)
							params.append(pfinal)

							vals.append(np.sum(llh(pfinal, newt, fity,fitsig)))
						
#						print vals
						
						ax2 = plt.subplot(311)
						
						
						
						best = min(vals)
						index = vals.index(best)
						time = maxtime - (offset[index] + 1e-5)

						lower=1.
						upper=(maxtime-lastt)
						mid=False
						
						for j in range(len(vals)):
							v = vals[j]
							if not mid:
								if v > (0.5+best):
									lower = offset[j]
								else:
									mid = True
							else:
								if v > (0.5+best):
									upper = offset[j]
									break
								
						print "Done!", lower, upper, offset[index], (lower+upper)/2, (upper-lower)
						
						plt.plot(offset, np.array(vals)-best)
						plt.xlabel("Log(Age of TDE at observed maximum brightness)")
						plt.ylabel(r"$\Delta$ Log Likelihood (Gaussian)")
						ax2.set_xscale("log")
						
						plt.axvspan(lower, upper, color="r", alpha=0.3)
						plt.scatter(offset[index], 0.0, color="r", marker="*")
						plt.axhline(0.5, color="g")

						bestparams = params[index]
						covar = covars[index]
						
#						print best, pfinal, offset[index]
#						raw_input("prompt")

						plottimes = np.linspace(fitt[0], fitt[-1], 200)
						shiftedtimes = plottimes-time
						
						def fit(x):
							return bestparams[0] * (x**bestparams[1])
						
						plt.title("Minimum occurs at " + str(int(time)) + " MJD (" + str(int(offset[index]))+" days before date of observed maximum, with 1 " r"$\sigma$ range of " + str(int(lower)) + " to "+ str(int(upper)) + " days before maximum)")							
						if type(covar) != type(None):							
							ax1.annotate("Best fit has power law with gradient " + "{0:.3g}".format(bestparams[1]) + " +/- " + str(math.sqrt(np.abs(covar[1][1]))), xy=(0.05, 0.1), xycoords="axes fraction")
						else:
							ax1.annotate("Best fit has power law with gradient " + "{0:.3g}".format(bestparams[1]) + " but failed to find covariance (fit did not converge correctly)", xy=(0.05, 0.1), xycoords="axes fraction")
							
						ax3 = plt.subplot(313)
						ax3.errorbar(fitt-time, fity, yerr=[fitup, fitdown], label="measurements", color="b", ecolor="r", fmt='o')
						ax3.plot(shiftedtimes, fit(shiftedtimes)) 
						ax1.plot(plottimes, fit(shiftedtimes))
						plt.ylabel("Log(Luminosity)")
						plt.xlabel("Log(t-"+ str(int(time))+" MJD)")
						
						ax3.set_yscale('log')
						ax3.set_xscale('log')
							
						fig.set_size_inches(15, 20)						
						
						path = "graphs/xrays/" + self.name+ ".pdf"						
						
						plt.savefig(path)
						plt.close()
						print "Saving to", path

				
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
			f=[]
			
			spec = self.spectra
			
			for entry in spec.data:
				f.append(entry[0])
				wl.append(entry[1])
				
			plt.plot(f, wl)
			
			plt.xlabel(spec.u_wavelengths)
			plt.ylabel(spec.u_fluxes)
			
			title = self.name 
			if hasattr(spec, "redshift"):
				title += " [z =" + str(spec.redshift)+ "]"
			if hasattr(self, "mjdmax") & hasattr(spec, "time"):
				time_after_max = int(float(spec.time) - float(self.mjdmax))
				if time_after_max < 0:
					title += " [" + str(-1 * time_after_max) +" days before max]"
				else:
					title += " [" + str(time_after_max) + " days after max]"
			
			plt.title(title)
			
			path = "graphs/spectra/"+self.name+".pdf"
			print "Saving to", path
			plt.savefig(path)
			plt.close()
				
	def return_photometry_dataset(self, band, var):
		if hasattr(self.photometry, band+"_"+var):
			data = getattr(self.photometry,band+"_"+var)
			
			time=[]
			y=[]
			upt = []
			upy =[]
			
			if hasattr(self, "mjdmax"):
				md = [self.mjdmax]
			else:
				md = [None]
				
			if hasattr(self, "mjdvismax"):
				md.append(self.mjdvismax)
			else:
				md.append(None)
			
			for entry in data:
				if "time" in entry.keys():
					t = entry["time"]
					if isinstance(t, list):
						t = t[0]
					if "upperlimit" in entry.keys():
						upy.append(float(entry[var]))
						upt.append(float(t))
					else:
						y.append(float(entry[var]))
						time.append(float(t))
			return time, y, upt, upy, md
		else:
			return None, None, None, None, None
			
	def plot_candidate(self):
		if self.has_photometry:
			folder = "candidates/"
			title = self.name
			
			allt = []
			ally = []
			allupt = []
			allupy = []
			labels =[]
			variables = []
			maxdate=[]
			
			for band in sorted(self.bandlist, key=str.lower):
				i = self.bandlist.index(band)
				
				var = self.varlist[i]
				t, y, upt, upy, md = self.return_photometry_dataset(band, var)
				
				labels.append(band+"_"+var)
				allt.append(t)
				ally.append(y)
				allupt.append(upt)
				allupy.append(upy)
				variables.append(var)
				maxdate.append(md)
				
			scatter_plot(folder, title, allt, ally, allupt, allupy, labels, variables, maxdate, sharex=True)
		else:
			"No photometry was found for", self.name
		
def scatter_plot(folder, title, allt, ally, allupt, allupy,labels, variables, maxdate, sharex=False):
	"""Produces scatter plots for a given set of data.
	""" 	
	fig = plt.figure()
	plt.suptitle(title)
	npoints = len(labels)
	nrows = int(float(npoints)/2.) + (npoints % 2)

	if npoints > 0:
		for i in range(0, npoints):
			t = allt[i]
			y = ally[i]
			upt = allupt[i]
			upy = allupy[i]
			
			label = labels[i]
			var = variables[i]
			
			md = maxdate[i]
			
			ax = plt.subplot(nrows,2, i+1)
			
			if md[0] != None:
				plt.axvline(md[0], label="max date")
				x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=md[0])
				ax.xaxis.set_major_formatter(x_formatter)
			
			if md[1] != None:
				plt.axvline(md[1], color = "r", label="max vis date")
			
			if len(upt) > 0:
				arrow = u'$\u2193$'
				plt.scatter(upt, upy, marker=arrow, s=300, label="upper limits")
			
			plt.scatter(t, y, label="measurements", color="orange")

			plt.title(label)

			if var == "magnitude":
				plt.gca().invert_yaxis()
				
			elif var == "luminosity":
				ax.set_yscale('log')

			
			
			plt.xlabel("Time (MJD)")
			plt.ylabel(var)
			plt.legend()
		
		fig.set_size_inches(25, nrows*5)
		fig.subplots_adjust(hspace=.5)
		
		path = "graphs/" + folder + title + ".pdf"
		
		plt.savefig(path)
		plt.close()
		print "Saving to", path
		
	else:
		print "Not Saved!"
		print allt, ally, labels, variables
		
def histograms_plot(titles, xlabels, values):
	"""Plots histograms for a given set of variables
	"""
	fig = plt.figure()
	npoints = len(titles)
	nrows = int(float(npoints)/2.) + (npoints % 2)

	if npoints > 0:
		for i in range(0, npoints):
			ax = plt.subplot(nrows,2, i+1)
			n, bins, _ = plt.hist(values[i], bins=20, label=xlabels[i])
	
			plt.title(titles[i])
			plt.xlabel(xlabels[i])
			
	fig.set_size_inches(25, nrows*5)
	fig.subplots_adjust(hspace=.5)
	
	path = "graphs/histogram.pdf"
	plt.savefig(path)
	print "Saving to", path
	plt.close()
				
def get_degrees(time_str):
	"""Converts a HH:MM:SS time string to a float (number of hours).
	"""
	if time_str.count(':') >1 :
		s, arcm, arcs = time_str.split(':')
		if float(s) < 0.0:
			return -(abs(float(s)) + (float(arcm)/60.)+(float(arcs)/3600.))
		else:
			return (float(s) + (float(arcm)/60.)+(float(arcs)/3600.))
	else:
		arcm, arcs = time_str.split(':')
		if float(arcm) < 0.0:
			return -(abs(float(arcm)/60.)+(float(arcs)/3600.))
		else:
			return ((float(arcm)/60.)+(float(arcs)/3600.))