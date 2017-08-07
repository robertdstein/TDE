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
from matplotlib.markers import CARETDOWN
import datetime as d
import os.path
import lightcurves as lc
from tqdm import tqdm
from sklearn.externals import joblib
import plotting as p

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
	exponent = (solarabsmag - absmag)/2.5
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
		self.update_time = datetime.datetime.now()
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
		self.list_best()
		
	def creationtime(self):
		"""Returns the time that the online catalogue was last updated
		"""
		print "Catalogue version created:", self.creation_time, "\n"
		print "Dataset last expanded or updated on:", self.update_time, "\n"
		
	def info_candidate(self, name):
		"""Plots all available light curves for a given candidate.
		"""
		while not hasattr(self.TDEs, name):
			print "No TDE with this name was found!"
			print "Availablie candidates are:"
			print self.TDEs.__dict__.keys()
			print "\n"
			name=raw_input("Please reenter candidate name: \t")
		tde = getattr(self.TDEs, name) 
		tde.plot_candidate(fit=True)
		self.update_time = datetime.datetime.now()
				
	def plot_all_candidates(self, fit=False):
		"""Plots all available light curves for each candidate.
		"""
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			tde.plot_candidate(fit)
			if fit:
				self.update_time = datetime.datetime.now()
			
	def plot_all_bands(self):
		"""Plots all available light curves for each band.
		"""
		for band in self.all_available_bands:
			title = band
			
			folder = "bands/"
			combifolder="combinedbands/"
			
			allresults = []
			combiresults=[]
			
			for name in vars(self.TDEs):
				tde = getattr(self.TDEs, name)
				if tde.has_photometry:
					if band in tde.bandlist:
						i = tde.bandlist.index(band)
						
						var = tde.varlist[i]
						res= tde.return_photometry_dataset(band, var)
						
						res["variable"]=var
						res["label"]=name
						
						allresults.append(res)
						if (var == "magnitude") & (len(res["time"]) > 0):
							combiresults.append(res)
				
			p.scatter_plot(folder, title, allresults)
			p.scatter_plot(combifolder, title, combiresults, combined=True, tshift=True)
		
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
		
	def list_fits(self):
		print "The following candidates have been fitted: \n"
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			if len(tde.fits) > 0:
				print tde.name, ":"
				for d in tde.fits:
					print "\t", d,"\t", list(tde.fits[d])
		
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
		p.histograms_plot(titles, xlabels, values)
		
	def plot_fit_parameters(self, names=None, bands=None, savepath="graphs/optical_fits/combined/all.pdf"):
		""" Plots binned distributions of candidates for each variable in 'variables'.
		"""
		full = []
		
		for i in range(len(lc.return_parameter_bounds())):
			full.append([])
			
		if names == None:
			names = vars(self.TDEs)
			
		for name in names:
			tde = getattr(self.TDEs, name)
			if len(tde.fits) > 0:
				path = "graphs/optical_fits/" + tde.name+ "/"
				if not os.path.exists(path):
					os.mkdir(path)
				filename = path+"histogram.pdf"
				data = tde.fit_parameter_histogram(filename, bands)
				for j in range(len(data)):
					full[j].extend(data[j])
					
		titles, xlabels = lc.return_parameter_names()
		print "\n"
		print "Parameters \t Max \t \t Min \n"
		
		for i, param in enumerate(xlabels):
			print param, "\t \t", "{0:.3}".format(min(full[i])), "\t \t","{0:.3}".format(max(full[i]))
		
		print "\n"
		p.histograms_plot(titles, xlabels, full, savepath)
		
	def plot_fit_parameters_with_peaks(self):
		names = ["PS1-11af", "PS1-10jh", "PTF09ge"]
		bands = ["r", "g", "i"]
		savepath="graphs/optical_fits/combined/with_peak.pdf"
		self.plot_fit_parameters(names, bands, savepath)
				
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
		self.fits=dict()
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
#		self.swiftdata()
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
			metadata = dict()
			print "Found file:", datapath, 
			f = open(datapath)
			line1 = f.readline().split(" ")
			
			#Finds Swift Trigger ID
			metadata["triggerID"] = int(line1[7][:-1])
			
			#Finds swift name, and check that this matches the name of the candidate
			#This ensures an incorrectly-named file, downoladed from swift, cannot be loaded by the wrong candidate
			catname=str(self.name)
			
			mdname = str(line1[9][:-1])			
			metadata["name"] = mdname
			
			for x in [" ", "-"]:
				mdname = mdname.replace(x, "")
				catname= catname.replace(x, "")
				
			if mdname.lower() in catname.lower():
				print "\t Successful name match of", metadata["name"], "and", self.name
			else:
				print "\t", line1, metadata
				print "\n \n"
				print "WARNING! The name found in this file's metadata,", metadata["name"], ", does not match the name of the candidate,", self.name
				print "\n \n"
				print "Continue anyway?"
				choice = raw_input("[y/n]").lower()
				if choice =="y":
					pass
				else:
					raise Exception("Mismatch of candidate name (from TDE Catalogue) and metadata name (from Swift)")
			
			f.readline()

			#Set trigger time
			#All times in dataset are given relative to the time of trigger
			line3 = f.readline().split(" ")
			metadata["year"] = line3[11]
			metadata["month"] = line3[12]
			metadata["day"]= line3[13]
			metadata["time"] = line3[15]+"000"
			date = d.datetime.strptime(metadata["year"]+metadata["month"]+ metadata["day"]+metadata["time"], '%Y%b%d%X.%f')		
			triggertime =Time(date)
			metadata["triggertime"] = triggertime
			
			#Records the energy range
			line4 = f.readline().split(" ")
			metadata["e_range"]=line4[4]
			metadata["e_units"]=line4[5]
			
			for i in range (4):
				f.readline()
				
			line9  = f.readline().split("\t")
			f.readline()
			entries=[]
			i=0
			j=0
			for unsplit in f.readlines():
				entry=dict()
				row = unsplit.split("\t")
				if len(row) > 6:
					for k in range(len(line9)):
						var=line9[k]
						val=row[k]
						entry[var.replace(" ", "")]=float(val)
					
					tst = datetime.timedelta(seconds=float(entry["!Time"]))
					entry["tst"]=tst				
					isot =date+tst
					entry["isotime"] = isot
					astrotime = Time(isot)				
					entry["time"] = astrotime.mjd
					entry["u_time"] = "MJD"
					entry["telescope"]= "Swift"
					entry["instrument"]="XRT"
					entry["countrate"]=entry["Rate"]
					entry["mode"]="PC"
					entry["energy"] = metadata["e_range"].split("-")
					entry["u_energy"] = metadata["e_units"]
					if float(entry["Ratepos"]) == 0.:
						entry["upperlimit"]=True
						j+=1
					else:
						entry["e_upper_countrate"] = entry["Ratepos"]
						entry["e_lower_countrate"] = np.abs(entry["Rateneg"])
						i +=1
					entries.append(entry)
			self.swiftmetadata = metadata
			
			band = "X-Ray ("+ line4[4] +") " + str(entry["u_energy"].replace("\n",""))
			variable="countrate"
			
			self.bandlist.append(band)
			self.varlist.append(variable)
			self.nbands += 1
			setattr(self.photometry, band+"_"+variable, entries)
			
	
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
#		self.fit_xray()
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
				
	def fit_data_directly(self, var, band, ax):
		fitdata = self.dataset_to_fit(band, var)						
		if "t" in fitdata.keys():			
			if len(fitdata["t"]) > 2:
				self.plot_llh_minimisation(fitdata, ax)
				
	def plot_llh_minimisation(self, fitdata, ax):
		fitt = np.array(fitdata["t"])
		fity = np.array(fitdata["y"])
		
		uncertaintyfrac = 0.1
		
		fitup = np.array(fitdata["up"])+uncertaintyfrac
		fitdown=np.array(fitdata["down"])+uncertaintyfrac
		fit_ul_t=np.array(fitdata["ul_t"])
		
		fit_ul_y=np.array(fitdata["ul_y"])

		if fitdata["var"] == "magnitude":
			fity = -fity
			plt.ylabel(fitdata["var"])
			maxlum = max(fity)
			maxtime = self.mjdvismax
			ax.axvline(maxtime, label="Peak Visual Magnitude",color="r", linestyle="--")
		else:
			plt.ylabel("Log(" +  fitdata["var"]+")")
			maxlum = max(fity)
			mask = fity == maxlum
			maxtime = fitt[mask]
			
		fitsig = np.array(zip(fitup, fitdown))
		
		ax.errorbar(fitt, fity, yerr=[fitup, fitdown], label="measurements", color="b", ecolor="r", fmt='o')
		
		plt.xlabel("t")
		
		params=[]
		covars = []
		vals = []
#		
		offset = np.linspace(-30.0, 120., num=6)
		A0vals = [maxlum]
		avals = [10**-3]
		ttvals = [10., 100.]
		cvals = [-0.2, -1.0, -2.0, -4.0]
		fvals = [0.0]
		N = len(avals)*len(ttvals)*len(cvals)*len(fvals)*len(offset)*len(A0vals)
		
		def llh(p, to_print=False):
			ll=0.0
	
			for k, t in enumerate(fitt):
				time = t - maxtime + p[5]
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
				time = ul_t - maxtime + p[5]
				y = fit_ul_y[k]
				model = lc.fitfunc(time, p)
				if y < model:
					ll += -np.log(0.05)				
			return ll
		
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
								for f in fvals:
									pinit = [A0, a, tt, c, f, i]
									bnds = lc.return_parameter_bounds(A0)
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
		
		self.fits[fitdata["band"]] = bestparams
		covar = covars[index]
		
		t_max_to_peak=bestparams[5]
		
		fittedpeaktime = maxtime-t_max_to_peak
		
		plottimes = np.linspace(fitt[0]-100, fitt[-1]+100, 5000)		
		shiftedtimes = plottimes-fittedpeaktime
				
		def fit(times):
			models=[]
			for t in times:
				models.append(lc.fitfunc(t, bestparams))
			return models
		
		ax.annotate(str(bestparams) +" " + str(best), xy=(0.05, 0.1), xycoords="axes fraction")

		ax.plot(plottimes, fit(shiftedtimes), label="Fit")
		
		ax.axvline(fittedpeaktime -t_max_to_peak + bestparams[2], label="Model Transition",color="purple", linestyle="--")
		
		plt.ylabel("Magnitude")
		
		if len(fit_ul_t) > 0:
			ax.scatter(fit_ul_t, -fit_ul_y, marker=CARETDOWN, s=150, label="Upper Limits")
			
		lowerlim = min(fity)-1
		upperlim = max(fity) + 1
#		
		ax.set_ylim(bottom=lowerlim, top=upperlim)
		ax.set_xlim(fitt[0]-(100+np.abs(t_max_to_peak)))
		ax.legend()		
		
		plt.legend()

	def fit_xray(self):
		"""Checks to see if there is a light curve visible in X Ray measurements.
		Fits this distribution to see if it is close to the charecteristic t^-5/3.
		"""
		xrband = "X-Ray (0.3-2.0) keV"
		var = "luminosity"
		if hasattr(self, "bandlist"):
			if xrband in self.bandlist:
				fitdata = self.dataset_to_fit(xrband, var)						
				if "t" in fitdata.keys():			
					if len(fitdata["t"]) > 2:
						fig = plt.figure()	
						ax = plt.subplot(111)
						plt.suptitle(self.name)
						self.plot_llh_minimisation(fitdata, ax)
				
						fig.set_size_inches(10, 7)						
				
						path = "graphs/xrays/" + self.name+ ".pdf"
						plt.savefig(path)
						plt.close()
						print "Saving to", path
				
	def fit_optical(self):
		"""Checks to see if there is a light curve visible in X Ray measurements.
		Fits this distribution to see if it is close to the charecteristic t^-5/3.
		"""
		var="magnitude"
		path = "graphs/optical_fits/" + self.name+ "/"
		if hasattr(self, "bandlist"):
			for band in self.bandlist:
				if hasattr(self.photometry, band+"_"+var):
					fitdata = self.dataset_to_fit(band, var)						
					if "t" in fitdata.keys():			
						if len(fitdata["t"]) > 2:
							fig = plt.figure()	
							ax = plt.subplot(111)
							plt.suptitle(self.name)
							
							print "Minimising for candidate", self.name, "in band", band							
							
							self.plot_llh_minimisation(fitdata, ax)
					
							fig.set_size_inches(10, 7)						

							if not os.path.exists(path):
								os.mkdir(path)
							file_name = path + band + ".pdf"
							plt.savefig(file_name)
							plt.close()
							print "Saving to", file_name

		filename = path + "histogram.pdf"
		if len(self.fits) > 0:
			self.fit_parameter_histogram(filename)				
							
	def fit_parameter_histogram(self, path, bands=None):
		"""Produces histograms to show the distribution of fit parameters
		"""
		data = self.fits
		ndata = len(data[data.keys()[0]])
		
		if bands == None:
			bands = data.keys()
		
		lists=[]
		
		titles, xlabels = lc.return_parameter_names()
		
		if len(xlabels) != ndata:
			raise Exception("Number of parameters does not match number of labels!")
		
		for i in range(ndata):
			lists.append([])
			
		for band in bands:
			if band in data.keys():
				entry = data[band]
				for i in range(len(entry)):
					lists[i].append(entry[i])			
		
		if len(lists[0]) > 0:
			p.histograms_plot(titles, xlabels, lists, path, bnds=True)
		return lists
				
	def dataset_to_fit(self, band, var):
		"""Checks to see if there is a light curve visible in swift X Ray measurements.
		Fits this distribution to see if it is close to the charecteristic t^-5/3.
		"""
		fitdata = dict()
		fitdata["var"]=var
		fitdata["band"]=band

		if hasattr(self, "bandlist"):
			if band in self.bandlist:
				
				results = self.return_photometry_dataset(band, var)
									
				if len(results["y"]) > 0:					
					fitdata["t"]=[]
					fitdata["y"]=[]
					fitdata["up"]=[]
					fitdata["down"]=[]
					fitdata["ul_t"]=results["upper_limit_time"]
					fitdata["ul_y"]=results["upper_limit_y"]
					
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
								
#					allt = results["time"]+results["upper_limit_time"]
					
#					fitdata["last_t"]=0.0			
#					
#					for j in range(len(allt)):
#						t = allt[j]
#						if fitdata["last_t"] < t < fitdata["maxtime"]:
#							fitdata["last_t"]=t
		
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
			
			results = dict()
			results["time"]=[]
			results["y"]=[]
			results["upper_limit_time"] = []
			results["upper_limit_y"] =[]
			results["sigma_upper"]=[]
			results["sigma_lower"]=[]											
			
			if hasattr(self, "mjdmax"):
				results["maxdate"] = [self.mjdmax]
			else:
				results["maxdate"]  = [None]
				
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
						
						if "e_upper_"+var in entry.keys():
							sigup = float(entry["e_upper_" + var])
							sigdown = float(entry["e_lower_" + var])
						elif "e_"+var in entry.keys():
							sigup = float(entry["e_" + var])
							sigdown=sigup
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
			opticalfolder="optical/"
			title = self.name
			
			allresults=[]
			opticalresults=[]
			
			for band in sorted(self.bandlist, key=str.lower):
				i = self.bandlist.index(band)
				
				var = self.varlist[i]
				res = self.return_photometry_dataset(band, var)
				
				res["variable"]=var
				res["label"]=band+"_"+var
				res["band"]=band
				
				allresults.append(res)
				if var == "magnitude":
					opticalresults.append(res)
				
			p.scatter_plot(folder, title, allresults)
			if len(opticalresults) > 0:
				p.scatter_plot(opticalfolder, title, opticalresults, combined=True)
		else:
			"No photometry was found for", self.name
						
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
