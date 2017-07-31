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
				
			scatter_plot(folder, title, allresults)
			scatter_plot(combifolder, title, combiresults, combined=True, tshift=True)
		
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
				
	def fit_data_directly(self, var, band):
		fitdata = self.dataset_to_fit(band, var)						
		if "t" in fitdata.keys():			
			if len(fitdata["t"]) > 2:
				self.plot_llh_minimisation_complex(fitdata)
				
	def plot_llh_minimisation_complex(self, fitdata):		
		fig = plt.figure()	
		ax1 = plt.subplot(312)
		
		if hasattr(self, "mjdmax"):
			plt.axvline(self.mjdmax, color="k", label="max date")
			
		if hasattr(self, "mjdmax"):
			plt.axvline(self.mjdvismax, color = "green", label="max vis date")
			
		plt.suptitle(self.name)
		ax1.set_yscale('log')
		
		fitt = np.array(fitdata["t"])
		fity = np.array(fitdata["y"])
		fitup = np.array(fitdata["up"])
		fitdown=np.array(fitdata["up"])
		fitsig = np.array(zip(fitup, fitdown))
		
		maxtime=float(int(fitdata["maxtime"]))
		maxlum=fitdata["maxlum"]
		
		ax1.errorbar(fitt, fity, yerr=[fitup, fitdown], label="measurements", color="b", ecolor="r", fmt='o')
		plt.xlabel("t")
		plt.ylabel("Log(" +  fitdata["var"]+")")
		plt.legend()
		
		params=[]
		covars = []
		vals = []			
		
		offset = np.linspace(-50, 150, num=200)
		acceptedoffset=[]
		
		for i in offset:
			starttime = maxtime - i
			
			newt = fitt-starttime
				
			def llh(p):
				ll=0.0
				for k in range(len(newt)):
					time = newt[k]
					y = fity[k]
					model = lc.fitfunc(time, p)
					if y > model:
						err=fitsig[k][0]
					else:
						err=fitsig[k][1]
					ll += ((y - model)/err)**2
				return ll
	
			pinit = [maxlum, 10**-3, 10, -5./3.]
			bnds = [(0.0, None), (0.001, 0.1), (1., 250), (-2,-0.01) ]
			out = optimize.minimize(llh, pinit, bounds=bnds, tol=10**-10)		
			
			pfinal = out.x
			covar = out.hess_inv
			covars.append(covar)
			params.append(pfinal)
			vals.append(np.sum(out.fun))
			acceptedoffset.append(i)
		
		ax2 = plt.subplot(311)
		
		best = min(vals)
		index = vals.index(best)
		time = maxtime - (acceptedoffset[index])
		
		bestparams = params[index]
		covar = covars[index]
		
		print bestparams

		lower=acceptedoffset[0]
		upper=acceptedoffset[-1]
		mid=False
			
		lowerindex=index
		upperindex=index
		lower=best
		upper = best
		
		while (lower < best+0.5) and (lowerindex >0):
			lowerindex-=1
			lower = vals[lowerindex]
		
		while (upper < best+0.5) and (upperindex +1 <len(vals)):
			upperindex += 1
			upper = vals[upperindex]
		
#		for j in range(len(vals)):
#			v = vals[j]
#			if not mid:
#				if v > (0.5+best):
#					lower = acceptedoffset[j]
#				else:
#					mid = True
#			elif v > (0.5+best):
#				upper = acceptedoffset[j]
#				break	
						
		plt.plot(offset, np.array(vals)-best)
		plt.xlabel("Log(Age of TDE at observed maximum brightness)")
		plt.ylabel(r"$\Delta$ Log Likelihood (Gaussian)")
#		ax2.set_xscale("log")
		
		plt.axvspan(acceptedoffset[lowerindex], acceptedoffset[upperindex], color="r", alpha=0.3)
		plt.scatter(acceptedoffset[index], 0.0, color="r", marker="*")
		plt.axhline(0.5, color="g")

		plottimes = np.linspace(fitt[0]-100, fitt[-1], 200)
		shiftedtimes = plottimes-time
				
		def fit(t):
			models=[]
			for k in range(len(t)):
				time = t[k]
				models.append(lc.fitfunc(time, bestparams))
			return models
		
		plt.title("Minimum occurs at " + str(int(time)) + " MJD (" + str(int(acceptedoffset[index]))+" days before date of observed maximum, with 1 " r"$\sigma$ range of " + str(int(acceptedoffset[lowerindex])) + " to "+ str(int(acceptedoffset[upperindex])) + " days before maximum")						
		ax3 = plt.subplot(313)
		
#		if type(covar) != type(None):							
#			ax1.annotate("Best fit has power law with gradient " + "{0:.3g}".format(bestparams[1]) + " +/- " + str(math.sqrt(np.abs(covar[1][1]))), xy=(0.05, 0.1), xycoords="axes fraction")
#			ax3.annotate("Best fit maximum luminosity is " "{0:.3g}".format(bestparams[0]) + " +/- " + str(math.sqrt(np.abs(covar[0][0]))), xy=(0.05, 0.1), xycoords="axes fraction")
		if True:
			ax1.annotate("Best fit has power law with gradient " + "{0:.3g}".format(bestparams[1]) + " but failed to find covariance (fit did not converge correctly)", xy=(0.05, 0.1), xycoords="axes fraction")
			ax3.annotate("Best fit maximum luminosity is " "{0:.3g}".format(bestparams[0]) + " but failed to find covariance (fit did not converge correctly)", xy=(0.05, 0.1), xycoords="axes fraction")
		
		newt = fitt-time

		ax3.errorbar(newt, fity, yerr=[fitup, fitdown], label="measurements", color="b", ecolor="r", fmt='o')
		ax3.plot(shiftedtimes, fit(shiftedtimes)) 
		ax1.plot(plottimes, fit(shiftedtimes))
		plt.ylabel("Log(Luminosity)")
		ax1.set_ylim(bottom=maxlum*10**-3)
		ax3.set_ylim(bottom=maxlum*10**-3)
		plt.xlabel("Log(t-"+ str(int(time))+" MJD)")
		
		
		ax3.set_yscale('log')
		ax3.set_xscale('log')
			
		fig.set_size_inches(15, 20)						
		
		path = "graphs/xrays/" + self.name+ ".pdf"						
		
		plt.savefig(path)
		plt.close()
		print "Saving to", path
		
#	def plot_llh_minimisation(self, fitdata):		
#		fig = plt.figure()	
#		ax1 = plt.subplot(312)
#		
#		if hasattr(self, "mjdmax"):
#			plt.axvline(self.mjdmax, color="k", label="max date")
#			
#		if hasattr(self, "mjdmax"):
#			plt.axvline(self.mjdvismax, color = "green", label="max vis date")
#			
#		plt.suptitle(self.name)
#		ax1.set_yscale('log')
#		
#		fitt = np.array(fitdata["t"])
#		fity = np.array(fitdata["y"])
#		fitup = np.array(fitdata["up"])
#		fitdown=np.array(fitdata["down"])
#		fitsig = np.array(zip(fitup, fitdown))
#		
#		maxtime=float(int(fitdata["maxtime"]))
#		maxlum=fitdata["maxlum"]
#		
#		ax1.errorbar(fitt, fity, yerr=[fitup, fitdown], label="measurements", color="b", ecolor="r", fmt='o')
#		plt.xlabel("t")
#		plt.ylabel("Log(" +  fitdata["var"]+")")
#		plt.legend()
#		
#		params=[]
#		covars = []
#		vals = []			
#		
#		offset = np.linspace(-150, 300, num=100)
#		acceptedoffset=[]
#		
#		for i in offset:
#			starttime = maxtime - (i + 1e-5)
#			
#			newt = fitt-starttime
#			check = newt>(0.0)
#			if len(newt[check]) > 2:
#			
#				def indices(y,model):
#					if y > model:
#						return 0
#					else:
#						return 1
#						
#				vfunc = vectorize(indices)
#				
#				def fitfunc(p, x):
#					return (p[0]*(x**p[1]))
#					
#				def llh(p, x, y, err2):
#					model = fitfunc(p, x)
#					err_indices = vfunc(y, model)
#					err = []
#					for j in range(len(err_indices)):
#						err.append(err2[j][err_indices[j]])
#					err= np.array(err)
#					llh = ((y - model)/err)**2
#					return llh
#	
#				pinit = [max(fity[check]), -1.0]
#				out = optimize.leastsq(llh, pinit,
#				                       args=(newt[check], fity[check], fitsig[check]), full_output=1)
#				
#				pfinal = out[0]
#				covar = out[1]
#				covars.append(covar)
#				params.append(pfinal)
#				vals.append(np.sum(llh(pfinal, newt[check], fity[check],fitsig[check])))
#				acceptedoffset.append(i)
#		
#		ax2 = plt.subplot(311)
#		
#		best = min(vals)
#		index = vals.index(best)
#		time = maxtime - (acceptedoffset[index])
#
#		lower=acceptedoffset[0]
#		upper=acceptedoffset[-1]
#		mid=False
#		
#		for j in range(len(vals)):
#			v = vals[j]
#			if not mid:
#				if v > (0.5+best):
#					lower = acceptedoffset[j]
#				else:
#					mid = True
#			elif v > (0.5+best):
#				upper = acceptedoffset[j]
#				break
#						
#		plt.plot(acceptedoffset, np.array(vals)-best)
#		plt.xlabel("Log(Age of TDE at observed maximum brightness)")
#		plt.ylabel(r"$\Delta$ Log Likelihood (Gaussian)")
##		ax2.set_xscale("log")
#		
#		plt.axvspan(lower, upper, color="r", alpha=0.3)
#		plt.scatter(acceptedoffset[index], 0.0, color="r", marker="*")
#		plt.axhline(0.5, color="g")
#
#		bestparams = params[index]
#		covar = covars[index]
#
#		plottimes = np.linspace(fitt[0], fitt[-1], 200)
#		shiftedtimes = plottimes-time
#		
#		def fit(x):
#			return bestparams[0] * (x**bestparams[1])
#		
#		plt.title("Minimum occurs at " + str(int(time)) + " MJD (" + str(int(acceptedoffset[index]))+" days before date of observed maximum, with 1 " r"$\sigma$ range of " + str(int(lower)) + " to "+ str(int(upper)) + " days before maximum")						
#		ax3 = plt.subplot(313)
#		
#		if type(covar) != type(None):							
#			ax1.annotate("Best fit has power law with gradient " + "{0:.3g}".format(bestparams[1]) + " +/- " + str(math.sqrt(np.abs(covar[1][1]))), xy=(0.05, 0.1), xycoords="axes fraction")
#			ax3.annotate("Best fit maximum luminosity is " "{0:.3g}".format(bestparams[0]) + " +/- " + str(math.sqrt(np.abs(covar[0][0]))), xy=(0.05, 0.1), xycoords="axes fraction")
#		else:
#			ax1.annotate("Best fit has power law with gradient " + "{0:.3g}".format(bestparams[1]) + " but failed to find covariance (fit did not converge correctly)", xy=(0.05, 0.1), xycoords="axes fraction")
#			ax3.annotate("Best fit maximum luminosity is " "{0:.3g}".format(bestparams[0]) + " but failed to find covariance (fit did not converge correctly)", xy=(0.05, 0.1), xycoords="axes fraction")
#		
#		newt = fitt-time
#		check = newt>(0.0)
#
#		ax3.errorbar(newt[check], fity[check], yerr=[fitup[check], fitdown[check]], label="measurements", color="b", ecolor="r", fmt='o')
#		ax3.plot(shiftedtimes, fit(shiftedtimes)) 
#		ax1.plot(plottimes, fit(shiftedtimes))
#		plt.ylabel("Log(Luminosity)")
#		plt.xlabel("Log(t-"+ str(int(time))+" MJD)")
#		
#		
#		ax3.set_yscale('log')
#		ax3.set_xscale('log')
#			
#		fig.set_size_inches(15, 20)						
#		
#		path = "graphs/xrays/" + self.name+ ".pdf"						
#		
#		plt.savefig(path)
#		plt.close()
#		print "Saving to", path
				
	def check_xray(self):
		"""Checks to see if there is a light curve visible in X Ray measurements.
		Fits this distribution to see if it is close to the charecteristic t^-5/3.
		"""
		xrband = "X-Ray (0.3-2.0) keV"
		var = "luminosity"
		if hasattr(self, "bandlist"):
			if xrband in self.bandlist:
				self.fit_data_directly(var, xrband)
				
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
					fitdata["maxlum"] = float(max(results["y"]))
					fitdata["maxtime"] = results["time"][results["y"].index(fitdata["maxlum"])]
					
					fitdata["t"]=[]
					fitdata["y"]=[]
					fitdata["up"]=[]
					fitdata["down"]=[]
					for i in range(len(results["time"])):
						t = results["time"][i]
						if float(t) >= fitdata["maxtime"]:
							fitdata["t"].append(float(t))
							fitdata["y"].append(float(results["y"][i]))
							if results["sigma_upper"][i] > 0:
								fitdata["up"].append(results["sigma_upper"][i])
								fitdata["down"].append(results["sigma_lower"][i])
							else:
								fitdata["up"].append(float(results["y"][i])*0.2)
								fitdata["down"].append(float(results["y"][i])*0.2)
			
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
			
	def plot_candidate(self):
		if self.has_photometry:
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
				
			scatter_plot(folder, title, allresults)
			if len(opticalresults) > 0:
				scatter_plot(opticalfolder, title, opticalresults, combined=True)
			
		else:
			"No photometry was found for", self.name
		
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
