import zipfile
import json
import datetime
import matplotlib
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt

minpoints = 3


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
		if hasattr(self.TDEs, name):
			tde = getattr(self.TDEs, name) 
			print "TDE found, with photometry =", tde.has_photometry
			tde.plot_candidate()
		else:
			raise Exception("No TDE with this name was found!")
				
	def plot_all_candidates(self):
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			tde.plot_candidate()
			
	def plot_all_bands(self):
		for band in self.all_available_bands:
			title = band
			#~ print "Producing graph for band", band
			
			
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
		best=[]
		print "\n The best candidates, in terms of available data, are: \n"
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			if tde.best:
				best.append(name)
				print name
		self.best_list = best
		print "\n \n"
	
	def plot_distibutions(self):
		variables = ["redshift", "comovingdist", "lumdist", "hostoffsetdist", "mjdmax", "nbands","maxabsmag",  "nhhost", "ra_float", "dec_float"]
		
		titles = ["Redshift", "Comoving Distance", "Luminosity Distance", "Host Offset Distance", "Date of Max Luminosity", "Number of observational bands", "Max Absolute Magnitude", "Something cool, probably (nnhost)", "Right ascension", "Declination"]
		xlabels = ["Redshift (z)", "Distance (MPc)", "Distance (MPc)", "Distance (?)", "Time(MJD)", "n", "", "???", "hours", "degrees"]
		
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
		for name in vars(self.TDEs):
			tde = getattr(self.TDEs, name)
			if tde.has_spectra:
				tde.plot_spectrum()
		
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
		if hasattr(self, "ra") and hasattr(self, "dec"):			
			self.ra_float = get_degrees(self.ra.value)
			self.dec_float = get_degrees(self.dec.value)
	
	def isbest(self):
		if self.has_photometry:
			if hasattr(self, "mjdmax") and hasattr(self, "mjdvismax"):
				if self.nbands > 2.0:
					self.best=True
	
	def mjd_time(self):
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
		#~ print self.name, self.bandlist, self.varlist
		self.photometry = allbands
		self.nbands = float(len(bandlist))
		#~ print vars(allbands)
		for band in vars(allbands):
			if len(getattr(self.photometry, band)) > minpoints:
				self.has_photometry = True
				
	def spectrum(self):
		if hasattr(self, "spectra"):
			if self.spectra.u_fluxes != "Uncalibrated":
				self.has_spectra = True
			
	def plot_spectrum(self):
		
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
			return None, None, None
			
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
			
			plt.title(label)

			if var == "magnitude":
				plt.gca().invert_yaxis()
				
			elif var == "luminosity":
				ax.set_yscale('log')

			plt.scatter(t, y, label="measurements")
			
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
	
#def get_sec(time_str):
#	h, m, s = time_str.split(':')
#	return int(h) * 3600 + int(m) * 60 + int(s)
				
def get_degrees(time_str):
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
	
	
	
	
		
		
		
