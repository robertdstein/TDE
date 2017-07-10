import zipfile
import json
import datetime
from pprint import pprint
import matplotlib
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
		self.problems=[]
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
		#~ print vars(allbands)
		for band in vars(allbands):
			if len(getattr(self.photometry, band)) > minpoints:
				self.has_photometry = True
				
	def return_photometry_dataset(self, band, var):
		if hasattr(self.photometry, band+"_"+var):
			data = getattr(self.photometry,band+"_"+var)
			
			time=[]
			y=[]
			upt = []
			upy =[]
			
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
				else:
					print "Time not in entry for some reason...."
					print entry
			return time, y, var, upt, upy
		else:
			return None, None, None
			
	def plot_candidate(self):
		if self.has_photometry:
			print "Producing individualised graph for", self.name, "covering the following bands:"
			print self.bandlist
			folder = "candidates/"
			title = self.name
			
			allt = []
			ally = []
			allupt = []
			allupy = []
			labels =[]
			variables = []
			
			for band in sorted(self.bandlist, key=str.lower):
				i = self.bandlist.index(band)
				
				var = self.varlist[i]
				t, y, var, upt, upy = self.return_photometry_dataset(band, var)
				
				labels.append(band+"_"+var)
				allt.append(t)
				ally.append(y)
				allupt.append(upt)
				allupy.append(upy)
				variables.append(var)
				
			plot(folder, title, allt, ally, allupt, allupy, labels, variables, sharex=True)
		else:
			"No photometry was found for", self.name		
		
def plot(folder, title, allt, ally, allupt, allupy,labels, variables, sharex=False):
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
			
			if sharex and i > 0:
				plt.subplot(nrows,2, i+1, sharex = ax)
			else:
				ax = plt.subplot(nrows,2, i+1)
				if len(t) > 0:
					starttime = t[0]
				else:
					starttime = upt[0]
				if sharex:
					pass
				else:
					x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=starttime)
					ax.xaxis.set_major_formatter(x_formatter)
					ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(11))
					
			plt.title(label)

			if var == "magnitude":
				plt.gca().invert_yaxis()

			plt.scatter(t, y)
			if len(upt) > 0:
				plt.scatter(upt, upy, marker="x")
			
			plt.xlabel("Time (MJD)")
			plt.ylabel(var)
		plt.gca().set_xlim(left=(starttime - 20))
		fig.set_size_inches(25, nrows*5)
		fig.subplots_adjust(hspace=.5)
		
		path = "graphs/" + folder + title + ".pdf"
		
		plt.savefig(path)
		plt.close()
		print "Saved to", path
		
	else:
		print "Not Saved!"
		print allt, ally, labels, variables
	
	
	
		
		
		
