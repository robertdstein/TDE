# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:25:22 2017

@author: steinrob
"""
from astropy.time import Time
import datetime as d
import numpy as np
import math

def run(datapath, catname):
	metadata = dict() 
	with open(datapath) as f:
		line1 = f.readline().split(" ")
		
		#Finds Swift Trigger ID
		metadata["triggerID"] = int(line1[7][:-1])
		
		#Finds swift name, and check that this matches the name of the candidate
		#This ensures an incorrectly-named file, downoladed from swift, cannot be loaded by the wrong candidate
		
		mdname = str(line1[9][:-1])			
		metadata["name"] = mdname
		
		for x in [" ", "-"]:
			mdname = mdname.replace(x, "")
			catname= catname.replace(x, "")
			
		if mdname.lower() in catname.lower():
			print "\t Successful name match of", metadata["name"], "and", catname
		else:
			print "\t", line1, metadata
			print "\n \n"
			print "WARNING! The name found in this file's metadata,", metadata["name"], ", does not match the name of the candidate,", catname
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
				
				tst = d.timedelta(seconds=float(entry["!Time"]))
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
		
		
		band = "X-Ray ("+ line4[4] +") " + str(entry["u_energy"].replace("\n",""))

		return band, metadata, entries
		
def condense(entries):
	new=[]
	t0 = 0.0
	entry = entries[0]
	times=[]
	counts=[]
	up_errors=0.0
	down_errors=0.0
	for entry in entries:
		if "upperlimit" not in entry.keys():
			if int(entry["time"]) > t0:
				if len(times) > 0:
					newentry = dict()
					newentry["time"] = np.mean(times)
					newentry["countrate"] = np.sum(counts)/float(len(counts))
					newentry["e_upper_countrate"] = math.sqrt(up_errors)/float(len(counts))
					newentry["e_lower_countrate"] = math.sqrt(down_errors)/float(len(counts))
					new.append(newentry)
					variables = ["u_energy", "energy", "instrument", "telescope"]
					for var in variables:
						newentry[var] = entry[var]
				t0 = int(entry["time"])
				times=[]
				counts=[]
				up_errors=0.0
				down_errors=0.0
			
			times.append(entry["time"])
			counts.append(entry["countrate"])
			up_errors += entry["e_upper_countrate"]**2
			down_errors += entry["e_lower_countrate"]**2
			
	print "Condensed data from", len(entries), "to", len(new)
	return new
			