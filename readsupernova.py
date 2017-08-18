# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 11:29:14 2017

@author: steinrob
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotting as p
import loglikelihoodminimisation as llh
import lightcurves as lc

path = "supernova/"

bands =["tUBVRIJHK", "tUBVRI", "tV"]

extensions = ["sn1a_lc.v1.2.dat", "sn91t_lc.v1.1.dat", "sn91bg_lc.v1.1.dat", "sn1bc_lc.v1.1.dat", "hyper_lc.v1.2.dat", "sn2p_lc.v1.2.dat", "sn2l_lc.v1.2.dat", "sn2n_lc.v2.1.dat"]
#extensions = ["sn1a_lc.v1.2.dat", "sn91t_lc.v1.1.dat", "sn91bg_lc.v1.1.dat", "sn1bc_lc.v1.1.dat", "hyper_lc.v1.2.dat", "sn2p_lc.v1.2.dat", "sn2l_lc.v1.2.dat"]

names = ["Type Ia (Normal)", "Type Ia (1991T-like)", "Type-Ia (1991bg-like)", "Type Ib+Ic", "Hypernova", "Type IIP", "Type IIL", "Type IIn"]
bandindexes=[0, 1, 1, 1, 2, 2, 2, 2]
to_plot=[0, 1, 2, 3, 4, 5, 6]
#to_plot=[0, 3]

def run():
	sn_data=[]
	for i in to_plot:
		extension=extensions[i]
		filename = path+extension
		name = names[i]
		band = bands[bandindexes[i]]
		
		data = dict()
		data["name"] = name
		data["band"] = band
		
		
		for letter in band:
			data[letter] = []
		
		with open(filename) as f:
			for unsplit in f.readlines():
				line = unsplit.split(" ")
				line = [x for x in line if (x != '') and (x != "\n")]
				if len(line) > 1:
					for j, value in enumerate(line):
						letter = band[j]
						data[letter].append(float(value))
					
		sn_data.append(data)
	return sn_data
		
def plot():
	sn_data = run()
	params=[]
	for data in sn_data:
		letters = data["band"][1:]
		for band in letters:
			fitdata = dict()
			y = np.array(data[band])
			mask = y < 5.0
			fitdata["y"] = y[mask]
			fitdata["t"] = np.array(data["t"])[mask]
			fitdata["maxtime"]=0.0
			fitdata["up"] = np.zeros_like(y[mask])
			fitdata["down"] = fitdata["up"]
			fitdata["ul_t"] =[]
			fitdata["ul_y"] =fitdata["ul_t"]
			fitdata["var"] = "magnitude"
			
			path = "graphs/optical_fits/" + data["name"]+ "/"
			fig = plt.figure()
									
			print "Minimising for candidate", data["name"], "in band", band							
									
			info = llh.run(fitdata, data["name"])
							
			fig.set_size_inches(7, 10)						
		
			if not os.path.exists(path):
				os.mkdir(path)
			file_name = path + band + ".pdf"
			plt.savefig(file_name)
			plt.close()
			print "Saving to", file_name
			params.append(info)
	
	return params
	
def fit_histogram(data):	
	savepath="graphs/optical_fits/combined/supernovae"
	
	for suffix in ["", "_loose"]:		
		
		titles, xlabels = lc.return_all_parameter_names()
		
		full=[[] for i in xlabels]
		
		for d in data:
			for j, label in enumerate(xlabels):
				full[j].append(d[label+suffix])

		print "\n"
		print "Supernovae" + suffix
		print "\n"
		print "Parameters \t Max \t \t Min \n"
	
		for j, param in enumerate(xlabels):
			print param, "\t \t", "{0:.3}".format(min(full[j])), "\t \t","{0:.3}".format(max(full[j]))
		
		print "\n"
		p.histograms_plot(titles, xlabels, full, savepath+suffix+"_")

def fit_scatter(data):
	folder = "misc/"
	title = "supernova_2d_fit_distribution"
	
	variables, variablenames, pairIDs = p.config_sd()
	
	n = [names[i] for i in to_plot]
	suffixes = ["", "_loose"]
	
	for suffix in suffixes:
		
		allvals=[]
		for pair in pairIDs:
			res=dict()
			
			letters = ["x", "y", "z"]
			allvars=[]
			
			for i, ID in enumerate(pair):
				letter = letters[i]
				res[letter+"var"]= variables[ID]
				res[letter+"_label"] = variablenames[ID]
				allvars.append(variables[ID])
			
			res["vars"] = allvars

			allvals.append(res)
		
		for res in allvals:
			x=[]
			y=[]
			z=[]
#			for name in n:
			for entry in data:
				vals=[]
				for var in res["vars"]:
					val = entry[var.split(".")[-1]+suffix]
			
				if  isinstance(val, float) or isinstance(val, list):
					vals.append(val)
					
				if len(vals) == len(res["vars"]):
					for k, v in enumerate(vals):
						if isinstance(v, list):
							[x,y,z][k].extend(v)
						elif isinstance(v, float):
							[x,y,z][k].append(v)
			res["x"] = x
			res["y"] = y
			if len(z) > 0:
				res["z"] = z
						
		p.scatter_distibution(folder, title+suffix, allvals)
