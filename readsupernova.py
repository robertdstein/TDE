# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 11:29:14 2017

@author: steinrob
"""

path = "supernova/"

bands =["tUBVRIJHK", "tUBVRI", "tV"]

extensions = ["sn1a_lc.v1.2.dat", "sn91t_lc.v1.1.dat", "sn91bg_lc.v1.1.dat", "sn1bc_lc.v1.1.dat", "hyper_lc.v1.2.dat", "sn2p_lc.v1.2.dat", "sn2l_lc.v1.2.dat", "sn2n_lc.v2.1.dat"]
#extensions = ["sn1a_lc.v1.2.dat", "sn91t_lc.v1.1.dat", "sn91bg_lc.v1.1.dat", "sn1bc_lc.v1.1.dat", "hyper_lc.v1.2.dat", "sn2p_lc.v1.2.dat", "sn2l_lc.v1.2.dat"]

names = ["Type Ia (Normal)", "Type Ia (1991T-like)", "Type-Ia (1991bg-like)", "Type Ib/c", "Type Ib/c High Velocity (Hypernova)", "Type IIP", "Type IIL", "Type IIn"]
bandindexes=[0, 1, 1, 1, 2, 2, 2, 2]
to_plot=[0, 1, 2, 3, 4, 5, 6]
to_plot=[0, 3]

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
		
#run()