# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 11:29:14 2017

@author: steinrob
"""

path = "supernova/"

bands =["UBVRIJHK", "UBVRI", "V"]

extensions = ["sn1a_lc.v1.2.dat", "sn91t_lc.v1.1.dat", "sn91bg_lc.v1.1.dat", "sn1bc_lc.v1.1.dat", "hyper_lc.v1.2.dat", "sn2p_lc.v1.2.dat", "sn2l_lc.v1.2.dat", "sn2n_lc.v2.1.dat"]
names = ["Type Ia (Normal)", "Type Ia (1991T-like)", "Type-Ia (1991bg-like)", "Type Ib/c", "Type Ib/c High Velocity (Hypernova)", "Type IIP", "Type IIL", "Type IIn"]
bandindexes=[0, 1, 1, 1, 2, 2, 2, 2]

def run():
	for i, extension in enumerate(extensions):
		filename = path+extension
		name = names[i]
		band = bands[bandindexes[i]]
		print filename, name, band
		with open()
		
run()