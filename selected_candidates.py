# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 17:58:12 2017

@author: steinrob
"""

def return_candidates_with_peaks():
	names = ["PS1-11af", "PS1-10jh", "PTF09ge"]
	bands = ["r", "g", "i"]
	return names, bands
	
def xray_names():
	xrbands = ["X-Ray (0.3-2.0) keV", "X-Ray (0.3-10) keV"]
	varnames = ["luminosity", "countrate"]
	return xrbands, varnames