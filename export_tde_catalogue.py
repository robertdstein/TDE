
# coding: utf-8

# In[16]:

import numpy as np
#from classes import Full_set

root = "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/Catalogue/"


def run(data):
	print data

	ra = []
	dec = []
	dist = []
	ddate = []
	names = []

	variablenames = ["ra_deg", "dec_deg", "lumdist", "mjddisc", "name"]

	for name in vars(data.TDEs):
		tde = getattr(data.TDEs, name)
		include = True
		for v in variablenames:
			if not hasattr(tde, v):
				include = False

		if include:
			ra.append(tde.ra_deg)
			dec.append(tde.dec_deg)
			dist.append(tde.lumdist)
			ddate.append(tde.mjddisc)
			names.append(name)

	n_sources = len(ra)

	sources = np.empty((n_sources), dtype=[("ra", np.float), ("dec", np.float),
		("flux", np.float), ("n_exp", np.float), ("weight", np.float),
		("weight_acceptance", np.float), ("weight_time", np.float),
		("weight_distance", np.float), ("norm_time", np.float),
		("global_weight_norm_time", np.float), ("discoverydate_mjd", np.float),
		("distance", np.float), ('name', 'a30'),
		])

	sources['ra'] = np.deg2rad(ra)
	sources['dec'] = np.deg2rad(dec)
	sources['distance'] = np.ones_like(sources['ra'])
	sources['flux'] = np.ones_like(sources['ra'])
	Norm = n_sources * 1.e-9 / np.sum(sources['flux'])
	sources['flux'] = sources['flux'] * Norm
	sources['weight'] = np.ones_like(sources['ra'])

	sources['discoverydate_mjd'] = np.array(ddate)
	sources['name'] = names
	path = root + "full_TDE_catalogue.npy"
	np.save(path, sources)

	print "Exporting catalogue to", path
