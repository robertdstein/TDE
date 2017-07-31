import argparse
import import_data as i
from classes import *

parser = argparse.ArgumentParser(description='Toggle online updating on or off')
parser.add_argument("-u", "--update", action="store_true")

cfg = parser.parse_args()

savepath = "tde_cat.zip"
varpath = "variables.csv"

plot = False

if cfg.update:
	i.run(savepath)

dataset = Full_set(savepath)
#dataset.plot_skymap()
#dataset.plot_spectra()
#dataset.plot_all_candidates()
#dataset.plot_all_bands()
#dataset.plot_distibutions()
#~ dataset.info_candidate("ASASSN-14ae")