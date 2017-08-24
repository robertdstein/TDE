import argparse
import os
import import_data as i
import add_metadata as md
from classes import FullSet
from sklearn.externals import joblib
import export_tde_catalogue

parser = argparse.ArgumentParser(description='Select modues to run')
parser.add_argument("-u", "--update", action="store_true")
parser.add_argument("-md", "--metadata", action="store_true")
parser.add_argument("-f", "--fit", action="store_true")
parser.add_argument("-c", "--candidate", nargs='?', const="OGLE16aaa", type=str)
parser.add_argument("-p", "--plot", action="store_true")
parser.add_argument("-sn", "--supernova", action="store_true")
parser.add_argument("-e", "--export", action="store_true")

cfg = parser.parse_args()

root_path = "/afs/ifh.de/user/s/steinrob/Desktop/python/TDE/"

sourcepath = root_path + "tde_cat.zip"
varpath = root_path + "variables.csv"
savepath = root_path + "pickle/dataset.pkl"

print ""

if not cfg.update and os.path.exists(savepath):
    dataset = joblib.load(savepath)
else:
    i.run(sourcepath)
    if cfg.metadata:
        md.run(sourcepath)
    dataset = FullSet(sourcepath)
    joblib.dump(dataset, savepath)

dataset.creationtime()
# dataset.list_fits()

if cfg.candidate is not None:
    dataset.info_candidate(cfg.candidate)
    joblib.dump(dataset, savepath)

# dataset.add_swift()
# joblib.dump(dataset, savepath)

if cfg.fit:
    dataset.plot_all_candidates(fit=True)
    joblib.dump(dataset, savepath)

if cfg.supernova:
    dataset.plot_supernova()
    joblib.dump(dataset, savepath)

if cfg.plot:
    dataset.plot_distibutions()
    dataset.plot_all_bands()
    dataset.plot_2d_distributions()
    dataset.plot_2d_fit_distributions()
    dataset.plot_2d_fit_with_peaks()
    dataset.plot_2d_fit_xray()
    dataset.plot_fit_parameters_with_peaks()
    dataset.plot_fit_parameters()
    dataset.plot_fit_parameters_xrays()
    dataset.plot_skymap()
    dataset.plot_spectra()

if cfg.export:
    export_tde_catalogue.run(dataset)
