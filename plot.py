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
dataset.plot_spectra()
#dataset.plot_all_candidates()
dataset.plot_all_bands()
dataset.plot_distibutions()
#~ dataset.info_candidate("ASASSN-14ae")


#~ if plot:
	#~ with open(varpath, 'rb') as csvfile:
		#~ reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		#~ for row in reader:
			
			#~ key, variable, channel, value = row
			
			#~ print "Finding", variable, "(", key, ") with requirement", channel, "=", value
	
			#~ accepted =[]
			#~ nreject = 0
			
			#~ alltimes=[]
			#~ allydata=[]
			
			#~ for name in filelist:
				#~ with zf.open(name) as f:  
					#~ data = f.read()
					#~ d = json.loads(data)
					
					#~ print d.values()
					
					#~ fullset =  d.values()[0]
					
					#~ time=[]
					#~ ydata=[]
					
					#~ passed = False
					
					#~ if key in fullset.keys():
						#~ for entry in fullset[key]:
							#~ if channel in entry.keys():
								#~ if entry[channel] == value:
									#~ if variable in entry.keys() and "time" in entry.keys():
										#~ passed= True
										#~ t = entry["time"]
										#~ if isinstance(t, list):
											#~ t = t[0]
										#~ time.append(float(t))
										#~ ydata.append(float(entry[variable]))		
						
					#~ if passed and len(time) > 2:
						#~ accepted.append(fullset['name'])
						#~ allydata.append(ydata)
						#~ alltimes.append(time)
					#~ else:
						#~ nreject += 1
			
			#~ print "In total", key, "data with requirement", channel, "=", value,"and at least three data points was provided for", str(len(accepted)), "TDE candidates." 
			#~ print "Not provided for", str(nreject), "further entries. \n"
			
			#~ def curve(x, a, b, c):
				#~ return a * ((x-x0) **b) + c
			
			#~ fig = plt.figure()
			#~ npoints = len(accepted)
			#~ nrows = int(float(npoints)/2.) + (npoints % 2)
			#~ if npoints > 0:
				#~ for i in range(0, npoints):
					#~ t = alltimes[i]
					#~ y = allydata[i]
					
					#~ starttime = alltimes[i][0]
					
					#~ ax = plt.subplot(nrows,2, i+1)
					#~ plt.title(accepted[i])
					#~ plt.gca().invert_yaxis()
					#~ ax.set_yscale('log')
					
					#~ x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=starttime)
					#~ ax.xaxis.set_major_formatter(x_formatter)
					
					#~ ax.set_xlim([min(t), max(t)+1])
					#~ ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(11))
					#~ plt.scatter(t, y)
					#~ plt.xlabel("Time (MJD)")
					#~ plt.ylabel("Magnitude")
					
					#~ def curve(x, a, b, c):
						#~ return a * ((x-starttime)**b) +c
					
					#~ popt, pcov = curve_fit(curve, t, y)
					#~ print popt, pcov, popt[0]
					
					#~ fity=[]
					#~ trange = np.linspace(min(t), max(t)+1, 100)
					#~ for time in trange:
						#~ fity.append(curve(time, popt[0], popt[1], popt[2]))
					
					#~ plt.plot(trange, fity)
					#~ line = "y = %.2G" % (popt[0]) + "(t ^ %.2G)" %(popt[1]) + " + %.2G" %(popt[2])
					#~ print line
					
				
				#~ fig.set_size_inches(25, nrows*5)
				#~ fig.subplots_adjust(hspace=.5)
				#~ plt.savefig("graphs/"+variable+"_"+channel+"="+value+".pdf")
			
			#~ plt.close()
			
			
	
