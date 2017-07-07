import json
import argparse
import zipfile
import datetime
import csv
import matplotlib.pyplot as plt
import import_data as i

parser = argparse.ArgumentParser(description='Toggle online updating on or off')
parser.add_argument("-u", "--update", action="store_true")

cfg = parser.parse_args()

savepath = "tde_cat.zip"
varpath = "variables.csv"

if cfg.update:
	i.run(savepath)

zf = zipfile.ZipFile(savepath, 'r')
print "Catalogue version created:", datetime.datetime(*zf.infolist()[0].date_time), "\n"

filelist = zf.namelist()[2:]

with open(varpath, 'rb') as csvfile:
	reader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in reader:
		
		key, variable, channel, value = row
		
		print "Finding", variable, "(", key, ") with requirement", channel, "=", value

		accepted =[]
		nreject = 0
		
		alltimes=[]
		allydata=[]
		
		for name in filelist:
			with zf.open(name) as f:  
				data = f.read()
				d = json.loads(data)
				fullset =  d.values()[0]
				
				time=[]
				ydata=[]
				
				passed = False
				
				if key in fullset.keys():
					for entry in fullset[key]:
						if channel in entry.keys():
							if entry[channel] == value:
								if variable in entry.keys() and "time" in entry.keys():
									passed= True
									time.append(entry["time"])
									ydata.append(entry[variable])		
					
				if passed:
					accepted.append(fullset['name'])
					allydata.append(ydata)
					alltimes.append(time)
				else:
					nreject += 1
		
		print "In total", key, "data with requirement", channel, "=", value,"was provided for", str(len(accepted)), "TDE candidates." 
		print "Not provided for", str(nreject), "further entries. \n"
		
		fig = plt.figure()
		npoints = len(accepted)
		nrows = int(float(npoints)/2.) + (npoints % 2)
		for i in range(0, npoints):
			plt.subplot(nrows,2, i+1)
			plt.title(accepted[i])
			plt.scatter(alltimes[i], allydata[i])
			plt.gca().invert_yaxis()
			plt.xlabel("time")
			plt.ylabel("Magnitude")
		
		print npoints, nrows
		fig.set_size_inches(25, 15)
		plt.savefig("graphs/"+variable+"_"+channel+"="+value+".pdf")
		
		
