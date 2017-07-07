import urllib2
 
url ='https://github.com/astrocatalogs/tde-1980-2025/archive/master.zip'

def run(savepath):
	print "Downloading latest TDE Catalogue from", url
	f = urllib2.urlopen(url)
	data = f.read()
	
	with open(savepath, "wb") as code:
		code.write(data)
		
	print "Catalogue is now updated to latest available version!"
	

