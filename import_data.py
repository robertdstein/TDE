import os

url ='https://github.com/astrocatalogs/tde-1980-2025/archive/master.zip'

def run(savepath):
    print "Downloading latest TDE Catalogue from", url
    # f = urllib2.urlopen(url)
    # data = f.read()

    cmd = "wget {0} -o {1}".format(url, savepath)
    print cmd
    os.system(cmd)

    # with open(savepath, "wb") as code:
    #     print "Saving to", savepath
    #     code.write(data)

    print "Catalogue is now updated to latest available version! \n"
