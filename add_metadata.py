import json
import zipfile
import os


def run(zip_file_name="/afs/ifh.de/user/s/steinrob/Desktop/python/TDE/tde_cat"
                      ".zip"):
    root = "/" .join(zip_file_name.split("/")[:-1]) + "/"

    metadata_path = root + "metadata.json"
    print "Loading additional metadata in", metadata_path

    zip_path = zip_file_name
    temp_path = root + "temp.zip"
    with zipfile.ZipFile(zip_path, 'r') as zf:
        with zipfile.ZipFile(temp_path, 'w') as tf:
            filelist = zf.namelist()[2:]
            with open(metadata_path) as metadata:
                md = json.load(metadata)
                for filename in filelist:
                    with zf.open(filename, "r") as f:
                        data = f.read()
                        d = json.loads(data)
                        name = d.keys()[0]
                        if name in md.keys():
                            print "Found additional metadata for", name, ": ",
                            new_data = md[name]
                            for key in md[name].keys():
                                if key not in d[name].keys():
                                    print key,
                                    d[name][key] = new_data[key]
                            print ""

                        g = json.dumps(d)
                        tf.writestr(filename, g)
    os.remove(zip_path)
    os.rename(temp_path, zip_path)
    print ""
