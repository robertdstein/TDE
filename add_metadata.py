import json
import zipfile
import os
import re


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

                            print type(new_data)

                            if "sources" in new_data.keys():

                                if "sources" in d[name].keys():
                                    source_no = len(d[name]["sources"]) + 1
                                else:
                                    source_no = 1
                                    d[name]["sources"] = []

                                for new_source in new_data["sources"]:
                                    if "NEW" in new_source["alias"]:
                                        print new_source
                                        print new_data["sources"]
                                        new_ID = str(source_no)
                                        old_json = json.dumps(new_data)
                                        new_json = re.sub(
                                            new_source["alias"], new_ID,
                                            old_json)
                                        new_data = json.loads(new_json)
                                        source_no += 1
                                        print new_data["sources"]
                                    else:
                                        print "FALSE!"
                                raw_input("prompt")
                            else:
                                print "No new source provided for metadata"

                            for key in md[name].keys():
                                if key not in d[name].keys():
                                    print key,
                                    d[name][key] = new_data[key]

                                elif key == "sources":
                                    added = [x for x in new_data[key]]
                                    print d[name][key]
                                    d[name][key].extend(added)
                                    print d[name][key]
                                    print added

                            print ""

                        g = json.dumps(d, indent=4)
                        tabbed_g = re.sub(
                            '\n +', lambda match: '\n' + '\t' * (
                                len(match.group().strip('\n')) / 3),
                            g)
                        tf.writestr(filename, tabbed_g)

    os.remove(zip_path)
    os.rename(temp_path, zip_path)
    print ""
