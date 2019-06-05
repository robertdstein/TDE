import simplejson
import zipfile
import os
import re


def run(zip_file_name="/afs/ifh.de/user/s/steinrob/scratch/TDE_output/tde_cat"
        ".zip", source_dir = "/afs/ifh.de/user/s/steinrob/Desltop/python/TDE"):
    """Function to modify the ZIP file downloaded from the Open TDE
    catalogue. It incorporates additional metadata stored in the same
    directory, under the name metadata.json

    Creates a new ZIP file under the name temp.zip, which is filled by
    updated files. The original zip catalogue is then removed, and replaced
    with temp.zip

    :param zip_file_name: Path of TDE Zip file
    """

    # Finds root of catalogue, giving metadata path
    root = "/" .join(zip_file_name.split("/")[:-1]) + "/"

    metadata_path = source_dir + "metadata.json"
    print "Loading additional metadata in", metadata_path

    # Opens both old and new zip files concurrently
    zip_path = zip_file_name
    temp_path = root + "temp.zip"
    with zipfile.ZipFile(zip_path, 'r') as zf:
        with zipfile.ZipFile(temp_path, 'w') as tf:

            ignore = [
                "tde-1980-2025-master/",
                "tde-1980-2025-master/.gitignore",
                "tde-1980-2025-master/---.json"
            ]

            file_list = [x for x in zf.namelist() if x not in ignore]


            with open(metadata_path) as metadata:
                md = simplejson.load(metadata)

                # Loops over each file in TDE zip

                for filename in file_list:
                    print filename
                    with zf.open(filename, "r") as f:
                        data = f.read()
                        d = simplejson.loads(data)
                        name = d.keys()[0]

                        # Checks for additional metadata

                        if name in md.keys():
                            print "Found additional metadata for", name, ": ",
                            new_data = md[name]

                            # If new sources are present, renumbers the sources
                            # in the metadata file (Which are named according
                            # to NEW1, NEW2 etc.), with the appropriate
                            # integers, following the existing sources.

                            if "sources" in new_data.keys():

                                if "sources" in d[name].keys():
                                    source_no = len(d[name]["sources"]) + 1
                                else:
                                    source_no = 1
                                    d[name]["sources"] = []

                                for new_source in new_data["sources"]:
                                    if "NEW" in new_source["alias"]:
                                        new_id = str(source_no)
                                        old_json = simplejson.dumps(new_data)
                                        new_json = re.sub(
                                            new_source["alias"], new_id,
                                            old_json)
                                        new_data = simplejson.loads(new_json)
                                        source_no += 1

                            for key in new_data.keys():
                                if key not in d[name].keys():
                                    print key,
                                    d[name][key] = new_data[key]

                                elif key in ["sources", "alias"]:
                                    added = [x for x in new_data[key]]
                                    d[name][key].extend(added)

                                else:
                                    print "OVERWRITING:", key,
                                    d[name][key] = new_data[key]

                            print ""

                        # Sets the spacing structure of the json file to be
                        # tabs, in order to match the pre-existing format
                        # used by the Open TDE Catalog

                        g = simplejson.dumps(d, indent=4)
                        tabbed_g = re.sub(
                            '\n +', lambda match: '\n' + '\t' * (
                                len(match.group().strip('\n')) / 3),
                            g)

                        tf.writestr(filename, tabbed_g)

    os.remove(zip_path)
    os.rename(temp_path, zip_path)
    print ""
