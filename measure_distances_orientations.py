#! /usr/bin/env python
"""Step three of the pipeline

Take a series of triange graphs in gt format from pycurv and measure 
distances within and between surfaces as well as orientations relative to
the growth plane and to each other.

Two usage options:
1. Measure distances and orientations on a single **tomogram** worth of graphs to a series - yes, supply a segmentation file, not the meshes.
python measure_distances.py config.yml segmentation.mrc
2. Measure distances and orientations for the graphs for every tomogram in a folder:
python measure_distances.py config.yml

Since these distance measurements can take a minute or two per measurement per surface, 
it may be worth doing this in parallel using option 1."""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

from sys import argv
import os
import glob

import yaml

import intradistance_verticality
import interdistance_orientation

# Check for a config file
if len(argv) < 2:
    print("Usage: python measure_distances_orientations.py config.yml [segmentation.mrc]")
    exit()

# Check for a data dir and a work dir
with open(argv[1]) as file:
    config = yaml.safe_load(file)
    if not config["data_dir"]:
        print("data_dir not specified in config.yml")
        exit()
    elif not config["data_dir"].endswith("/"):
        config["data_dir"] += "/"
    if not config["work_dir"]:
        print("work_dir not specified in config.yml - data_dir will be used for output")
        config["work_dir"] = config["data_dir"]
    elif not config["work_dir"].endswith("/"):
        config["work_dir"] += "/"
# See if a specific file was specified
if len(argv) == 2:
    print("No input file specified - will run on meshes for all segmentation files in the data directory")
    print("This may take a bit of time - each measurement will take a minute or two per surface, and they can add up.")
    print("You may prefer to run in parallel with a cluster submission script for individual files")
    print("The example config and tutorial data takes about 12-15 minutes on a laptop.")
    print("Recommended usage: measure_distances_orientations.py config.yml <segmentation.mrc>")
    answer = input("Continue? [y/n]")
    if answer != "y":
        exit()
    print("Pattern Matched: "+config["data_dir"]+"*.mrc")
    segmentation_files = glob.glob(config["data_dir"]+"*.mrc")
    segmentation_files = [os.path.basename(f) for f in segmentation_files]
    print(segmentation_files)

else:
    print("Input file specified - will run on meshes associated with this segmentation only")
    segmentation_files = [argv[2]]

dist_settings = config["distance_and_orientation_measurements"]
print("Distance and orientation settings:")
print("Will measure intra-surface distances for:")
if dist_settings["intra"]:
    print(dist_settings["intra"])
    if dist_settings["verticality"]:
        print("Will also measure verticality for those surfaces")
if dist_settings["inter"]:
    print("Will make (bi-directional) inter-surface distance measurements for:")
    for surface, comparison in dist_settings["inter"].items():
        print(surface+f": {comparison}")
    if dist_settings["relative_orientation"]:
        print("Will also measure relative orientation for those surfaces")
print("---------")

for file in segmentation_files:
    file_base = file[:-4]
    print(file_base)
    # Intra-surface distances
    if dist_settings["intra"]:
        for label in dist_settings["intra"]:
            graphname = config["work_dir"]+file_base+"_"+label+".AVV_rh"+str(config["curvature_measurements"]["radius_hit"])+".gt"
            if not os.path.isfile(graphname):
                print("No file found for "+graphname)
                print("Skipping this label for this tomogram")
                continue
            print("Intra-surface distances for "+graphname)
            surfacename = graphname[:-3]+".vtp"
            if dist_settings["verticality"]:
                intradistance_verticality.surface_verticality(graphname)
            intradistance_verticality.surface_self_distances(graphname, surfacename,
                                                            dist_min=dist_settings["mindist"],
                                                            dist_max=dist_settings["maxdist"],
                                                            tolerance=dist_settings["tolerance"],
                                                            exportcsv=True)
    # Inter-surface distances
    if dist_settings["inter"]:
        for label1, comparison in dist_settings["inter"].items():
            graphname1 = config["work_dir"]+file_base+"_"+label1+".AVV_rh"+str(config["curvature_measurements"]["radius_hit"])+".gt"
            if not os.path.isfile(graphname1):
                print("No file found for "+graphname1)
                print("Skipping all intersurface measurements for label "+label1)
                continue
            for label2 in comparison:
                graphname2 = config["work_dir"]+file_base+"_"+label2+".AVV_rh"+str(config["curvature_measurements"]["radius_hit"])+".gt"
                if not os.path.isfile(graphname2):
                    print("No file found for "+graphname2)
                    print("Skipping comparison with "+label2)
                    continue
                print("Inter-surface distances for "+label1 + " and " + label2)
                interdistance_orientation.surface_to_surface(graphname1, label1,
                                                                        graphname2, label2,
                                                                        orientation=dist_settings["relative_orientation"],
                                                                        save_neighbor_index=True,
                                                                        exportcsv=True)





