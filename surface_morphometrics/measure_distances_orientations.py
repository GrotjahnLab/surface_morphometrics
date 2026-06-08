#! /usr/bin/env python
"""Step three of the pipeline

Take a series of triange graphs in gt format from pycurv and measure
distances within and between surfaces as well as orientations relative to
the growth plane and to each other.

Two usage options:
1. Measure distances and orientations on a single **tomogram** worth of graphs - yes, supply a segmentation file, not the meshes.
  morphometrics distances_orientations config.yml segmentation.mrc
2. Measure distances and orientations for the graphs for every tomogram in a folder:
  morphometrics distances_orientations config.yml

Since these distance measurements can take a minute or two per measurement per surface,
it may be worth doing this in parallel using option 1."""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import glob
import sys

import click
import yaml

from . import intradistance_verticality
from . import interdistance_orientation


@click.command(name="distances_orientations")
@click.argument("configfile", type=click.Path(exists=True))
@click.argument("segmentation", required=False, default=None)
@click.option("-f", "--force", is_flag=True, default=False,
              help="Skip the interactive confirmation prompt when running on all files.")
def distances_orientations_cli(configfile, segmentation, force):
    """Measure intra- and inter-surface distances and orientations.

    CONFIGFILE: path to config.yml.
    SEGMENTATION: optional single segmentation .mrc whose graphs to process; if
    omitted, all segmentations in seg_dir are processed.
    """
    with open(configfile) as file:
        config = yaml.safe_load(file)
    if not config["seg_dir"]:
        print("seg_dir not specified in config.yml")
        sys.exit(1)
    elif not config["seg_dir"].endswith("/"):
        config["seg_dir"] += "/"
    if not config["work_dir"]:
        print("work_dir not specified in config.yml - seg_dir will be used for output")
        config["work_dir"] = config["seg_dir"]
    elif not config["work_dir"].endswith("/"):
        config["work_dir"] += "/"

    # See if a specific file was specified
    if segmentation is None:
        print("No input file specified - will run on meshes for all segmentation files in the data directory")
        print("This may take a bit of time - each measurement will take a minute or two per surface, and they can add up.")
        print("You may prefer to run in parallel with a cluster submission script for individual files")
        print("The example config and tutorial data takes about 12-15 minutes on a laptop.")
        print("Recommended usage: morphometrics distances_orientations config.yml <segmentation.mrc>")
        if not force:
            answer = input("Continue? [y/n]")
            if answer != "y":
                sys.exit(1)
        print("Pattern Matched: " + config["seg_dir"] + "*.mrc")
        segmentation_files = glob.glob(config["seg_dir"] + "*.mrc")
        segmentation_files = [os.path.basename(f) for f in segmentation_files]
        print(segmentation_files)
    else:
        print("Input file specified - will run on meshes associated with this segmentation only")
        segmentation_files = [segmentation]

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
            print(surface + f": {comparison}")
        if dist_settings["relative_orientation"]:
            print("Will also measure relative orientation for those surfaces")
    print("---------")

    radius_hit = config["curvature_measurements"]["radius_hit"]
    for file in segmentation_files:
        file_base = file[:-4]
        print(file_base)
        # Intra-surface distances
        if dist_settings["intra"]:
            for label in dist_settings["intra"]:
                graphname = config["work_dir"] + file_base + "_" + label + ".AVV_rh" + str(radius_hit) + ".gt"
                if not os.path.isfile(graphname):
                    print("No file found for " + graphname)
                    print("Skipping this label for this tomogram")
                    continue
                print("Intra-surface distances for " + graphname)
                surfacename = graphname[:-3] + ".vtp"
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
                graphname1 = config["work_dir"] + file_base + "_" + label1 + ".AVV_rh" + str(radius_hit) + ".gt"
                if not os.path.isfile(graphname1):
                    print("No file found for " + graphname1)
                    print("Skipping all intersurface measurements for label " + label1)
                    continue
                for label2 in comparison:
                    graphname2 = config["work_dir"] + file_base + "_" + label2 + ".AVV_rh" + str(radius_hit) + ".gt"
                    if not os.path.isfile(graphname2):
                        print("No file found for " + graphname2)
                        print("Skipping comparison with " + label2)
                        continue
                    print("Inter-surface distances for " + label1 + " and " + label2)
                    interdistance_orientation.surface_to_surface(graphname1, label1,
                                                                 graphname2, label2,
                                                                 orientation=dist_settings["relative_orientation"],
                                                                 save_neighbor_index=True,
                                                                 exportcsv=True)

    print("---------")
    print("Done with distance and orientation measurements. Check out the output vtp and csv files and then start getting to work doing stats!")
    print("Citation: Barad BA*, Medina M*, Fuentes D, Wiseman RL, Grotjahn DA. Quantifying organellar ultrastructure in cryo-electron tomography using a surface morphometrics pipeline. J Cell Biol 2023.")


if __name__ == "__main__":
    distances_orientations_cli()
