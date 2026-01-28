#! /usr/bin/env python
"""Step one of the pipeline for surface morphometrics
Convert a semantic segmentation MRC file to a series of membrane meshes for each segment.
Takes a label MRC file and outputs a point cloud xyz file, a mesh ply file, and a mesh vtp file.
Uses a config yml file to configure the conversion process.

Two usage options:
1. Convert a single segmentation file to a series of meshes:
  segmentation_to_meshes.py config.yml segmentation.mrc
2. Convert a series of segmentation files to a series of meshes:
 segmentation_to_meshes.py config.yml"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import sys
import subprocess
import glob
from sys import argv

import yaml

import mrc2xyz
import ply2vtp

# Check for a config file
if len(argv) < 2:
    print("Usage: segmentation_to_meshes.py config.yml [segmentation.mrc]")
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

# See if a specific file was specified
if len(argv) == 2:
    print("No input file specified - will run on all MRC files in the data directory")
    print("Pattern Matched: "+config["data_dir"]+"*.mrc")
    segmentation_files = glob.glob(config["data_dir"]+"*.mrc")
    segmentation_files = [os.path.basename(f) for f in segmentation_files]
    print(segmentation_files)
else:
    print("Input file specified - will run on this file only")
    segmentation_files = [argv[2]]

# Load segmentation labels from config file
seg_values = config["segmentation_values"]
print("Configured segmentation labels:")
print(seg_values)

# Check that work dir exists for outputs
if not os.path.isdir(config["work_dir"]):
    os.mkdir(config["work_dir"])

def run_xyz_to_ply(xyz_file, ply_file, surface_config):
    """Run xyz2ply in a subprocess.

    Running in a subprocess avoids loading pymeshlab (Qt5) in the same
    process as vtk (Qt6).
    """
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "xyz2ply.py")
    sg = surface_config
    cmd = [
        sys.executable, script, xyz_file, ply_file,
        "--pointweight", str(sg["point_weight"]),
        "--simplify", str(sg["simplify"]),
        "--num_faces", str(sg["max_triangles"]),
        "--k_neighbors", str(sg["neighbor_count"]),
        "--deldist", str(sg["extrapolation_distance"]),
        "--smooth_iter", str(sg["smoothing_iterations"]),
        "--depth", str(sg["octree_depth"]),
        "--isotropic_remesh", str(sg.get("isotropic_remesh", True)),
        "--target_area", str(sg.get("target_area", 1.0)),
    ]
    result = subprocess.run(cmd)
    if os.path.exists(ply_file) and os.path.getsize(ply_file) > 0:
        return 0
    return 1


# Run the conversion
for file in segmentation_files:
    print(f"Processing segmentation {file}")
    basename = file[:-4]
    input = config["data_dir"]+file
    for key,value in seg_values.items():
        # Convert the segmentation file to xyz files
        xyz_file = config["work_dir"]+basename+"_"+str(key)+".xyz"
        print(f"Generating xyz file: {xyz_file}")
        ret_val = mrc2xyz.mrc_to_xyz(input, xyz_file, value, config["surface_generation"]["angstroms"]) # Convert the segmentation file to xyz files
        if ret_val != 0:
            print("Error converting segmentation file to xyz")
            continue
        #Generate the membrane mesh ply file from the xyz file
        ply_file =   f"{config['work_dir']}{basename}_{key}.ply"
        print(f"Generating a ply mesh with Screened Poisson: {ply_file}")
        ret_val = run_xyz_to_ply(xyz_file, ply_file, config["surface_generation"])
        if ret_val != 0:
            print("Error converting xyz file to ply")
            continue
        # Convert the ply file to a vtp file
        vtp_file = config["work_dir"]+basename+"_"+str(key)+".surface.vtp"
        print(f"Converting the ply file to a vtp file: {vtp_file}")
        ply2vtp.ply_to_vtp(ply_file, vtp_file)


print("-------------------------------------------------------")
print("Conversions complete. It is highly recommended to check the intermediate files for correctness before proceeding. All files are in the work directory.")
print("If you are happy with the results, you can move on to `run_pycurv.py`.")
print("Citation: Barad BA*, Medina M*, Fuentes D, Wiseman RL, Grotjahn DA. Quantifying organellar ultrastructure in cryo-electron tomography using a surface morphometrics pipeline. J Cell Biol 2023.")