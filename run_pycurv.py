"""Step two of the pipeline for surface morphometrics
Take a VTP style mesh generatead through the screen and perform curvature analysis using AVV.

Citation: Barad BA, Medina M et al. A surface morphometrics toolkit to quantify organellar membrane ultrastructure using cryo-electron tomography. Biorxiv 2022.
Pycurv Citation: Salfer M et al. Reliable estimation of membrane curvature for cryo-electron tomography. PLOS Computational Biology 2020.

Two usage options:
1. Run pycurv on a single mesh (Recommended!):
  run_pycurv.py config.yml mesh.vtp
2. Run pycurv on all meshes in the working directory:
 segmentation_to_meshes.py config.yml
 
Because pycurv is quite resource intensive, it is recommended to use option 1 with a cluster submission script in parallel, rather than sequentially running all vtp files.
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

from sys import argv
import os
import glob

import yaml

import curvature

# Check for a config file
if len(argv) < 2:
    print("Usage: run_pycurv.py config.yml [filename.surface.vtp]")
    exit()

# Check for a data dir and a work dir
with open(argv[1]) as file:
    config = yaml.safe_load(file)
    if not config["work_dir"]:
        if not config["data_dir"]:
            print("No working directory is specified in the config file. Please specify a working directory or a data directory.")
            exit()
        else:
            print("No working directory is specified in the config file. The data directory will be used for input and output.")
            config["work_dir"] = config["data_dir"]

# Figure out what files will be run
if len(argv) == 2:
    print("No input file specified - will run on all VTP files in the data directory")
    print("This may take a very long time - pycurv can take over an hour to run on a single mesh")
    print("It is recommended to run in parallel with a cluster submission script for individual files")
    print("Recommended usage: run_pycurv.py config.yml <meshname.surface.vtp>")
    answer = input("Continue? [y/n]")
    if answer != "y":
        exit()
    mesh_files = glob.glob(config["work_dir"]+"*.surface.vtp")
    mesh_files = [os.path.basename(f) for f in mesh_files]
else:
    print("Input file specified - will run on this file only")
    mesh_files = [argv[2]]

# Check that work dir exists for outputs
if not os.path.isdir(config["work_dir"]):
    os.mkdir(config["work_dir"])

for surface in mesh_files:
    print("Processing "+surface)
    curvature.run_pycurv(surface, config["work_dir"],
                        scale=1.0,
                        radius_hit=config["curvature_measurements"]["radius_hit"],
                        min_component=config["curvature_measurements"]["min_component"],
                        exclude_borders=config["curvature_measurements"]["exclude_borders"],
                        cores=config["cores"])

print("-------------------------------------------------------")
print("Pycurv complete. It is highly recommended to check the AVV vtp file with paraview to confirm good results.")
print("If you are happy with the results, you can move on to `distances_and_orientations.py`.")
print("Pycurv Citation: Salfer M, Collado JF, Baumeister W, Fernández-Busnadiego R, Martínez-Sánchez A. Reliable estimation of membrane curvature for cryo-electron tomography. PLOS Comp Biol 2020.")
print("Pipeline Citation: Barad BA*, Medina M*, Fuentes D, Wiseman RL, Grotjahn DA. A surface morphometrics toolkit to quantify organellar membrane ultrastructure using cryo-electron tomography. Biorxiv 2022.")


    





