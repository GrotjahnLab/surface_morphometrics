#! /usr/bin/env python
"""Step one of the pipeline for surface morphometrics
Convert a semantic segmentation MRC file to a series of membrane meshes for each segment.
Takes a label MRC file and outputs a point cloud xyz file, a mesh ply file, and a mesh vtp file.
Uses a config yml file to configure the conversion process.

Two usage options:
1. Convert a single segmentation file to a series of meshes:
  morphometrics make_meshes config.yml segmentation.mrc
2. Convert a series of segmentation files to a series of meshes:
  morphometrics make_meshes config.yml"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import sys
import subprocess
import glob

import click
import yaml

from .config_utils import load_config, require_keys

from . import mrc2xyz
from . import ply2vtp


def run_xyz_to_ply(xyz_file, ply_file, surface_config):
    """Run xyz2ply in a subprocess.

    Running in a subprocess avoids loading pymeshlab (Qt5) in the same
    process as vtk (Qt6).
    """
    sg = surface_config
    cmd = [
        sys.executable, "-m", "surface_morphometrics.xyz2ply", xyz_file, ply_file,
        "--pointweight", str(sg["point_weight"]),
        "--simplify", str(sg["simplify"]),
        "--num_faces", str(sg["simplify_max_triangles"]),
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


@click.command(name="make_meshes")
@click.argument("configfile", type=click.Path(exists=True))
@click.argument("segmentation", required=False, default=None)
def make_meshes_cli(configfile, segmentation):
    """Convert semantic segmentation MRC files into membrane meshes.

    CONFIGFILE: path to config.yml.
    SEGMENTATION: optional single segmentation .mrc to process; if omitted, all
    MRC files in seg_dir are processed.
    """
    # Check for a seg dir and a work dir
    config = load_config(configfile)

    # Warn about renamed config keys (see README). Older config files used
    # different names for these settings and will no longer work as-is.
    renamed_keys = []
    if "data_dir" in config:
        renamed_keys.append("  'data_dir' has been renamed to 'seg_dir' (the directory of segmentation MRC files).")
    if "max_triangles" in config.get("surface_generation", {}):
        renamed_keys.append("  'max_triangles' has been renamed to 'simplify_max_triangles' (only used when simplify: true).")
    if renamed_keys:
        print("=" * 70)
        print("WARNING: Your config.yml uses outdated key names that have changed:")
        for msg in renamed_keys:
            print(msg)
        print("Please update your config.yml to use the new key names.")
        print("See the README for details on the current configuration format.")
        print("=" * 70)

    # Validate after the rename warning above, so a user still on the old key names
    # sees the migration hint rather than only a bare "missing seg_dir".
    require_keys(config, ("seg_dir", "work_dir", "segmentation_values"), configfile)

    # See if a specific file was specified
    if segmentation is None:
        print("No input file specified - will run on all MRC files in the data directory")
        print("Pattern Matched: " + config["seg_dir"] + "*.mrc")
        segmentation_files = glob.glob(config["seg_dir"] + "*.mrc")
        segmentation_files = [os.path.basename(f) for f in segmentation_files]
        print(segmentation_files)
    else:
        print("Input file specified - will run on this file only")
        segmentation_files = [segmentation]

    # Load segmentation labels from config file
    seg_values = config["segmentation_values"]
    print("Configured segmentation labels:")
    print(seg_values)

    # Check that work dir exists for outputs
    if not os.path.isdir(config["work_dir"]):
        os.mkdir(config["work_dir"])

    # Run the conversion
    for file in segmentation_files:
        print(f"Processing segmentation {file}")
        basename = file[:-4]
        input = config["seg_dir"] + file
        for key, value in seg_values.items():
            # Convert the segmentation file to xyz files
            xyz_file = config["work_dir"] + basename + "_" + str(key) + ".xyz"
            print(f"Generating xyz file: {xyz_file}")
            ret_val = mrc2xyz.mrc_to_xyz(input, xyz_file, value, config["surface_generation"]["angstroms"])
            if ret_val != 0:
                print("Error converting segmentation file to xyz")
                continue
            # Generate the membrane mesh ply file from the xyz file
            ply_file = f"{config['work_dir']}{basename}_{key}.ply"
            print(f"Generating a ply mesh with Screened Poisson: {ply_file}")
            ret_val = run_xyz_to_ply(xyz_file, ply_file, config["surface_generation"])
            if ret_val != 0:
                print("Error converting xyz file to ply")
                continue
            # Convert the ply file to a vtp file
            vtp_file = config["work_dir"] + basename + "_" + str(key) + ".surface.vtp"
            print(f"Converting the ply file to a vtp file: {vtp_file}")
            ply2vtp.ply_to_vtp(ply_file, vtp_file)

    print("-------------------------------------------------------")
    print("Conversions complete. It is highly recommended to check the intermediate files for correctness before proceeding. All files are in the work directory.")
    print("Citation: Barad BA*, Medina M*, Fuentes D, Wiseman RL, Grotjahn DA. Quantifying organellar ultrastructure in cryo-electron tomography using a surface morphometrics pipeline. J Cell Biol 2023.")


if __name__ == "__main__":
    make_meshes_cli()
