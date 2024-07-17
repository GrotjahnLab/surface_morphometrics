#! /usr/bin/env python
"""Generate a cropped mesh from a point cloud using pymeshlab and the screened poisson method.

Usage: xyz2ply.py cloud.xyz mesh.ply"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

from sys import argv
import importlib.metadata
import glob
import pymeshlab as pm
import click

# Ensure compatibility with pymeshlab 2022.2.post3
PML_VER = importlib.metadata.version('pymeshlab')
if PML_VER == '2022.2.post3':
    pm.PercentageValue = pm.Percentage

@click.command()
@click.argument('xyzfile', type=str)
@click.argument('plyfile', type=str)
@click.option('--pointweight', type=float, default=1)
@click.option('--simplify', type=bool, default=False)
@click.option('--ultrafine', type=bool, default=True)
@click.option('--num_faces', type=int, default=300000)
@click.option('--k_neighbors', type=int, default=40)
@click.option('--deldist', type=int, default=1)
@click.option('--smooth_iter', type=int, default=1)
@click.option('--depth', type=int, default=9)
def xyz_to_ply_from_CLI(xyzfile, plyfile, pointweight, ultrafine, simplify, num_faces, k_neighbors, deldist, smooth_iter, depth):
    """Click wrapper for xyz_to_ply"""
    xyz_to_ply(xyzfile, plyfile, pointweight=pointweight, simplify=simplify, num_faces=num_faces, k_neighbors=k_neighbors, deldist=deldist, smooth_iter=smooth_iter, depth=depth)


def xyz_to_ply(xyzfile, plyfile, pointweight=0.1, simplify=False, num_faces=150000, k_neighbors=500, deldist=1.5, smooth_iter=1, depth=9, ultrafine=True, remesh_sampling=0.1):
    """Convert an xyz file to a ply file using pymeshlab

    Arguments:
    xyzfile {str} -- Input xyz filename
    plyfile {str} -- Output ply filename
    pointweight {float} -- Screening weight (0 for max smoothness, 1 to 4 for beter fit to points). Default 0.7.
    simplify {bool} -- If True, will simplify the mesh to a set number of faces. Default False. Leave False for ultrafine meshes.
    num_faces {int} -- Maximal number of allowed faces after decimation. Default 150000, use more for finer sampling but with greater computational cost.
    k_neighbors {int} -- Number of neighbors for point cloud normal estimation - default 400
    deldist {int} -- Max distance to extrapolate. Default 1.5; distances are in the point cloud distance unit (default nm).
    smooth_iter {int} -- Number of smoothing iterations. Default 1.
    ultrafine {bool} -- If True, will use a new routine to generate finer, higher quality meshes. Default True.
    simplify {bool} -- If True, will simplify the mesh to a set number of faces. Default True.
    remesh_sampling {int} -- Sampling rate in world units for the remeshing. Default 1.
    """
    print(f"Processing {xyzfile} into {plyfile}")
    ms = pm.MeshSet()
    ms.load_new_mesh(xyzfile)
    ms.compute_normal_for_point_clouds(k=k_neighbors, smoothiter=smooth_iter) # Predict smooth normals
    ms.generate_surface_reconstruction_screened_poisson(depth=depth, pointweight=pointweight, samplespernode=5.,iters=10, scale=1.2) # Screened Poisson
    if ultrafine:
        ms.meshing_surface_subdivision_loop(iterations=18)
        ms.generate_resampled_uniform_mesh(cellsize=pm.PercentageValue(0.05))
        ms.meshing_surface_subdivision_loop(iterations=6)
        ms.compute_scalar_by_distance_from_another_mesh_per_vertex(measuremesh=2, refmesh=0, maxdist = pm.PercentageValue(0.2), signeddist=False)
    else:
        ms.compute_scalar_by_distance_from_another_mesh_per_vertex(measuremesh = 1, refmesh=0 , maxdist=pm.PercentageValue(0.2), signeddist=False) # Delete points that are too far from the reference mesh
    ms.compute_selection_by_condition_per_vertex(condselect = f'(q>{deldist})') # Select only the best quality vertices
    ms.compute_selection_by_condition_per_face(condselect = f'(q0>{deldist} || q1>{deldist} || q2>{deldist})') # Select only the best quality vertices
    ms.meshing_remove_selected_vertices_and_faces()
    if simplify:
        ms.meshing_decimation_quadric_edge_collapse(targetfacenum=num_faces, qualitythr=0.6, preserveboundary=True, preservenormal=True, optimalplacement=True, planarquadric=True) # Simplify
    ms.save_current_mesh(plyfile)
    ms.clear()
    return 0
    

if __name__ == "__main__":
    xyzfiles = sorted(glob.glob("*.xyz"))
    for point_cloud in xyzfiles:
        print(f"Processing {point_cloud}")
        plyfile = f"{point_cloud[:-4]}.smooth.ply"
        xyz_to_ply(point_cloud, plyfile)

