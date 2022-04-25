#! /usr/bin/env python
"""Generate a cropped mesh from a point cloud using pymeshlab and the screened poisson method.

Usage: xyz2ply.py cloud.xyz mesh.ply"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

from sys import argv
import glob
import pymeshlab as pm
import click

@click.command()
@click.argument('xyzfile', type=str)
@click.argument('plyfile', type=str)
@click.option('--pointweight', type=float, default=1)
@click.option('--simplify', type=bool, default=True)
@click.option('--num_faces', type=int, default=150000)
@click.option('--k_neighbors', type=int, default=70)
@click.option('--deldist', type=int, default=2)
@click.option('--smooth_iter', type=int, default=1)
@click.option('--depth', type=int, default=9)
def xyz_to_ply_from_CLI(xyzfile, plyfile, pointweight, simplify, num_faces, k_neighbors, deldist, smooth_iter, depth):
    """Click wrapper for xyz_to_ply"""
    xyz_to_ply(xyzfile, plyfile, pointweight=pointweight, simplify=simplify, num_faces=num_faces, k_neighbors=k_neighbors, deldist=deldist, smooth_iter=smooth_iter, depth=depth)


def xyz_to_ply(xyzfile, plyfile, pointweight=1, simplify=True, num_faces=150000, k_neighbors=70, deldist=2, smooth_iter=1, depth=9):
    """Convert an xyz file to a ply file using pymeshlab

    Arguments:
    xyzfile {str} -- Input xyz filename
    plyfile {str} -- Output ply filename
    pointweight {float} -- Screening weight (0 for max smoothness, 1 to 4 for beter fit to points). Default 0.7.
    num_faces {int} -- Maximal number of allowed faces after decimation. Default 150000, use more for finer sampling but with greater computational cost.
    k_neighbors {int} -- Number of neighbors for point cloud normal estimation - default 70
    deldist {int} -- Max distance to extrapolate. Default 4; distances are in the point cloud distance unit (default nm).
    smooth_iter {int} -- Number of smoothing iterations. Default 1.
    """
    print(f"Processing {xyzfile} into {plyfile}")
    ms = pm.MeshSet()
    ms.load_new_mesh(xyzfile)
    ms.compute_normals_for_point_sets(k=k_neighbors, smoothiter=smooth_iter) # Predict smooth normals
    ms.surface_reconstruction_screened_poisson(depth=depth, pointweight=pointweight, iters=10, scale=1.2) # Screened Poisson
    ms.distance_from_reference_mesh(measuremesh = 1, refmesh=0 , maxdist=pm.Percentage(20), signeddist=False) # Delete points that are too far from the reference mesh
    ms.conditional_vertex_selection(condselect = f'(q>{deldist})') # Select only the best quality vertices
    ms.conditional_face_selection(condselect = f'(q0>{deldist} || q1>{deldist} || q2>{deldist})') # Select only the best quality vertices
    ms.delete_selected_faces_and_vertices()
    if simplify:
        ms.simplification_quadric_edge_collapse_decimation(targetfacenum=num_faces, qualitythr=0.6, preserveboundary=True, preservenormal=True, optimalplacement=True, planarquadric=True) # Simplify
    ms.save_current_mesh(plyfile)
    ms.clear()
    return 0
    

if __name__ == "__main__":
    xyzfiles = sorted(glob.glob("*.xyz"))
    for point_cloud in xyzfiles:
        print(f"Processing {point_cloud}")
        plyfile = f"{point_cloud[:-4]}.ply"
        xyz_to_ply(point_cloud, plyfile)

