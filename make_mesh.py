# Generate a mesh from a point cloud using pymeshlab and the screened poisson method
# I intentionally undercrop here with the understanding that I will crop again later
# to get the final mesh.
# Author: Benjamin Barad <benjamin.barad@gmail.com>
from sys import argv
import glob
import pymeshlab as pm

pts_files = sorted(glob.glob("*.pts"))
pointweight = 2 # 0 for max smoothness, 1 to 4 for beter fit to points
num_faces = 150000
# percent = .5 # percent of points to use
k_neighbors = 50 # number of neighbors for point cloud normal estimation
deldist = 1.5 # max distance to delete
ms = pm.MeshSet()

for point_cloud in pts_files:
    print(f"Processing {point_cloud}")
    # point_cloud = argv[1]
    ms.load_new_mesh(point_cloud)
    ms.compute_normals_for_point_sets(k=k_neighbors, smoothiter=2) # Predict smooth normals
    ms.surface_reconstruction_screened_poisson(depth=10, pointweight=pointweight, iters=10, scale=1.2) # Screened Poisson
    ms.distance_from_reference_mesh(measuremesh = 1, refmesh=0 , maxdist=pm.Percentage(20), signeddist=False) # Delete points that are too far from the reference mesh
    ms.conditional_vertex_selection(condselect = f'(q>{deldist})') # Select only the best quality vertices
    ms.conditional_face_selection(condselect = f'(q0>{deldist} || q1>{deldist} || q2>{deldist})') # Select only the best quality vertices

    ms.delete_selected_faces_and_vertices()
    ms.simplification_quadric_edge_collapse_decimation(targetfacenum=num_faces, qualitythr=0.6, preserveboundary=True, preservenormal=True, optimalplacement=True, planarquadric=True) # Simplify
    ms.save_current_mesh(f"{point_cloud[:-4]}.ply")
    ms.clear()
