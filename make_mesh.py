# Generate a mesh from a point cloud using pymeshlab and the screened poisson method
# I intentionally undercrop here with the understanding that I will crop again later
# to get the final mesh.
# Author: Benjamin Barad <benjamin.barad@gmail.com>
from sys import argv
import glob
import pymeshlab as pm



def reconstruct_surface(xyzfile, plyfile, pointweight=0.7, num_faces=150000, k_neighbors=70, deldist=4):

def xyz_to_ply(xyzfile, plyfile, pointweight=0.7, num_faces=150000, k_neighbors=70, deldist=4):
    """Convert an xyz file to a ply file using pymeshlab

    Arguments:
    xyzfile {str} -- Input xyz filename
    plyfile {str} -- Output ply filename
    pointweight {float} -- Screening weight (0 for max smoothness, 1 to 4 for beter fit to points). Default 0.7.
    num_faces {int} -- Maximal number of allowed faces after decimation. Default 150000, use more for finer sampling but with greater computational cost.
    k_neighbors {int} -- Number of neighbors for point cloud normal estimation - default 70
    deldist {int} -- Max distance to extrapolate. Default 4; distances are in the point cloud distance unit (default nm).
    """
    print(f"Processing {xyzfile}")
    ms = pm.MeshSet()
    ms.load_new_mesh(point_cloud)
    ms.compute_normals_for_point_sets(k=k_neighbors, smoothiter=1) # Predict smooth normals
    ms.surface_reconstruction_screened_poisson(depth=10, pointweight=pointweight, iters=10, scale=1.2) # Screened Poisson
    ms.distance_from_reference_mesh(measuremesh = 1, refmesh=0 , maxdist=pm.Percentage(20), signeddist=False) # Delete points that are too far from the reference mesh
    ms.conditional_vertex_selection(condselect = f'(q>{deldist})') # Select only the best quality vertices
    ms.conditional_face_selection(condselect = f'(q0>{deldist} || q1>{deldist} || q2>{deldist})') # Select only the best quality vertices
    ms.delete_selected_faces_and_vertices()
    ms.simplification_quadric_edge_collapse_decimation(targetfacenum=num_faces, qualitythr=0.6, preserveboundary=True, preservenormal=True, optimalplacement=True, planarquadric=True) # Simplify
    ms.save_current_mesh(plyfile)
    ms.clear()
    

if __name__ == "__main__":
    xyzfiles = sorted(glob.glob("*.xyz"))
    for point_cloud in xyzfiles:
        print(f"Processing {point_cloud}")
        plyfile = f"{point_cloud[:-4]}.ply"
        xyz_to_ply(point_cloud, plyfile)

