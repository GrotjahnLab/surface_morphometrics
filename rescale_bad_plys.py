from re import A
import pymeshlab as pm
from sys import argv
import glob
import os
import vtk


def fix_ply(ply, xyz, plyout, deldist = 1.5, num_faces = 150000, simplify = True):
    """Fix plys that have been scaled badly and not properly masked or simplified

    Arguments:
    ply {str} -- Input ply filename
    xyz {str} -- Input xyz filename
    plyout {str} -- Output ply filename
    num_faces {int} -- Maximal number of allowed faces after decimation. Default 150000, use more for finer sampling but with greater computational cost.
    deldist {int} -- Max distance to extrapolate. Default 4; distances are in the point cloud distance unit (default nm).
    smooth_iter {int} -- Number of smoothing iterations. Default 1.
    """
    print(f"Processing {ply} with {xyz}")
    ms = pm.MeshSet()
    ms.load_new_mesh(xyz)
    ms.load_new_mesh(ply)
    ms.compute_matrix_from_scaling_or_normalization(axisx=0.1, axisy=0.1, axisz=0.1, freeze=True)
    ms.distance_from_reference_mesh(measuremesh = 1, refmesh=0 , maxdist=pm.Percentage(20), signeddist=False) # Delete points that are too far from the reference mesh
    ms.conditional_vertex_selection(condselect = f'(q>{deldist})') # Select only the best quality vertices
    ms.conditional_face_selection(condselect = f'(q0>{deldist} || q1>{deldist} || q2>{deldist})') # Select only the best quality vertices
    ms.delete_selected_faces_and_vertices()
    if simplify:
        ms.simplification_quadric_edge_collapse_decimation(targetfacenum=num_faces, qualitythr=0.6, preserveboundary=True, preservenormal=True, optimalplacement=True, planarquadric=True) # Simplify
    ms.save_current_mesh(plyout)
    ms.clear()
    return 0

if __name__ == "__main__":
    folder = argv[1]
    os.chdir(folder)
    files = glob.glob("*_scaled.ply")
    for file in files:
        fileout = file.replace("_scaled.ply", "_rescaled.ply")
        xyzfile = "_".join(file.split("_")[0:2])+".xyz"
        fix_ply(file, xyzfile, fileout)
        vtpout = fileout.replace("_rescaled.ply", ".surface.vtp")
        plyfile = vtk.vtkPLYReader()
        plyfile.SetFileName(fileout)
        plyfile.Update()
        surf = plyfile.GetOutput() 
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(vtpout)
        writer.SetInputData(surf)
        if writer.Write() != 1:
            raise pexceptions.PySegInputError(
                expr='save_vtp', msg='Error writing the file {}.'.format(fname))



