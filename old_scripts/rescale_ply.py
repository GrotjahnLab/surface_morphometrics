import os
from sys import argv
from glob import glob
import pymeshlab as pm
import mrcfile

mrcfiles = glob("*_labels.mrc")
deldist = 3
for mrc in mrcfiles:
    print(mrc)
    with mrcfile.mmap(mrc, 'r', permissive=True) as file:
        voxel_size = file.voxel_size.x
        print(voxel_size)
    ply_base = mrc[:-11]
    endings = ["IMM", "ER", "OMM"]
    for ending in endings:
        
        ptsfile = ply_base+"_"+ending+".pts"
        filename = ply_base + "_" + ending + ".ply"
        outfile = ply_base + "_" + ending + "_scaled.ply"
        ms = pm.MeshSet()
        if os.path.isfile(filename):
            print(filename)
            ms.load_new_mesh(ptsfile)
            ms.load_new_mesh(filename)
            ms.distance_from_reference_mesh(measuremesh = 1, refmesh=0 , maxdist=pm.Percentage(20), signeddist=False) # Delete points that are too far from the reference mesh
            ms.conditional_vertex_selection(condselect = f'(q>{deldist})') # Select only the best quality vertices
            ms.conditional_face_selection(condselect = f'(q0>{deldist} || q1>{deldist} || q2>{deldist})') # Select only the best quality vertices
            ms.delete_selected_faces_and_vertices()

            ms.transform_scale_normalize(axisx=voxel_size, uniformflag=True, scalecenter=0)
            ms.save_current_mesh(outfile)
            ms.clear()
        else:
            print("File not found: {}".format(filename))
