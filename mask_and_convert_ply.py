import gc
import glob
import os

import click
import numpy as np
from pycurv import pycurv_io as io
from scipy.ndimage.morphology import distance_transform_edt
import mrcfile
import vtk

TOLERANCE_DISTANCE_NM=2 # Heuristic to get complete surfaces but avoid spurious junk
# Closed in for TgGSK because it was spurious junk!

@click.command()
@click.argument('input_surfaces', type=str, nargs=-1)
@click.argument('input_mrc', type=str, nargs=1)
# @click.argument('output_base', type=str, nargs=1)
@click.option('-l', '--labels', type=list, default=[1])
@click.option('-d', '--mask_tolerance_distance', type=int, default=3)
def convert_from_CLI(input_surfaces, input_mrc, mask_tolerance_distance, labels):
    """Click wrapper for convert and mask script"""
    assert(len(labels)==len(input_surfaces))
    convert(input_surfaces, input_mrc, mask_tolerance_distance, labels)

def convert(input_surfaces, input_mrc, mask_tolerance_distance, labels):
    print("open")
    with mrcfile.mmap(input_mrc, mode="r+", permissive=True) as mrc:
        header = mrc.header
        maximum = (header.nx, header.ny, header.nz)
        minimum = (0,0,0)
        # print(maximum)
        # print(minimum)
        voxel_size = mrc.voxel_size.x/10 # nm
        # maximum = [int(i) for i in maximum]
        print(maximum)
        # Since I am keeping everything to voxel spacing until the scale step in pycurv, I will just reduce mask_tolerance_distance according to voxel size
        mask_tolerance_distance = mask_tolerance_distance/voxel_size
        print(voxel_size, mask_tolerance_distance)
        # print([round(i) for i in (maximum/voxel_size)])
        for index,input_surface in enumerate(input_surfaces):
            label=labels[index]
            print(label, input_surface)

            output_base = input_surface[:-4]
            plyfile = vtk.vtkPLYReader()
            plyfile.SetFileName(input_surface)
            plyfile.Update()
            surf = plyfile.GetOutput() 
            writer = vtk.vtkXMLPolyDataWriter()
            fname = output_base + "_unmasked.vtp"
            writer.SetFileName(fname)
            writer.SetInputData(surf)
            if writer.Write() != 1:
                raise pexceptions.PySegInputError(
                    expr='save_vtp', msg='Error writing the file {}.'.format(fname))
            print("Making the mask")
            mask = (mrc.data == label).astype(int)
            # print(mask)
            dist_from_mask = distance_transform_edt(mask == 0)
            # print(dist_from_mask)
            for i in range(surf.GetNumberOfCells()):
                # Check if all points which made up the polygon are in the mask
                points_cell = surf.GetCell(i).GetPoints()
                # print(points_cell.GetPoint(0))
                count = 0
                for j in range(0, points_cell.GetNumberOfPoints()):
                    x,y,z = [int(round(i)) for i in points_cell.GetPoint(j)]
                    # print(x,y,z)
                    # print(dist_from_mask[z, y, x])

                    # x = int(round(x/voxel_size))
                    # y = int(round(y/voxel_size))
                    # z = int(round(z/voxel_size))
                    # print(x,y,z)
                    if z>=maximum[2] or z<=minimum[2] or x>=maximum[0] or x<=minimum[0] or y>=maximum[1] or y<=minimum[1]:
                        count += 1
                        continue
                    if (dist_from_mask[z, y, x] >mask_tolerance_distance):
                        count += 1
                    # else: print(x,y,z)
                # Mark cells that are not completely in the mask for deletion
                if count > 0:
                    surf.DeleteCell(i)

            # Delete
            surf.RemoveDeletedCells()
            
            writer = vtk.vtkXMLPolyDataWriter()
            fname = output_base + ".surface.vtp"
            writer.SetFileName(fname)
            writer.SetInputData(surf)
            if writer.Write() != 1:
                raise pexceptions.PySegInputError(
                    expr='save_vtp', msg='Error writing the file {}.'.format(fname))
            del surf
            del plyfile
            del writer
            del mask
            del dist_from_mask
            gc.collect()






if __name__=="__main__":
    print("test")
    # convert_from_CLI()
    # mrcfiles = glob.glob("*_labels.mrc")
    mrcfiles = ["TT9_labels.mrc"]
    structures = ["OMM", "IMM", "ER"]
    # structures=["IMM"]
    print(mrcfiles)
    for mrc in mrcfiles:
        print(mrc)
        basename = mrc[:-11]
        print(basename)
        surface_names = []
        labels = []
        for index,structure in enumerate(structures):
            label = index+1
            # print(label)
            filename = basename+"_"+structure+".ply" 
            if os.path.isfile(filename):
                labels.append(label)
                surface_names.append(filename)
            else:
                print(f"{filename} does not exist - skipping")
        convert(surface_names, mrc, mask_tolerance_distance=TOLERANCE_DISTANCE_NM, labels=labels) # 4nm tolerance to close holes




