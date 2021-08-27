import click
import numpy as np
from pycurv import pycurv_io as io
from scipy.ndimage.morphology import distance_transform_edt
import mrcfile
import vtk


@click.command()
@click.argument('input_surface', type=str)
@click.argument('input_mrc', type=str)
@click.argument('output_base', type=str)
@click.option('-l', '--label', type=int, default=1)
@click.option('-d', '--mask_tolerance_distance', type=int, default=3)
def convert(input_surface, input_mrc, output_base, mask_tolerance_distance, label):
    mrc = mrcfile.mmap(input_mrc, mode="r+", permissive=True)
    header = mrc.header
    maximum = (header.nx, header.ny, header.nz)
    minimum = (0,0,0)
    # print(maximum)
    # print(minimum)
    voxel_size = mrc.voxel_size.x/10 # nm
    # print([round(i) for i in (maximum/voxel_size)])
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

    mask = (mrc.data == label).astype(int)
    dist_from_mask = distance_transform_edt(mask == 0)
    for i in range(surf.GetNumberOfCells()):
        # Check if all points which made up the polygon are in the mask
        points_cell = surf.GetCell(i).GetPoints()
        count = 0
        for j in range(0, points_cell.GetNumberOfPoints()):
            x,y,z = points_cell.GetPoint(j)
            x = int(round(x/voxel_size))
            y = int(round(y/voxel_size))
            z = int(round(z/voxel_size))
            # print(x,y,z)
            if z>=maximum[2] or z<=minimum[2] or x>=maximum[0] or x<=minimum[0] or y>=maximum[1] or y<=minimum[1]:
                count += 1
                continue
            if (dist_from_mask[
                    z, y, x] >
                    mask_tolerance_distance):
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






if __name__=="__main__":
    convert()
