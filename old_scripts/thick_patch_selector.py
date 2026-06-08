import pandas as pd
import numpy as np
import starfile as sf
from scipy import spatial
from sys import argv
from glob import glob
import subprocess


from pycurv import  TriangleGraph, io
from graph_tool import load_graph

PATCHSIZE = 7.5
# PATCHTYPE = "high"
NPATCHES = 50
PIXSIZE = 6.65

def get_patches(csv, patchsize, patchtype, npatches):
    """
    Function to get patches from a csv file
    """
    # Read the csv file
    df = pd.read_csv(csv)
    if "thickness" not in df.columns:
        return ([],[],[])
    # Sort the dataframe by thickness from high to low
    peak_triangles = []
    patches = []
    if patchtype == "high":
        ascendingtruth = False
    elif patchtype == "low":
        ascendingtruth = True
    else:
        raise ValueError("Patch type must be either 'high' or 'low'")

    df = df.sort_values(by="thickness", ascending=ascendingtruth)
    while len(peak_triangles)<npatches:
        # Grab the top current triangle
        current_triangle = df.iloc[0]
        peak_triangles.append(current_triangle)
        # Extract X,y,z coordinates and use a kd tree to find the nearest neighbors
        xyz = current_triangle[["xyz_x","xyz_y","xyz_z"]]
        xyzset = df[["xyz_x","xyz_y","xyz_z"]].to_numpy()
        xyztree = spatial.cKDTree(xyzset)
        l,neighbors = xyztree.query(xyz,k=500, distance_upper_bound=patchsize, workers=-1)
        # Weight the neighbors by distance
        # Calculate weights as 1/(1+l)
        neighbors = neighbors[np.where(l != np.inf)]
        if len(neighbors) == 0:
            print("No neighbors found")
            break
        # Add the neighbors to the patches dataframe
        patches.extend(df.index[neighbors].to_list())
        # Remove the neighbors from the dataframe
        df = df.drop(df.index[neighbors])
    
    indices_to_remove = df.index.to_list()
    # indices_to_remove.remove(patches)
    return peak_triangles, patches, indices_to_remove   

def save_particle_picks(triangles, starfile, newxyzfile, modfile, pixsize ):
    """
    Function to convert a list of triangles to a star file
    """
    # return "Not implemented yet"
    with open(newxyzfile, "w") as f:
        for triangle in triangles:
            f.write(f"{triangle['xyz_x']/(pixsize/10):0.2f}\t{triangle['xyz_y']/(pixsize/10):0.2f}\t{triangle['xyz_z']/(pixsize/10):0.2f}\n")
    subprocess.run(["point2model", "-input", newxyzfile, "-output", modfile, "-scat", "-sp", f"{int(PATCHSIZE*20/pixsize)}"])
    star = sf.read("test.star")[0:0]
    for triangle in triangles:
        star.loc[len(star)] = [triangle['xyz_x']/(pixsize/10), triangle['xyz_y']/(pixsize/10), triangle['xyz_z']/(pixsize/10), 0.00, 0.00, 0.00]
    sf.write(star, starfile)
    return


def make_patch_surface(graph_file, outfile, indices_to_remove):
    """
    Function to make a surface from a list of patches

    Parameters
    ----------
    graph_file : str
        The graph file to read
    outfile : str
        The output file to write
    indices_to_remove : list
        The list of indexes to remove
    """
    # Load the triangle graph
    tg = TriangleGraph()
    tg.graph=load_graph(graph_file)
    # Set the graph to remove all but the patches
    tg.graph.remove_vertex(indices_to_remove, fast=True)

    # Write the graph to a file
    surf_clean = tg.graph_to_triangle_poly()
    io.save_vtp(surf_clean, outfile)

if __name__ == "__main__":
    basefolder = argv[1]
    csvs = glob(basefolder+"/*IMM*8.csv")
    for csv in csvs:
        for patch_variety in ["high","low"]:
            print(csv)
            peak_triangles, patches, indices_to_remove = get_patches(csv, patchtype = patch_variety, npatches=NPATCHES, patchsize=PATCHSIZE, )
            if len(peak_triangles) == 0:
                continue
            gtfile = csv[:-4]+".gt"
            vtpfile = csv[:-4]+f"_{patch_variety}_patches.vtp"
            starfile = csv[:-4]+f"_{patch_variety}_patches.star"
            newxyzfile = csv[:-4]+f"_{patch_variety}_patches.txt"
            modfile = csv[:-4]+f"_{patch_variety}_patches.mod"
            plyfile = csv[:-4]+"_patches.ply"
            save_particle_picks(peak_triangles, starfile, newxyzfile, modfile, PIXSIZE)
            make_patch_surface(gtfile, vtpfile, indices_to_remove)
        


    