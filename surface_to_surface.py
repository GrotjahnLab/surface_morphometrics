from graph_tool import load_graph
import numpy as np
from pycurv import  TriangleGraph, io
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import time
import pandas
import vtk

from sys import argv
from os import path
from time import sleep

folder = "/gpfs/group/grotjahn/bbarad/Final_Dataset/"
rh_level = 12
basename = argv[1][:-11]

print(basename)

# folder = "/Users/benjaminbarad/Dropbox (Scripps Research)/Surface_analysis/"
# surface1_base= "UF1_OMM.AVV_rh10"
# surface2_base = "UF1_ER.AVV_rh10"
# surface1_name = "omm"
# surface2_name = "er"
save_neighbor_index = True
export_csv=True
# Additional properties to add to CSV (in addition to distance and neighbor index)
# double_properties_to_write = ["kappa_1","kappa_2","curvedness_VV", "gauss_curvature_VV", "mean_curvature_VV", "shape_index_VV", "area"]
# vector_properties_to_write = ["n_v"]

def save_ply(poly, fname):
    writer = vtk.vtkPLYWriter()    
    writer.SetFileName(fname)
    writer.SetInputData(poly)
    writer.Write()

def get_dist_two_directions(point, normal, locator, xyz, dist_min=3, dist_max=400, tolerance=0.1, ):
    """Returns the distance and cell ID from a certain point along both the 
    positive and normal axis

    point: array of the point
    normal: array of the normal - recommended to use voted normals.
    locator: vtk cell locator object
    dist_min: minimum distance to consider (to avoid self-intersection)
    dist_max: maximum distance to consider
    tolerance: tolerance for the distance calculation
    """
    positions = []
    distances = []
    cell_ids = []
    for direction in [-normal, normal]:
        
        p0 = point + direction*dist_min
        pmax = point + direction*dist_max
        p1 = [0.,0.,0.]
        pcoords = [0,0,0]
        sub_id = vtk.mutable(-1)
        cell1_id = vtk.mutable(-1)
        t = vtk.mutable(0.) # parametric coordinate of the intersection - 0 to 1 if intersection is found
        locator.IntersectWithLine(p0, pmax, tolerance, t, p1, pcoords, sub_id,
                              cell1_id)
        distance = np.linalg.norm(p1-point)
        if cell1_id == -1:
            distances.append(np.nan)
            cell_ids.append(-1)
        elif distance > dist_max+1 or distance < dist_min-1: # Tolerance adjustment
            distances.append(np.nan)
            cell_ids.append(-1)       
        else:
            distances.append(np.linalg.norm(p1-point)) # distance in the direction of the normal
            cell_ids.append(int(cell1_id))
    # switch orders if the first distance is larger:
    if np.isnan(distances[0]):
        distances = distances[::-1]
        cell_ids = cell_ids[::-1]
    elif np.isnan(distances[1]):
        # Only one distance
        pass
    elif distances[0] > distances[1]:
        # print("Distances misordered, flipping")
        distances = distances[::-1]
        cell_ids = cell_ids[::-1]
    
    # close dist, close id, far dist, far id
    assert np.isnan(distances[1]) or distances[0]<=distances[1]
    
    return distances[0], cell_ids[0], distances[1], cell_ids[1]

def surface_self_distances(graph_file, surface_file, dist_min=6, dist_max=200, tolerance=0.1, exportcsv=True):
    """Returns the distances between all vertices of two surfaces - 
    inspired by find_1_distance in pycurv

    graph_file: graph-tool graph filename of the surface
    surface_file: vtk surface filename of the first surface
    dist_min: minimum distance to consider (to avoid self-intersection)
    dist_max: maximum distance to consider (stick to physiological relevant distances)
    tolerance: tolerance for the distance calculation (recommended to be 0.1)
    """
    # Initialize stuff
    print("Loading graph and surface")
    surface = io.load_poly(surface_file)
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()
    tg = TriangleGraph()
    tg.graph = load_graph(graph_file)
    xyz = tg.graph.vp.xyz.get_2d_array([0,1,2]).transpose()
    normal = tg.graph.vp.n_v.get_2d_array([0,1,2]).transpose() # use n_v for voted normals
    # Initialize variables
    print("Initializing variables")
    close_distances = tg.graph.new_vertex_property("float")
    close_id = tg.graph.new_vertex_property("long")
    far_distances = tg.graph.new_vertex_property("float")
    far_id = tg.graph.new_vertex_property("long")
    # Vectorized distance processor
    # Calculate distances
    print("Calculating distances")
    for i in range(len(close_distances.a)):
        close_distances.a[i], close_id.a[i], far_distances.a[i], far_id.a[i] = get_dist_two_directions(xyz[i], normal[i], locator, xyz, dist_min, dist_max, tolerance=0.001)
    # Write out distances
    print(np.nanmin(close_distances.a), np.nanmax(close_distances.a))
    print("Writing out distances")
    tg.graph.vp.self_dist_min = close_distances
    tg.graph.vp.self_id_min = close_id
    tg.graph.vp.self_dist_far = far_distances
    tg.graph.vp.self_id_far = far_id
    # Save graph
    tg.graph.save(graph_file)
    # Save CSV with all features
    if exportcsv:
        csvname = graph_file[:-3]+".csv"
        print(csvname)
        export_csv(tg, csvname)
    return tg


def surface_to_surface_distance(graph1, surface1, surface1_name, graph2, surface2, surface2_name, save_neighbor_index=True, exportcsv=True):
    """ Measure the nearest distance to graph2 for every point on graph1.
    Bidirectional - data is written to both graphs, and new VTP files are output.
    


    """
    print(f"Loading {surface1_name}")
    tg1 = TriangleGraph()
    tg1.graph = load_graph(graph1)
    xyz1 = tg1.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    print(f"Loading {surface2_name}")
    tg2 = TriangleGraph()
    tg2.graph = load_graph(graph2)
    xyz2 = tg2.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    print("Calculating Distances")
    vprop1 = tg1.graph.new_vertex_property("double")
    vprop2 = tg2.graph.new_vertex_property("double")
    datasize = len(xyz2)*len(xyz1)*64 # approximate size of a distnace matrix
    print(f"Estimated memory usage for cdist matrix: {datasize/1e9} GB")

    ## If the size of the matrix is small, use a cdist matrix.
    if datasize<128e9: # 128GB of memory
        print(f"Using cdist matrix")
        dist_matrix = cdist(xyz1, xyz2) # Efficient calculation of distances.
        vprop1.a = dist_matrix.min(axis=1)
        tg1.graph.vp[surface2_name+"_dist"] = vprop1

        vprop2.a = dist_matrix.min(axis=0)
        tg2.graph.vp[surface1_name+"_dist"] = vprop2
        if save_neighbor_index: # Argmin gets the first occurence of the minimum
            min_index_1 = np.argmin(dist_matrix, axis=1)
            min_index_2 = np.argmin(dist_matrix, axis=0)
    ## Otherwise, make your calculations using a cKDTree
    else:
        print("Using KDTree")
        tree1 = cKDTree(xyz2)
        mindist1, min_index_1 = tree1.query(xyz1)
        vprop1.a = mindist1
        tg1.graph.vp[surface2_name+"_dist"] = vprop1

        tree2 = cKDTree(xyz1)
        mindist2, min_index_2 = tree2.query(xyz2)
        vprop2.a = mindist2
        tg2.graph.vp[surface1_name+"_dist"] = vprop2
    
    if save_neighbor_index: # Argmin gets the first occurence of the minimum
        neighbor1=tg1.graph.new_vertex_property("int")
        neighbor1.a = min_index_1
        tg1.graph.vp[surface2_name+"_neighbor_index"] = neighbor1
        neighbor2=tg2.graph.new_vertex_property("int")
        neighbor2.a = min_index_2
        tg2.graph.vp[surface1_name+"_neighbor_index"] = neighbor2
    print(f"Distances from {surface1_name} to {surface2_name} (Min, Mean, Median, Max):")
    print(np.min(vprop1.a), np.mean(vprop1.a),np.median(vprop1.a), np.max(vprop1.a))

    print(f"Distances from {surface2_name} to {surface1_name} (Min, Mean, Median, Max):")
    print(np.min(vprop2.a), np.mean(vprop2.a),np.median(vprop2.a), np.max(vprop2.a))

    # Save updated surface and tranglegraph
    # curvedness = tg1.graph.vp["curvedness_VV"]
    print("Saving out files")
    surf1 = tg1.graph_to_triangle_poly()
    io.save_vtp(surf1, surface1)
    ply_name = graph1[:-3]+".ply"
    save_ply(surf1, ply_name)
    tg1.graph.save(graph1)
    # Save CSV
    if exportcsv:
        csvname = graph1[:-3]+".csv"
        export_csv(tg1,csvname)
    # Save updated surface and tranglegraph
    surf2 = tg2.graph_to_triangle_poly()
    io.save_vtp(surf2, surface2)
    ply_name = graph2[:-3]+".ply"
    save_ply(surf2, ply_name)
    tg2.graph.save(graph2)

    if exportcsv:
        csvname = graph2[:-3]+".csv"
        export_csv(tg2, csvname)

def export_csv(tg, csvname):
    """Export all properties in a triangle graph to a CSV with indices for quick visualization with other programs"""
    print(f"Exporting graph data to CSV {csvname}")
    df = pandas.DataFrame()
    vector_properties = []
    scalar_properties = []
    properties = tg.graph.vertex_properties.keys()
    
    ## defined scalar types from the graph-tool documentation
    scalars = ["bool", "int16_t", "int32_t", "int64_t", "unsigned long","double", "long double"]
    
    ## iterative over properties and classify them as scalars or vectors
    for property in properties:
        type = tg.graph.vertex_properties[property].value_type()
        if type in scalars:
            scalar_properties.append(property)
        elif type[:6] == "vector":
            vector_properties.append(property)
    print("Scalar Properties to save: ", scalar_properties)
    print("Vector Properties to save: ", vector_properties)
    for property in scalar_properties:
        df[property] = tg.get_vertex_property_array(property)
    for vector_property in vector_properties:
        x,y,z = tg.graph.vp[vector_property].get_2d_array([0,1,2])
        df[vector_property+"_x"] = x
        df[vector_property+"_y"] = y
        df[vector_property+"_z"] = z
    df.to_csv(csvname, index_label="index")     

# def main():
if __name__=="__main__":
    omm_graph = folder+basename+"_OMM.AVV_rh{}.gt".format(rh_level)
    omm_surface = folder+basename+"_OMM.AVV_rh{}.vtp".format(rh_level)
    imm_graph = folder+basename+"_IMM.AVV_rh{}.gt".format(rh_level)
    imm_surface = folder+basename+"_IMM.AVV_rh{}.vtp".format(rh_level)
    er_graph = folder+basename+"_ER.AVV_rh{}.gt".format(rh_level)
    er_surface = folder+basename+"_ER.AVV_rh{}.vtp".format(rh_level)
    if path.isfile(imm_graph):
        print("Calculating IMM-IMM distances")
        surface_self_distances(imm_graph, imm_surface, dist_min=4, dist_max=400, tolerance=0.001)
    if path.isfile(omm_graph):
        if path.isfile(imm_graph):
            print("Measuring OMM-IMM Distances")
            surface_to_surface_distance(omm_graph, omm_surface, "OMM", imm_graph, imm_surface, "IMM")
        else:
            print("No IMM Found - Skipping OMM-IMM Distances")
        
        if path.isfile(er_graph):
            print("Measuring OMM-ER Distances")
            surface_to_surface_distance(omm_graph, omm_surface, "OMM", er_graph, er_surface, "ER")
        else:
            print("No ER Found - Skipping OMM-ER Distances")
    else:
        print("No OMM Found - skipping all measurements")
    # graph1 = folder+surface1_base+".gt"
    # surface1 = folder+surface1_base+".vtp"
    # graph2 = folder+surface2_base+".gt"
    # surface2 = folder+surface2_base+".vtp"
    # surface_to_surface_distance(graph1, surface1, surface1_name, graph2, surface2, surface2_name, save_neighbor_index=save_neighbor_index)
