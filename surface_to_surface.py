from graph_tool import load_graph
import numpy as np
from pycurv import  TriangleGraph, io
from scipy.spatial.distance import cdist
import time
import pandas
import vtk


folder = "/Users/benjaminbarad/Dropbox (Scripps Research)/Surface_analysis/"
surface1_base= "UF1_OMM.AVV_rh10"
surface2_base = "UF1_ER.AVV_rh10"
surface1_name = "omm"
surface2_name = "er"
save_neighbor_index = True
export_csv=True
# Additional properties to add to CSV (in addition to distance and neighbor index)
double_properties_to_write = ["kappa_1","kappa_2","curvedness_VV", "gauss_curvature_VV", "mean_curvature_VV", "shape_index_VV", "area"]
vector_properties_to_write = ["n_v"]

def save_ply(poly, fname):
    writer = vtk.vtkPLYWriter()    
    writer.SetFileName(fname)
    writer.SetInputData(poly)
    writer.Write()

def surface_to_surface_distance(graph1, surface1, surface1_name, graph2, surface2, surface2_name, save_neighbor_index=True):
    """ Measure the nearest distance to graph2 for every point on graph1.
    Bidirectional - data is written to both graphs, and new VTP files are output.
    


    """
    print("Loading TG1")
    tg1 = TriangleGraph()
    tg1.graph = load_graph(graph1)
    xyz1 = tg1.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    print("Loading TG2")
    tg2 = TriangleGraph()
    tg2.graph = load_graph(graph2)
    xyz2 = tg2.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    # xyz1 = tg1.graph.vp.xyz.get_vertices() # get as array
    print("Calculating Distances")
    dist_matrix = cdist(xyz1, xyz2) # Efficient calculation of distances.
    print("Distances for Object 1 (Min, Mean, Median, Max):")
    vprop1 = tg1.graph.new_vertex_property("double")
    vprop1.a = dist_matrix.min(axis=1)
    print(np.min(vprop1.a), np.mean(vprop1.a),np.median(vprop1.a), np.max(vprop1.a))
    tg1.graph.vp[surface2_name+"_dist"] = vprop1
    if save_neighbor_index: # Argmin gets the first occurence of the minimum
        neighbor=tg1.graph.new_vertex_property("int")
        neighbor.a = dist_matrix.argmin(axis=1)
        tg1.graph.vp[surface2_name+"_neighbor_index"] = neighbor

    print("Distances for Object 2 (Min, Mean, Median, Max):")
    vprop2 = tg2.graph.new_vertex_property("double")
    vprop2.a = dist_matrix.min(axis=0)
    print(np.min(vprop2.a), np.mean(vprop2.a),np.median(vprop2.a), np.max(vprop2.a))
    tg2.graph.vp[surface1_name+"_dist"] = vprop2
    if save_neighbor_index: # Argmin gets the first occurence of the minimum
        neighbor=tg2.graph.new_vertex_property("int")
        neighbor.a = dist_matrix.argmin(axis=0)
        tg2.graph.vp[surface1_name+"_neighbor_index"] = neighbor

    # Save updated surface and tranglegraph
    curvedness = tg1.graph.vp["curvedness_vv"]
    surf1 = tg1.graph_to_triangle_poly()
    io.save_vtp(surf1, surface1)
    ply_name = graph1[:-3]+".ply"
    save_ply(surf1, ply_name, curvedness)
    tg1.graph.save(graph1)
    # Save CSV
    if export_csv:
        csvname = graph1[:-3]+".csv"
        df = pandas.DataFrame()
        for vertex_property in double_properties_to_write:
            df[vertex_property] = tg1.get_vertex_property_array(vertex_property)
        df[surface2_name+"_dist"] = tg1.get_vertex_property_array(surface2_name+"_dist")
        df[surface2_name+"_neighbor_index"] = tg1.get_vertex_property_array(surface2_name+"_neighbor_index")
        for vector_property in vector_properties_to_write:
            x,y,z = tg1.graph.vp[vector_property].get_2d_array([0,1,2])
            df[vector_property+"_x"] = x
            df[vector_property+"_y"] = y
            df[vector_property+"_z"] = z
        df.to_csv(csvname, index_label="index")
    # Save updated surface and tranglegraph
    surf2 = tg2.graph_to_triangle_poly()
    io.save_vtp(surf2, surface2)
    ply_name = graph2[:-3]+".ply"
    save_ply(surf2, ply_name)
    tg2.graph.save(graph2)

    if export_csv:
        csvname = graph2[:-3]+".csv"
        df = pandas.DataFrame()
        for vertex_property in double_properties_to_write:
            df[vertex_property] = tg2.get_vertex_property_array(vertex_property)
        df[surface1_name+"_dist"] = tg2.get_vertex_property_array(surface1_name+"_dist")
        df[surface1_name+"_neighbor_index"] = tg2.get_vertex_property_array(surface1_name+"_neighbor_index")
        for vector_property in vector_properties_to_write:
            x,y,z = tg2.graph.vp[vector_property].get_2d_array([0,1,2])
            df[vector_property+"_x"] = x
            df[vector_property+"_y"] = y
            df[vector_property+"_z"] = z
        df.to_csv(csvname, index_label="index")

# def main():

if __name__=="__main__":
    graph1 = folder+surface1_base+".gt"
    surface1 = folder+surface1_base+".vtp"
    graph2 = folder+surface2_base+".gt"
    surface2 = folder+surface2_base+".vtp"
    surface_to_surface_distance(graph1, surface1, surface1_name, graph2, surface2, surface2_name, save_neighbor_index=save_neighbor_index)
