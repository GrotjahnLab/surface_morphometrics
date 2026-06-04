# Load a triangle graph, and calculate the distance to the edge of the graph for every triangle. 
# Use the pycurv triangle graph tools to identify the edge, the calculate minimum geodesic distance
# from each triangle to the edge.


from pycurv import TriangleGraph, io
from graph_tool import load_graph
from scipy.spatial import cKDTree
import numpy as np
import os

from intradistance_verticality import export_csv


gtfile_folder = "/Users/bbarad/Downloads/TE/morphometrics/"
gtfiles = [os.path.join(gtfile_folder, f) for f in os.listdir(gtfile_folder) if f.endswith("AVV_rh8.gt")]
print(gtfiles[0])
filter_distance = 8

for gtfile in gtfiles:
    print(gtfile)
    # Load the triangle graph
    tg = TriangleGraph()
    tg.graph = load_graph(gtfile)
    border_vertices_indices = tg.find_graph_border()
    xyz = tg.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    xyz_edge  = xyz[border_vertices_indices]
    # Calculate the distance from each triangle to the edge using a KD tree
    tree = cKDTree(xyz_edge)
    distances, indices = tree.query(xyz)
    distprop = tg.graph.new_vertex_property("double")
    distprop.a = distances
    tg.graph.vp["edge_dist"] = distprop
    tg.graph.save(gtfile)

    # export other files
    vtpfile = gtfile.replace(".gt", ".vtp")
    surf = tg.graph_to_triangle_poly()
    io.save_vtp(surf, vtpfile)

    csvfile = gtfile.replace(".gt", ".csv")
    export_csv(tg,csvfile)


    # Filter the graph to remove triangles that are too close to the edge
    print(np.where(distances > filter_distance))
    tg.graph.remove_vertex(np.where(distances < filter_distance))
    
    gtfile2 = gtfile.replace(".gt", "_edgefiltered.gt")
    tg.graph.save(gtfile2)
    surf = tg.graph_to_triangle_poly()
    vtpfile2 = gtfile2.replace(".gt", ".vtp")
    io.save_vtp(surf, vtpfile2)
    csvfile2 = gtfile2.replace(".gt", ".csv")
    export_csv(tg,csvfile2)

