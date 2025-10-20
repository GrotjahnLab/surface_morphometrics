#! /usr/bin/env python
"""Measure inter-surface distances and relative orientations between two surfaces.

Usage: python interdistance_orientation.py surface1.gt surface2.gt [options]
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

from cProfile import label
import click
import numpy as np
from pycurv import  TriangleGraph, io
from graph_tool import load_graph
from scipy.spatial import cKDTree


from intradistance_verticality import export_csv
 

def surface_to_surface(graph1, label1, graph2, label2, orientation=True, save_neighbor_index=True, exportcsv=True):
    """ Measure the nearest distance to graph2 for every point on graph1.
    Optionally also measure the nearest neighbor relative orientation angle.
    Bidirectional - data is written to both graphs, and new VTP files are output.
    
    graph1 (str): filename of triangle graph object (.gt file)
    label1 (str): Segmentation label for graph 1
    graph2 (str): filename of triangle graph object (.gt file)
    label2 (str): Segmentation label for graph 2
    orientation (bool): If True, measure and save the relative orientation of the nearest neighbor
    save_neighbor_index (bool): If True, save the index of the nearest neighbor to the graph
    exportcsv (bool): If True, update the CSV file with the new data.
    """
    surface1 = graph1[:-3] + ".vtp" # remove .gt from filename and replace with vtp
    surface2 = graph2[:-3] + ".vtp" # remove .gt from filename and replace with vtp
    print(f"Loading {graph1}")
    tg1 = TriangleGraph()
    tg1.graph = load_graph(graph1)
    xyz1 = tg1.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    print(f"Loading {graph2}")
    tg2 = TriangleGraph()
    tg2.graph = load_graph(graph2)
    xyz2 = tg2.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    print("Calculating Distances")
    vprop1 = tg1.graph.new_vertex_property("double")
    vprop2 = tg2.graph.new_vertex_property("double")

    tree1 = cKDTree(xyz2)
    mindist1, min_index_1 = tree1.query(xyz1)
    vprop1.a = mindist1
    tg1.graph.vp[label2+"_dist"] = vprop1

    tree2 = cKDTree(xyz1)
    mindist2, min_index_2 = tree2.query(xyz2)
    vprop2.a = mindist2
    tg2.graph.vp[label1+"_dist"] = vprop2

    if save_neighbor_index: # Argmin gets the first occurence of the minimum
        print("Saving Neighbor Index")
        neighbor1=tg1.graph.new_vertex_property("int")
        neighbor1.a = min_index_1
        tg1.graph.vp[label2+"_neighbor_index"] = neighbor1
        neighbor2=tg2.graph.new_vertex_property("int")
        neighbor2.a = min_index_2
        tg2.graph.vp[label1+"_neighbor_index"] = neighbor2

    if orientation:
        print("Calculating Relative Orientations")
        angles1 =  tg1.graph.vp.n_v.get_2d_array([0,1,2]).transpose()
        angles2 = tg2.graph.vp.n_v.get_2d_array([0,1,2]).transpose()

        neighbor_angles_1to2 = angles2[min_index_1]
        relative_angles_1 = [np.arccos(np.abs(np.dot(angles1[i], neighbor_angles_1to2[i])))*180/np.pi for i in range(len(angles1))]
        ang1 = tg1.graph.new_vertex_property("float")
        ang1.a = relative_angles_1
        tg1.graph.vp[label2+"_orientation"] = ang1
        neighbor_angles_2to1 = angles1[min_index_2]
        relative_angles_2 = [np.arccos(np.abs(np.dot(angles2[i], neighbor_angles_2to1[i])))*180/np.pi for i in range(len(angles2))]
        ang2 = tg2.graph.new_vertex_property("float")
        ang2.a = relative_angles_2
        tg2.graph.vp[label1+"_orientation"] = ang2



        
    print(f"Distances from {label1} to {label2} (Min, Mean, Median, Max):")
    print(np.min(vprop1.a), np.mean(vprop1.a),np.median(vprop1.a), np.max(vprop1.a))

    print(f"Distances from {label2} to {label1} (Min, Mean, Median, Max):")
    print(np.min(vprop2.a), np.mean(vprop2.a),np.median(vprop2.a), np.max(vprop2.a))

    # Save updated surface and tranglegraph
    # curvedness = tg1.graph.vp["curvedness_VV"]
    print("Saving out files")
    surf1 = tg1.graph_to_triangle_poly()
    io.save_vtp(surf1, surface1)
    tg1.graph.save(graph1)
    # Save CSV
    if exportcsv:
        csvname = graph1[:-3]+".csv"
        export_csv(tg1,csvname)
    # Save updated surface and tranglegraph
    surf2 = tg2.graph_to_triangle_poly()
    io.save_vtp(surf2, surface2)
    tg2.graph.save(graph2)

    if exportcsv:
        csvname = graph2[:-3]+".csv"
        export_csv(tg2, csvname)

@click.command()
@click.argument('graph1', type=click.Path(exists=True))
@click.argument('label1', type=click.STRING)
@click.argument('graph2', type=click.Path(exists=True))
@click.argument('label2', type=click.STRING)
@click.option('--orientation', is_flag=True, default=True, help='Measure and save the relative orientation of the nearest neighbor')
@click.option('--save_neighbor_index', is_flag=True, default=True, help='Save the index of the nearest neighbor to the graph')
@click.option('--exportcsv', is_flag=True, default=True, help='Update the CSV file with the new data')
def inter_cli(graph1, label1, graph2, label2, orientation, save_neighbor_index, exportcsv):
    """Run intra-surface distance and optionally verticality measurements for a surface"""
    surface_to_surface(graph1, label1, graph2, label2,orientation=orientation, save_neighbor_index=save_neighbor_index, exportcsv=exportcsv)



if __name__ == "__main__":
    inter_cli()

 