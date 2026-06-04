# Takes a folder with gt files and loads each one, calculates the connected components, and saves the results in a new graph-tool label, then saves the graph and exports a csv and a vtp file.

import glob
from pycurv import  TriangleGraph, io
from graph_tool import load_graph
from graph_tool.topology import label_components
from intradistance_verticality import export_csv


folder = "/Users/bbarad/Downloads/OneDrive_1_5-12-2025/morphometrics2/"

fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

graph_files = glob.glob(folder+"*AVV*.gt")
print(graph_files)
for graph_file in graph_files:
    print(f"Processing {graph_file}")
    tg = TriangleGraph()
    tg.graph=load_graph(graph_file)## Extract the set of connected components
    unique_components,_ = label_components(tg.graph, directed=False)
    tg.graph.vp.unique_component = unique_components
    # Save files
    tg.graph.save(graph_file)
    surf = tg.graph_to_triangle_poly()
    surface_file = graph_file.replace(".gt", ".vtp")
    csv_outfile = graph_file.replace(".gt", ".csv")
    io.save_vtp(surf, surface_file)
    export_csv(tg, csv_outfile)
