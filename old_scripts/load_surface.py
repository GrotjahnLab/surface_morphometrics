

from pycurv import run_gen_surface, THRESH_SIGMA1, TriangleGraph, MAX_DIST_SURF, io
import numpy as np
from scipy import ndimage
from graph_tool import load_graph
import time

fold = r"/Users/benjaminbarad/Dropbox (Scripps Research)/Grotjahn Lab/Data/Surfaces/Lam2_TS1/" # output will be also written there
base_filename = "test"
pixel_size = 2.16  # pixel size of the (underlying) segmentation
radius_hit = 10  # radius of the smallest feature of interest (neighborhood)

# alternative or optional:
# for step 1.:
  # for segmentation input:
seg_file = "TS1_surf.labels.mrc" # MRC in this example
label = 1
cube_size = 5 # try 3 or 5
# filled_label = <lumen_label>  # if compartment segmentation
#   for surface input:
# surf_file = "TS1_surf_100k.obj"  # VTP in this example
# for step 2.:
# to remove small disconnected surface components within this size (default 100)
min_component = 100
# for step 3.:
methods = ["VV"]  # list of algorithms to run (default "VV")
area2 = True  # if method "VV": True for AVV (default), False for RVV
cores = 6  # number of cores to run VV in parallel (default 6)
print("Loading Segmentation at {}".format(time.time()))
seg = io.load_tomo(fold + seg_file)
data_type = seg.dtype
print("Making Binary at {}".format(time.time()))
binary_seg = (seg == label).astype(data_type)
cube = np.ones((cube_size, cube_size, cube_size))
print(np.sum(binary_seg > 0))
print("Closing Binary at {}".format(time.time()))
binary_seg = ndimage.binary_closing(
    binary_seg, structure=cube, iterations=1).astype(data_type)
print(np.sum(binary_seg > 0))
surf = run_gen_surface(binary_seg, fold + base_filename, lbl=1)
try:
	io.write_stl(surf,"{}{}.stl".format(fold, base_filename))
except:
	print("Failed to write STL after surf gen")

# surf = io.load_poly(fold + surf_file)
# print(surf)
print("Making TriangleGraph at {}".format(time.time()))
tg = TriangleGraph()
scale = [pixel_size, pixel_size, pixel_size]
tg.build_graph_from_vtk_surface(surf, scale)
print(tg)
tg.find_vertices_near_border(MAX_DIST_SURF * pixel_size, purge=True)
tg.find_small_connected_components(
    threshold=min_component, purge=True, verbose=True)
print('The graph has {} vertices and {} edges'.format(
    tg.graph.num_vertices(), tg.graph.num_edges()))

print("Writing out TriangleGraph at {}".format(time.time()))
clean_graph_file = '{}.scaled_cleaned.gt'.format(base_filename)
clean_surf_file = '{}.scaled_cleaned.vtp'.format(base_filename)
tg.graph.save(fold + clean_graph_file)
surf_clean = tg.graph_to_triangle_poly()
io.save_vtp(surf_clean, fold + clean_surf_file)

try:
	io.write_stl(surf,"{}{}.stl".format(fold, base_filename))
except:
	print("Failed to write STL after graph gen")