from curvature_calculation import from_ply_workflow, new_workflow, extract_curvatures_after_new_workflow
import time

t_begin = time.time()
fold = r"/gpfs/group/grotjahn/bbarad/poisson_test/" # output will be also written there
base_filename = "fission.ply"
surf_file = "fission.ply"
#seg_file = ""
# surf_file = 
lbl=1
holes=5
pixel_size = 1.42
radius_hit=8
min_component = 100
exclude_borders=1
runtimes_file = "{}{}_runtimes.csv".format(fold, base_filename)
print("\nCalculating curvatures for {}".format(base_filename))
#new_workflow(
#    base_filename, seg_file, fold, pixel_size, radius_hit, methods=['VV'],
#    label=lbl, holes=holes, min_component=min_component, runtimes=runtimes_file)
#print("\nExtracting curvatures for all surfaces")
#extract_curvatures_after_new_workflow(
#    fold, base_filename, radius_hit, methods=['VV'],
#    exclude_borders=exclude_borders, categorize_shape_index=True)
#binary_seg = (seg == label).astype(data_
from_ply_workflow(fold+surf_file, radius_hit)
t_end = time.time()
duration = t_end - t_begin
minutes, seconds = divmod(duration, 60)
print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))
