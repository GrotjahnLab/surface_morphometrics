from curvature_calculation import new_workflow, extract_curvatures_after_new_workflow
import time

t_begin = time.time()
fold = r"/Users/benjaminbarad/Dropbox (Scripps Research)/Grotjahn Lab/Data/Surfaces/Lam2_TS1/" # output will be also written there
base_filename = "all"
# seg_file = "TS1_surf.labels.mrc"
seg_file = ""
# surf_file = 
lbl=3
holes=0
pixel_size = 1.42
radius_hit=10
min_component = 100
runtimes_file = "{}{}_runtimes.csv".format(fold, base_filename)
print("\nCalculating curvatures for {}".format(base_filename))
new_workflow(
    base_filename, seg_file, fold, pixel_size, radius_hit, methods=['VV'],
    label=lbl, holes=0, min_component=min_component, runtimes=runtimes_file)
print("\nExtracting curvatures for all surfaces")
extract_curvatures_after_new_workflow(
    fold, base_filename, radius_hit, methods=['VV'],
    exclude_borders=1, categorize_shape_index=True)




t_end = time.time()
duration = t_end - t_begin
minutes, seconds = divmod(duration, 60)
print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))