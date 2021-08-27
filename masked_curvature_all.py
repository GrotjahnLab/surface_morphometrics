from curvature_calculation import from_ply_workflow, new_workflow, extract_curvatures_after_new_workflow
import time
from sys import argv

t_begin = time.time()
fold = r"/gpfs/group/grotjahn/bbarad/Surface_analysis/" # output will be also written there

index = int(argv[1])
filelist = ["TF1_IMM","TF1_OMM","TF2_IMM","TF2_OMM","TT2_IMM","TT2_OMM","TT3_IMM","TT3_OMM","TE1_IMM","TE1_OMM","TE2_IMM","TE2_OMM","UF1_IMM","UF1_OMM","UF2_IMM","UF2_OMM","UT1_IMM","UT1_OMM","UT2_IMM","UT2_OMM","UE1_IMM","UE1_OMM","UE2_IMM","UE2_OMM"]
base_filename = filelist[index]
#surf_file = "TE2_OMM.ply"
seg_file = ""
# surf_file = 
#lbl=1
#holes=5
pixel_size = 1
radius_hit=10
min_component = 20
exclude_borders=1
runtimes_file = "{}{}_runtimes.csv".format(fold, base_filename)
print("\nCalculating curvatures for {}".format(base_filename))
new_workflow(
    base_filename, seg_file, fold, pixel_size, radius_hit, methods=['VV'],  min_component=min_component, runtimes=runtimes_file)
print("\nExtracting curvatures for all surfaces")
extract_curvatures_after_new_workflow(
    fold, base_filename, radius_hit, methods=['VV'],
    exclude_borders=exclude_borders, categorize_shape_index=True)
#binary_seg = (seg == label).astype(data_
#from_ply_workflow(fold+surf_file, radius_hit)
t_end = time.time()
duration = t_end - t_begin
minutes, seconds = divmod(duration, 60)
print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))
