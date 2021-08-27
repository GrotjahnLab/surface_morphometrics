from distances_calculation import *

fold = r"/Users/benjaminbarad/Dropbox (Scripps Research)/Grotjahn Lab/Data/Surfaces/Lam2_TS1/" # output will be also written there
base_filename = "distance_calculation"
segmentation_file = "TS1_surf.labels.mrc"
lbl_mem1 = 1
lbl_mem2 = 3

log_file = '{}{}.distances_and_thicknesses_calculation.log'.format(
        fold, base_filename)
    sys.stdout = open(log_file, 'a')

distances_and_thicknesses_calculation(
        fold, segmentation_file, base_filename, 
        lbl_mem1=lbl_mem1, lbl_mem2=lbl_mem2, lbl_between_mem1_mem2=None,
        pixel_size=1.44, radius_hit=10, maxdist=50, maxthick=1,
        offset_voxels=0, both_directions=True, reverse_direction=False,
        mem1="IMM", mem2="OMM", smooth=False)