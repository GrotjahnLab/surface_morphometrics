from curvature_calculation import from_ply_workflow, new_workflow, extract_curvatures_after_new_workflow
import time
from sys import argv
import mrcfile

files = [("Treated", "TF1_IMM"), # Gotta rerun this bad boy
("TgGSK","GE1_ER"),  
("TgGSK","GE2_OMM"),  
("TgGSK","GF2_IMM"),  
("TgGSK","GT2_ER"),  
("TgGSK","GT3_OMM"),  
("TgGSK","GT5_IMM"),
("TgGSK","GE1_IMM"),  
("TgGSK","GF1_ER"),  
("TgGSK","GF2_OMM"),  
("TgGSK","GT2_IMM"),  
("TgGSK","GT4_ER"),  
("TgGSK","GT5_OMM"),
("TgGSK","GE1_OMM"),  
("TgGSK","GF1_IMM"),  
("TgGSK","GT1_ER"),  
("TgGSK","GT2_OMM"),  
("TgGSK","GT4_IMM"),
("TgGSK","GE2_ER"),  
("TgGSK","GF1_OMM"),  
("TgGSK","GT1_IMM"),  
("TgGSK","GT3_ER"),  
("TgGSK","GT4_OMM"),
("TgGSK","GE2_IMM"),  
("TgGSK","GF2_ER"),  
("TgGSK","GT1_OMM"),  
("TgGSK","GT3_IMM"),  
("TgGSK","GT5_ER"),
("Treated", "TT9_ER"),
("Treated", "TT9_IMM"),
("Treated", "TT9_OMM")]
# files = [("Treated","TT1_ER"),("Treated","TT1_IMM"),("Treated","TT1_OMM"),("Treated","TE2_ER"),
# ("Treated","TE2_IMM"),("Treated","TE2_OMM"),("Treated","TF5_ER"),("Treated","TF5_IMM"),("Treated","TF5_OMM"),
# ("Treated","TE5_ER"),("Treated","TE5_IMM"),("Treated","TE5_OMM"),("Treated","TF2_ER"),("Treated","TF2_IMM"),
# ("Treated","TF2_OMM"),("Treated","TT6_ER"),("Treated","TT6_IMM"),("Treated","TT6_OMM"),("Treated","TT3_IMM"),
# ("Treated","TT3_OMM"),("Treated","TE3_IMM"),("Treated","TE3_OMM"),("Treated","TT5_ER"),("Treated","TT5_IMM"),
# ("Treated","TT5_OMM"),("Treated","TF1_ER"),("Treated","TF1_IMM"),("Treated","TF1_OMM"),("Treated","TE6_ER"),
# ("Treated","TE6_IMM"),("Treated","TE6_OMM"),("Treated","TT8_IMM"),("Treated","TT8_OMM"),("Treated","TF6_ER"),
# ("Treated","TF6_IMM"),("Treated","TF6_OMM"),("Treated","TE1_ER"),("Treated","TE1_IMM"),("Treated","TE1_OMM"),
# ("Treated","TT2_IMM"),("Treated","TT2_OMM"),("Treated","TT7_ER"),("Treated","TT7_IMM"),("Treated","TT7_OMM"),
# ("Treated","TE4_ER"),("Treated","TE4_IMM"),("Treated","TE4_OMM"),("Treated","TF3_IMM"),("Treated","TF3_OMM"),
# ("Untreated","UT1_ER"),("Untreated","UT1_IMM"),("Untreated","UT1_OMM"),("Untreated","UE2_ER"),
# ("Untreated","UE2_IMM"),("Untreated","UE2_OMM"),("Untreated","UF5_IMM"),("Untreated","UF5_OMM"),
# ("Untreated","UT4_ER"),("Untreated","UT4_IMM"),("Untreated","UT4_OMM"),("Untreated","UT3_ER"),
# ("Untreated","UT3_IMM"),("Untreated","UT3_OMM"),("Untreated","UE5_ER"),("Untreated","UE5_IMM"),
# ("Untreated","UE5_OMM"),("Untreated","UF2_ER"),("Untreated","UF2_IMM"),("Untreated","UF2_OMM"),
# ("Untreated","UT5_ER"),("Untreated","UT5_IMM"),("Untreated","UT5_OMM"),("Untreated","UF1_ER"),
# ("Untreated","UF1_IMM"),("Untreated","UF1_OMM"),("Untreated","UE3_ER"),("Untreated","UE3_IMM"),
# ("Untreated","UE3_OMM"),("Untreated","UF4_IMM"),("Untreated","UF4_OMM"),("Untreated","UE4_ER"),
# ("Untreated","UE4_IMM"),("Untreated","UE4_OMM"),("Untreated","UF3_ER"),("Untreated","UF3_IMM"),
# ("Untreated","UF3_OMM"),("Untreated","UF6_ER"),("Untreated","UF6_IMM"),("Untreated","UF6_OMM"),
# ("Untreated","UE1_IMM"),("Untreated","UE1_OMM"),("Untreated","UT2_ER"),("Untreated","UT2_IMM"),
# ("Untreated","UT2_OMM")] # yowza
argument = int(argv[1])
treatment, basename = files[argument]
t_begin = time.time()
fold = r"/gpfs/group/grotjahn/bbarad/Final_Dataset/"+treatment+"/" # Treated or Untreated
mrc_filename = basename.split("_")[0]+"_labels.mrc"
with mrcfile.open(fold+"/"+mrc_filename, mode="r", permissive=True) as mrc:
    pixel_size =  mrc.voxel_size.x/10

#surf_file = "TE2_OMM.ply"
seg_file = ""

# surf_file = 
#lbl=1
#holes=5
# pixel_size = 1
radius_hit=12
min_component = 20
exclude_borders=1
runtimes_file = "{}{}_runtimes.csv".format(fold, basename)
print("\nCalculating curvatures for {}".format(basename))
new_workflow(
    basename, seg_file, fold, pixel_size, radius_hit, methods=['VV'],  min_component=min_component, runtimes=runtimes_file)
print("\nExtracting curvatures for all surfaces")
extract_curvatures_after_new_workflow(
    fold, basename, radius_hit, methods=['VV'],
    exclude_borders=exclude_borders, categorize_shape_index=True)
#binary_seg = (seg == label).astype(data_
#from_ply_workflow(fold+surf_file, radius_hit)
t_end = time.time()
duration = t_end - t_begin
minutes, seconds = divmod(duration, 60)
print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))
