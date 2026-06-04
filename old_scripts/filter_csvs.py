import os
import mrcfile
import pandas as pd
import shutil
from scipy.ndimage import distance_transform_edt
import numpy as np


directory = "/gpfs/group/grotjahn/bbarad/Final_Dataset/TgGSK/"
datasets = ["GE1", "GE2", "GF1", "GF2", "GT1", "GT2", "GT3","GT4","GT5"]
# datasets = ["GE1"]
postfix = ".AVV_rh12.csv"
for dataset in datasets:
    mrc = directory+dataset+"_labels.mrc"
    with mrcfile.open(mrc, 'r', permissive=True) as file:
        voxel_size = file.voxel_size.x/10
        data = file.data
        print(data.shape)
    
    dist_from_data = distance_transform_edt(data == 0)
    for surface in ["OMM", "IMM", "ER"]:
        filename = directory+dataset+"_"+surface+postfix
        print(filename)
        shutil.copyfile(filename, filename+".bak")
        df = pd.read_csv(filename)
        print(df)
        xyz = np.array([df["xyz_z"], df["xyz_y"], df["xyz_x"]]).transpose()
        xyz = xyz/voxel_size
        xyz = np.around(xyz).astype(int)
        mask = np.zeros(xyz.shape[0])
        for i in range(xyz.shape[0]):
            x,y,z = xyz[i]
            distance = dist_from_data[x,y,z]
            mask[i] = distance < 1

        new_df = df[mask.astype(bool)]
        new_df.to_csv(filename)

