import pandas as pd
import numpy as np
import mrcfile
import scipy.interpolate as interp
from glob import glob
from pathlib import Path
from multiprocessing import Pool

# Constants
base_folder = "/Users/bbarad/Downloads/TE/"
mrcfolder = base_folder + "tomo2/"
workfolder = base_folder + "morphometrics/"
# components = ["OMM", "IMM", "ER"]


## Load data from a CSV file using pandas and extract the xyz coordinates and n_v normal values
#  @param filename Name of the CSV file
#  @param voxsize The voxel size of the mrc data
#  @return x,y,z,n_v
def load_csv(filename, voxsize, origin=(0,0,0)):
    df = pd.read_csv(filename)
    origin = [0,0,0]
    x = (np.array(df['xyz_x'])-origin[0])/voxsize
    y = (np.array(df['xyz_y'])-origin[1])/voxsize
    z = (np.array(df['xyz_z'])-origin[2])/voxsize
    xyz = np.array([x,y,z])
    n_v = np.array([df['n_v_x']/voxsize,df['n_v_y']/voxsize,df['n_v_z']/voxsize])
    return xyz,n_v

## Load mrc data from an mrc file using mrcfile and extract the data to a numpy array and the pixel size
#  @param filename Name of the mrc file
#  @return data
def load_mrc(filename, angstroms=False):
    with mrcfile.open(filename, permissive=True) as mrc:
        print(mrc.header.origin.x, mrc.header.origin.y, mrc.header.origin.z)
        if angstroms:
            origin = (mrc.header.origin.x, mrc.header.origin.y, mrc.header.origin.z)
            voxsize = mrc.voxel_size.x

        else:
            origin = (mrc.header.origin.x/10, mrc.header.origin.y/10, mrc.header.origin.z/10)
            voxsize = mrc.voxel_size.x/10 # Convert from Angstroms to nm

        print(voxsize)
        data = mrc.data
        data = np.swapaxes(data,0,2)
        # data = np.flip(data, axis=2)
        print(data.shape)
        data_matrix = (np.arange(data.shape[0]),np.arange(data.shape[1]),np.arange(data.shape[2]))
    return data,data_matrix, voxsize, origin

## Stepping from -10 nm to 10 nm in parameterized nm steps, interpolate the values of the mrc data along each normal vector using scipy.interpn
# @param data The mrc data
# @param voxsize The voxel size of the mrc data
# @param xyz The xyz coordinates of the faces
# @param n_v The normal vectors 
# @param nm The number of nm steps
# @return value_array
def interpolate(data,data_matrix, xyz,n_v,nsamples=81, angstroms=False):
    averages = []
    # Create an array of nm steps from -10 to 10
    samples = np.linspace(-10,10, nsamples)
    if angstroms:
        samples = samples*10.
    # Create an empty array to store the interpolated values
    value_array = np.empty((len(n_v[0]),len(samples)))
    # Iterate through the normal vectors
    for i in range(len(n_v[0])):
        # print(i)
        # Create an empty array to store the interpolated values for each normal vector
        value_array_temp = np.array((samples))
        # Iterate through the nm steps
            # Interpolate the mrc data along the normal vector
            # skip = False

        locindices = [xyz[:,i]+j*n_v[:,i] for j in samples]
            # for k in [0,1,2]:
            #     if locindex[k]>(data.shape[k]-1):
            #         print(f"Out of bounds: {locindex}")
            #         value_array_temp[idx] = np.nan
            #         skip = True
            # if not skip:
        value_array_temp = interp.interpn(data_matrix,data,locindices, method="linear", bounds_error=False, fill_value=None)
        averages.append(value_array_temp[int((nsamples+1)/2)])
        # Store the interpolated values for each normal vector
        # print(value_array_temp)
        value_array[i] = value_array_temp
    print(np.mean(averages))
    print(value_array.shape)
    print(xyz.shape)
    return value_array


def run_mrc(filename):
    # components = ["OMM"]
    # components = ["ER"]
    # mrcbase = Path(filename).stem
    # print(mrcbase)
    # Load the mrc data
    mrcbase = filename.split(".mrc")[0].split("/")[-1]
    print(mrcbase)
    files = glob(workfolder+mrcbase+f"*.AVV_rh*.csv")
    print(files)
    data,data_matrix,voxsize, origin = load_mrc(filename)
    print(data.shape)

    # Load the xyz coordinates and normal vectors
    for file in files:
        print(file)
        xyz,n_v = load_csv(file,voxsize, origin)
        # Interpolate the mrc data along the normal vectors
        value_array = interpolate(data,data_matrix,xyz,n_v)
        # Save the interpolated values to a csv file
        print(file[:-4] + f"_sampling.csv")
        np.savetxt(file[:-4] + f"_sampling.csv",value_array,delimiter=",")

if __name__ == "__main__":
    # mrc = mrcfolder+"lam6_ts_003.mrc_13.30Apx_flipx.mrc"
    mrcs = glob(mrcfolder+"*.mrc")
    print(mrcs)
    # pool = Pool()
    # a = pool.map(run_mrc, mrcs)
    # pool.close()
    for mrc in mrcs:
        run_mrc(mrc)

    # run_mrc(mrc)



    
