"""
Density sampling for surface morphometrics pipeline.

This script samples tomogram density values along normal vectors for each triangle
in a surface mesh. The resulting density profiles can be used for thickness measurements,
membrane analysis, and other downstream analyses.

Pipeline Citation: Barad BA*, Medina M*, Fuentes D, Wiseman RL, Grotjahn DA.
Quantifying organellar ultrastructure in cryo-electron tomography using a surface morphometrics pipeline.
J Cell Biol 2023.

Thickness Citation: Medina M, Chang Y-T, Rahmani H, Frank M, Khan Z, Fuentes D, Heberle FA, Waxham MN,
Barad BA, Grotjahn DA. Surface Morphometrics reveals local membrane thickness variation in organellar
subcompartments. J Cell Biol 2025.
"""

import numpy as np
import mrcfile
import scipy.interpolate as interp
from glob import glob
from pathlib import Path
from sys import argv
import yaml
import click
from graph_tool import load_graph


def load_graph_data(filename, voxsize):
    """
    Load xyz coordinates and normal vectors from a graph-tool .gt file.

    Parameters
    ----------
    filename : str
        Path to the .gt file
    voxsize : float
        The voxel size of the mrc data (for converting surface coords to voxel coords)

    Returns
    -------
    tuple
        xyz coordinates (3, n_vertices), n_v normal vectors (3, n_vertices), and graph object
    """
    graph = load_graph(filename)

    # Get xyz coordinates and convert from surface units (nm) to voxel units
    xyz = graph.vp.xyz.get_2d_array([0, 1, 2]) / voxsize

    # Get normal vectors and scale to convert nm steps to voxel steps
    # n_v is a unit vector; dividing by voxsize converts nm distance to voxel distance
    n_v = graph.vp.n_v.get_2d_array([0, 1, 2]) / voxsize

    return xyz, n_v, graph

def load_mrc(filename, angstroms=False):
    """
    Load mrc data from an mrc file using mrcfile.

    Parameters
    ----------
    filename : str
        Name of the mrc file
    angstroms : bool
        If True, keep angstrom units; if False, convert to nm

    Returns
    -------
    tuple
        data array, data_matrix for interpolation, voxel size, and origin
    """
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

def interpolate(data, data_matrix, xyz, n_v, sample_spacing=0.25, angstroms=False, scan_range=10):
    """
    Interpolate the values of the mrc data along each normal vector.

    Stepping from -scan_range to +scan_range in parameterized steps, interpolate
    the values of the mrc data along each normal vector using scipy.interpn.

    Parameters
    ----------
    data : ndarray
        The mrc data
    data_matrix : tuple
        Coordinate arrays for interpolation
    xyz : ndarray
        The xyz coordinates of the faces
    n_v : ndarray
        The normal vectors
    sample_spacing : float
        Distance in nm between samples (default: 0.25)
    angstroms : bool
        If True, scale samples to angstroms
    scan_range : float
        Half-range in nm to scan along normal vectors (default: 10)

    Returns
    -------
    ndarray
        Interpolated values array (n_triangles x nsamples)
    """
    averages = []
    # Calculate number of samples from spacing and range
    nsamples = int(2 * scan_range / sample_spacing) + 1
    # Create an array of nm steps from -scan_range to +scan_range
    samples = np.linspace(-scan_range, scan_range, nsamples)
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
        averages.append(value_array_temp[nsamples // 2])
        # Store the interpolated values for each normal vector
        # print(value_array_temp)
        value_array[i] = value_array_temp
    print(np.mean(averages))
    print(value_array.shape)
    print(xyz.shape)
    return value_array


def sample_density_single(mrc_file, graph_file, sample_spacing=0.25, scan_range=10, angstroms=False):
    """
    Sample density values for a single graph file.

    This is a reusable function that can be called by other scripts.

    Parameters
    ----------
    mrc_file : str
        Path to the tomogram MRC file
    graph_file : str
        Path to the graph-tool .gt file
    sample_spacing : float
        Distance in nm between samples (default: 0.25)
    scan_range : float
        Half-range in nm to scan along normal vectors (default: 10)
    angstroms : bool
        If True, use angstrom units (default: False)

    Returns
    -------
    tuple
        (value_array, x_positions, voxsize) where value_array is (n_triangles x nsamples),
        x_positions is the array of sample positions, and voxsize is the voxel size
    """
    # Load MRC data
    data, data_matrix, voxsize, origin = load_mrc(mrc_file, angstroms=angstroms)

    # Load graph data
    xyz, n_v, graph = load_graph_data(graph_file, voxsize)

    # Sample density along normals
    value_array = interpolate(data, data_matrix, xyz, n_v,
                              sample_spacing=sample_spacing, angstroms=angstroms,
                              scan_range=scan_range)

    # Generate x positions
    nsamples = int(2 * scan_range / sample_spacing) + 1
    x_positions = np.linspace(-scan_range, scan_range, nsamples)
    if angstroms:
        x_positions = x_positions * 10.

    return value_array, x_positions, voxsize


def sample_density_for_tomogram(filename, work_dir, angstroms=False, sample_spacing=0.25, scan_range=10, radius_hit=None):
    """
    Sample density values from a tomogram along surface normal vectors.

    Parameters
    ----------
    filename : str
        Path to the tomogram MRC file
    work_dir : str
        Working directory containing the graph-tool .gt files
    angstroms : bool
        If True, use angstrom units
    sample_spacing : float
        Distance in nm between samples (default: 0.25)
    scan_range : float
        Half-range in nm to scan along normal vectors
    radius_hit : int or None
        If specified, only process files with this radius_hit value
    """
    mrcbase = filename.split(".mrc")[0].split("/")[-1]
    print(f"Processing {mrcbase}")

    # Build glob pattern based on radius_hit - look for .gt files
    if radius_hit is not None:
        files = glob(work_dir + mrcbase + f"*.AVV_rh{radius_hit}.gt")
    else:
        files = glob(work_dir + mrcbase + f"*.AVV_rh*.gt")

    if not files:
        print(f"No graph files (.gt) found for {mrcbase}")
        return

    print(f"Found {len(files)} files to process")

    # Process each graph file
    for file in files:
        print(f"Processing {file}")
        value_array, positions, voxsize = sample_density_single(
            filename, file,
            sample_spacing=sample_spacing, scan_range=scan_range, angstroms=angstroms
        )
        # Save the interpolated values to a csv file (same basename as .gt file)
        header = ",".join([f"{p:.4f}" for p in positions])
        output_file = file[:-3] + "_sampling.csv"
        print(f"Saving to {output_file}")
        np.savetxt(output_file, value_array, delimiter=",", header=header, comments="")

@click.command()
@click.argument('configfile', type=click.Path(exists=True))
@click.argument('mrcfile', required=False, default=None)
@click.option('--sample_spacing', type=float, default=None, help='Distance in nm between samples (overrides config)')
@click.option('--scan_range', type=float, default=None, help='Half-range in nm to scan (overrides config)')
def sample_density_cli(configfile, mrcfile, sample_spacing, scan_range):
    """
    Sample tomogram density along surface normal vectors.

    For each triangle in the surface meshes, samples density values from the
    tomogram along the normal vector. Output can be used for thickness
    measurements and other downstream analyses.

    CONFIGFILE: Path to the config.yml file

    MRCFILE: Optional path to a specific tomogram MRC file to process.
             If not provided, processes all MRC files in tomo_dir.
    """
    run_sample_density(configfile, mrcfile, sample_spacing_override=sample_spacing,
                       scan_range_override=scan_range)


def run_sample_density(configfile, mrcfile=None, sample_spacing_override=None, scan_range_override=None):
    """
    Main function to run density sampling from config file.

    Parameters
    ----------
    configfile : str
        Path to config.yml
    mrcfile : str or None
        Optional specific tomogram MRC file to process
    sample_spacing_override : float or None
        Override sample_spacing from config
    scan_range_override : float or None
        Override scan_range from config
    """
    # Load config
    with open(configfile) as file:
        config = yaml.safe_load(file)

    # Validate and set up directories
    # tomo_dir contains the tomogram MRC files
    if not config.get("tomo_dir"):
        print("tomo_dir not specified in config.yml")
        return
    tomo_dir = config["tomo_dir"]
    if not tomo_dir.endswith("/"):
        tomo_dir += "/"

    # work_dir contains the curvature CSV files and will store sampling output
    if not config.get("work_dir"):
        print("work_dir not specified in config.yml")
        return
    work_dir = config["work_dir"]
    if not work_dir.endswith("/"):
        work_dir += "/"

    # Get density sampling settings from config (falls back to thickness_measurements for compatibility)
    density_config = config.get("density_sampling", config.get("thickness_measurements", {}))
    angstroms = config.get("surface_generation", {}).get("angstroms", False)
    radius_hit = config.get("curvature_measurements", {}).get("radius_hit", None)

    # Set sampling parameters (CLI overrides take precedence)
    sample_spacing = sample_spacing_override if sample_spacing_override is not None else density_config.get("sample_spacing", 0.25)
    scan_range = scan_range_override if scan_range_override is not None else density_config.get("scan_range", 10)
    nsamples = int(2 * scan_range / sample_spacing) + 1

    # Warn if sample_spacing doesn't divide evenly into scan_range
    if (scan_range % sample_spacing) != 0:
        print(f"WARNING: sample_spacing ({sample_spacing}) does not divide evenly into scan_range ({scan_range}).")
        print(f"         Actual spacing will be {2 * scan_range / (nsamples - 1):.4f} nm.")

    print(f"Density sampling settings:")
    print(f"  Tomogram directory: {tomo_dir}")
    print(f"  Work directory: {work_dir}")
    print(f"  Angstroms: {angstroms}")
    print(f"  Sample spacing: {sample_spacing} nm")
    print(f"  Scan range: {scan_range} nm")
    print(f"  N samples: {nsamples} (computed from spacing and range)")
    if radius_hit:
        print(f"  Radius hit: {radius_hit}")

    # Determine which MRC files to process
    if mrcfile:
        mrcs = [mrcfile]
    else:
        mrcs = glob(tomo_dir + "*.mrc")

    if not mrcs:
        print(f"No MRC files found in {tomo_dir}")
        return

    print(f"Found {len(mrcs)} MRC file(s) to process")

    # Process each MRC file
    for mrc in mrcs:
        sample_density_for_tomogram(mrc, work_dir, angstroms=angstroms, sample_spacing=sample_spacing,
                                    scan_range=scan_range, radius_hit=radius_hit)


if __name__ == "__main__":
    sample_density_cli()
