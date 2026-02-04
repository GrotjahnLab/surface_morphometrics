# Surface Morphometrics Pipeline
![Workflow Figure](https://raw.githubusercontent.com/GrotjahnLab/surface_morphometrics/master/Workflow_title.png)
### Quantification of Membrane Surfaces Segmented from Cryo-ET or other volumetric imaging.  
Author: __Benjamin Barad__/*<benjamin.barad@gmail.com>*. 

Developed in close collaboration with Michaela Medina

A pipeline of tools to generate robust open mesh surfaces from voxel segmentations of biological membranes
using the Screened Poisson algorithm, calculate morphological features including curvature and membrane-membrane distance
using pycurv's vector voting framework, and tools to convert these morphological quantities into morphometric insights.

📖 **[Complete Quantifications Documentation](quantifications_documentation.md)** - Reference guide for all morphological measurements and their interpretations.


## Installation Options:
## Option 1: Native Installation through Conda
_This is the fastest and easiest starting point for most linux boxes, and now for mac as well_
1. Clone this git repository: `git clone https://github.com/grotjahnlab/surface_morphometrics.git`
2. Install the conda environment: `conda env create -f environment.yml`
3. Activate the conda environment: `conda activate morphometrics`

Note: Older ubuntu installs (and probably some other linux distributions!) have some known issues with graph-tool. If you run into issues with the main environment file, try using `environment-ubuntu.yml` instead. This worked on my Ubuntu 22.04 LTS box.



## Option 1: Docker Container
_This allows operation on windows and other operating systems where graph-tool or pymeshlab do not play nice with each other or with conda._
1. Clone this git repository: 
```bash
git clone https://github.com/grotjahnlab/surface_morphometrics.git
cd surface_morphometrics
```
2. Start the containerized environment:
```bash
cd docker
./sm-up.sh
```
This will pull and start the Docker container with all dependencies pre-installed.
3. When finished, exit the container and stop it:
```bash
exit # Exit the container
./sm-down.sh # Stop the container
```


## Working with Docker Container

1. Starting a Session:
```bash
cd surface_morphometrics # Go to project directory
cd docker                # Enter docker directory
./sm-up.sh               # Start and enter container
```
2. Inside Container:
- The environment is pre-configured
- All dependencies are installed
- You can directly run the pipeline commands
3. Ending a Session:
```bash
exit           # Exit the container
./sm-down.sh   # Stop the container
```

## Example data

There is tutorial data available in the `example_data` folder. Uncompress the tar file with:
```bash
cd example_data
tar -xzvf examples.tar.gz
```

There are two example datasets: `TE1.mrc` and `TF1.mrc`.
You can open them with `mrcfile`, like so:

```python
import mrcfile

with mrcfile.open('TE1.mrc', permissive=True) as mrc:
    print(mrc.data.shape)  # TE1.mrc has shape (312, 928, 960)
```

In some cases the file header may be non-standard (for example, mrc files exported from Amira software). In these cases, the `permissive=True` keyword argument is required, and you can ignore the warning that the file may be corrupt. All the surface morphometrics toolkit scripts will still run correctly.

## Running the configurable pipeline

Running the full pipeline on a 4 core laptop with the tutorial datasets takes about 8 hours (3 
for TE1, 5 for TF1), mostly in steps 3 and 4. With cluster parallelization, the full pipeline 
can run in 2 hours for as many tomograms as desired.

1. Edit the `config.yml` file for your specific project needs. **New users:** we recommend starting with `isotropic_remesh: true` and `simplify: false` in the `surface_generation` section for higher quality meshes with near-equilateral triangles. A `target_area` between 1.0 and 3.0 nm^2 generally yields good results, but smaller triangles significantly increase computation time. The default config uses `simplify: true` for faster processing.
2. Run the surface reconstruction for all segmentations: `python segmentation_to_meshes.py config.yml`
3. Run pycurv for each surface (recommended to run individually in parallel with a cluster): `python 
run_pycurv.py config.yml ${i}.surface.vtp`

    You may see warnings aobut the curvature, this is normal and you do not need to worry.

4. Measure intra- and inter-surface distances and orientations (also best to run this one in parallel for each original segmentation): `python measure_distances_orientations.py config.yml ${i}.mrc`
5. Combine the results of the pycurv analysis into aggregate Experiments and generate statistics and plots. This requires some manual coding using the Experiment class and its associated methods in the `morphometrics_stats.py`. Everything is roughly organized around working with the CSVs in pandas dataframes. Running  `morphometrics_stats.py` as a script with the config file and a filename will output a pickle file with an assembled "experiment" object for all the tomos in the data folder. Reusing a pickle file will make your life way easier if you have dozens of tomograms to work with, but it doesn't save too much time with just the example data...

### Examples of generating statistics and plots:
* `python single_file_histogram.py filename.csv -n feature` will generate an area-weighted histogram for a feature of interest in a single tomogram. I am using a variant of this script to respond to reviews asking for more per-tomogram visualizations!
* `python single_file_2d.py filename.csv -n1 feature1 -n2 feature2` will generate a 2D histogram for 2 features of interest for a single surface.
* `mitochondria_statistics.py` shows analysis and comparison of multiple experiment objects for different sets of tomograms (grouped by treatment in this case). Every single plot and statistic in the preprint version of the paper gets generated by this script.


## Running individual steps without pipelining
Individual steps are available as click commands in the terminal, and as functions

1. Robust Mesh Generation
    1. `mrc2xyz.py` to prepare point clouds from voxel segmentation
    2. `xyz2ply.py` to perform screened poisson reconstruction and mask the surface
    3. `ply2vtp.py` to convert ply files to vtp files ready for pycurv
2. Surface Morphology Extraction
    1. `curvature.py` to run pycurv in an organized way on pregenerated surfaces
    2. `intradistance_verticality.py` to generate distance metrics and verticality measurements within a surface.
    3. `interdistance_orientation.py` to generate distance metrics and orientation measurements between surfaces.
    4. Outputs: gt graphs for further analysis, vtp files for paraview visualization, and CSV files for         pandas-based plotting and statistics
3. Morphometric Quantification - there is no click function for this, as the questions answered depend on the biological system of interest!
    1. `morphometrics_stats.py` is a set of classes and functions to generate graphs and statistics with pandas.
    2. [Paraview](https://www.paraview.org/) for 3D surface mapping of quantifications.
    3. **[Quantifications Documentation](quantifications_documentation.md)** - Complete reference for all morphological measurements and their interpretations.

## File Descriptions:
* Files with.xyz extension are point clouds converted, in nm or angstrom scale. This is a flat text file with `X Y Z` coordinates in each line.
* Files with .ply extension are the surface meshes (in a binary format), which will be scaled in nm or angstrom scale, and work in many different softwares, including [Meshlab](https://www.meshlab.net/). 
* Files with .vtp extension are the same surface meshes in the [VTK](https://vtk.org/) format.
        * The .surface.vtp files are a less cross-compatible format, so you can't use them with as many types of software, but they are able to store all the fun quantifications you'll do!. [Paraview](https://www.paraview.org/) or [pyvista](https://docs.pyvista.org/) can load this format. This is the format pycurv reads to build graphs.
        * The .AVV_rh8.vtp files are those output from downstream components of the pipeline, and generally have the most available visualizations in paraview and pyvista. 
* Files with .gt extension are triangle graph files using the `graph-tool` python toolkit. These graphs enable rapid neighbor-wise operations such as tensor voting, but are not especially useful for manual inspection.
* Files with .csv extension are quantification outputs per-triangle. These are the files you'll use to generate statistics and plots.
* Files with .log extension are log files, mostly from the output of the pycurv run.
* Quantifications (plots and statistical tests) are output in csv, svg, and png formats. 

## Troubleshooting
1. Warnings of the type `Gaussian or Mean curvature of X has a large computation error`... can be ignored, as they get cleaned up by pycurv. These warnings are now suppressed by default.
2. MRC files that are output by AMIRA don't have proper machine stamps by default. They need to be imported with `mrcfile.open(filename, permissive=True)`. This is also true for many other softwares, including Dragonfly.
3. Pycurv has recently undergone significant performance improvements and has more feedback; if it seems to be hanging indefinitely, try setting cores to 1 in the config file.

## Dependencies
1. Numpy
2. Scipy
3. Pandas
4. mrcfile
5. Click
6. Matplotlib
7. Pymeshlab
8. Pycurv   
    1. Pyto
    2. Graph-tool


## Citation
The development of this toolkit and examples of useful applications can be found in the following manuscript. Please cite it if you use this software in your research, or extend it to make improvements!

> **Quantifying organellar ultrastructure in cryo-electron tomography using a surface morphometrics pipeline.**
> Benjamin A. Barad<sup>†</sup>, Michaela Medina<sup>†</sup>, Daniel Fuentes, R. Luke Wiseman, Danielle A. Grotjahn
> *Journal of Cell Biology* 2023, 222(4), e202204093; doi: https://doi.org/10.1083/jcb.202204093

All scientific software is dependent on other libraries, but the surface morphometrics toolkit is particularly dependent on [PyCurv](https://github.com/kalemaria/pycurv), which provides the vector voted curvature measurements and the triangle graph framework. As such, please also cite the pycurv manuscript:

> **Reliable estimation of membrane curvature for cryo-electron tomography.**  
> Maria Salfer,Javier F. Collado,Wolfgang Baumeister,Rubén Fernández-Busnadiego,Antonio Martínez-Sánchez  
> *PLOS Comp Biol* August 2020; doi: https://doi.org/10.1371/journal.pcbi.1007962  

