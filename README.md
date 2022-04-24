# Surface Morphometrics Pipeline
![Workflow Figure](https://raw.githubusercontent.com/GrotjahnLab/surface_morphometrics/master/Workflow_title.png)
### Quantification of Membrane Surfaces Segmented from Cryo-ET or other volumetric imaging.  
Author: __Benjamin Barad__/*<benjamin.barad@gmail.com>*. 

Developed in close collaboration with Michaela Medina

A pipeline of tools to generate robust open mesh surfaces from voxel segmentations of biological membranes
using the Screened Poisson algorithm, calculate morphological features including curvature and membrane-membrane distance
using pycurv's vector voting framework, and tools to convert these morphological quantities into morphometric insights.


## Installation:
1. Clone this git repository: `git clone https://github.com/grotjahnlab/surface_morphometrics.git`
2. Install the conda environment: `conda env create -f environment.yml`
3. Activate the conda environment: `source activate morphometrics`
4. Install additional dependencies: `pip install -r pip_requirements.txt`


## Running the configurable pipeline
1. Edit the `config.yml` file for your specific project needs.
2. Run the surface reconstruction for all segmentations: `python segmentation_to_meshes.py config.yml`
3. Run pycurv for each surface (recommended to run individually in parallel with a cluster):`python pycurv_analysis.py config.yml ${i}.surface.vtp`
4. Measure intra- and inter-surface distances and orientations (also best to run this one in parallel for each original segmentation): `python measure_distances.py config.yml ${i}.mrc`
5. Combine the results of the pycurv analysis into aggregate Experiments and generate statistics and plots. This requires some manual coding using the Experiment class and its associated methods in the `morphometrics_stats.py`. `mitochondria_statistics.py` is an example of how to do this.


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
    1. `morphometrics_stats.py` to generate graphs and statistics with pandas.
    2. [Paraview](https://www.paraview.org/) for 3D surface mapping of quantifications.

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

> **A surface morphometrics toolkit to quantify organellar membrane ultrastructure using cryo-electron tomography.**  
> Benjamin A. Barad<sup>†</sup>, Michaela Medina<sup>†</sup>, Daniel Fuentes, R. Luke Wiseman, Danielle A. Grotjahn  
> *bioRxiv* 2022.01.23.477440; doi: https://doi.org/10.1101/2022.01.23.477440

All scientific software is dependent on other libraries, but the surface morphometrics toolkit is particularly dependent on [PyCurv](https://github.com/kalemaria/pycurv), which provides the vector voted curvature measurements and the triangle graph framework. As such, please also cite the pycurv manuscript:

> **Reliable estimation of membrane curvature for cryo-electron tomography.**  
> Maria Salfer,Javier F. Collado,Wolfgang Baumeister,Rubén Fernández-Busnadiego,Antonio Martínez-Sánchez  
> *PLOS Comp Biol* August 2020; doi: https://doi.org/10.1371/journal.pcbi.1007962  

