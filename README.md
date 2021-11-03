# Surface Morphometrics Toolkit
### Quantification of Membrane Surfaces Segmented from Cryo-ET or other volumetric imaging.  
Author: __Benjamin Barad__/*<benjamin.barad@gmail.com>*. 

Developed in close collaboration with Michaela Medina

A collection of tools to generate robust open mesh surfaces from voxel segmentations of biological membranes
using the Screened Poisson algorithm, calculate morphological features including curvature and membrane-membrane distance
using pycurv's vector voting framework, and tools to convert these morphological quantities into morphometric insights.


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


## Organization
Scripts are still in a rough state tailored to mitochondrial workflow. They will be generalized in the coming weeks.
1. Robust Mesh Generation
    1. `mrc2pts.py` to prepare point clouds
    2. `make_mesh.py` to perform screened poisson reconstruction and associated processes
    3. `mask_and_convert_ply.py` to convert ply files to vtp files ready for pycurv
2. Surface Morphology Extraction
    1. `curvature_calculation_args.py` to run pycurv in an organized way on pregenerated surfaces
    2. `surface_to_surface.py` to generate distance metrics between and within membranes on pycurv-processed graphs 
    3. Outputs: gt graphs for further analysis, vtp files for paraview visualization, and CSV files for pandas-based plotting and statistics
3. Morphometric Quantification
    1. `csv_quantifications.py` to generate graphs and statistics with pandas.
    2. Paraview for structured mapping.

