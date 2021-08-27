import sys
import time
from os.path import isfile
from graph_tool import load_graph
import gzip
from os import remove
import pandas as pd
import numpy as np
from scipy import ndimage
import os
from pathlib import Path
import pathos.pools as pp
from functools import partial

from pycurv import (
    pexceptions, normals_directions_and_curvature_estimation, run_gen_surface,
    TriangleGraph, PointGraph, curvature_estimation, merge_vtp_files,
    split_segmentation, MAX_DIST_SURF, THRESH_SIGMA1, rescale_surface)
from pycurv import pycurv_io as io

"""
A script with example applications of the PyCurv package for estimation of
membrane curvature.

Author: Maria Salfer (Max Planck Institute for Biochemistry)
"""

__author__ = 'Maria Salfer'


def convert_vtp_to_stl_surface_and_mrc_curvatures(
        surf_vtp_file, outfile_base, scale, size):
    """
    Converts the '.vtp' surface file to '.stl' file and converts selected
    vtkPolyData cell arrays from the '.vtp' file as 3-D volumes in '.mrc' files.
    The selected arrays are: "kappa_1", "kappa_2", "curvedness_VV".

    Args:
        surf_vtp_file (str): surface .vtp file, should contain the final surface
            with curvatures
        outfile_base (str): base name for the output .mrc and .log files
        scale (tuple): pixel size (X, Y, Z) of the membrane mask in units of the
            surface
        size (tuple): size (X, Y, Z) of the membrane mask

    Returns:
        None
    """
    # Scaling the '.vtp' surface file back to voxels:
    surf_voxels_vtp_file = (surf_vtp_file[0:-4] + '_voxels.vtp')
    if not isfile(surf_voxels_vtp_file):
        poly_nm = io.load_poly(surf_vtp_file)
        reverse_scale = (1/scale[0], 1/scale[1], 1/scale[2])
        poly_voxels = rescale_surface(poly_nm, reverse_scale)
        io.save_vtp(poly_voxels, surf_voxels_vtp_file)
        print("The '.vtp' file was scaled back to voxels: {}".format(
            surf_voxels_vtp_file))

    # Converting the '.vtp' surface voxels file to '.stl' file:
    surf_voxels_stl_file = (surf_voxels_vtp_file[0:-4] + '.stl')
    if not isfile(surf_voxels_stl_file):
        io.vtp_file_to_stl_file(surf_voxels_vtp_file, surf_voxels_stl_file)
        print("The '.vtp' file {} was converted to .stl format".format(
            surf_voxels_vtp_file))

    # Converting vtkPolyData selected cell arrays from the '.vtp' file as 3-D
    # volumes in '.mrc' files (and saving them as '.mrc.gz' files).
    # max voxel value & .log files:
    _vtp_arrays_to_mrc_volumes(
        surf_voxels_vtp_file, outfile_base, size, log_files=True)
    # mean voxel value & no .log files:
    _vtp_arrays_to_mrc_volumes(
        surf_voxels_vtp_file, outfile_base, size, mean=True)


def _vtp_arrays_to_mrc_volumes(
        surf_vtp_file, outfile_base, size, mean=False, log_files=False,
        compress=False):
    """
    This function converts selected vtkPolyData cell arrays from the '.vtp' file
    as 3-D volumes in '.mrc' files. The selected arrays are: "kappa_1",
    "kappa_2", "curvedness_VV".

    Args:
        surf_vtp_file (str): surface .vtp file, should contain the final surface
            with curvatures
        outfile_base (str): base name for the output .mrc (and .log) files
        size (tuple): size (X, Y, Z) of the membrane mask
        mean (boolean, optional): if True (default False), in case multiple
            triangles map to the same voxel, takes the mean value, else the
            maximal value
        log_files (boolean, optional): if True (default False), writes the log
            files for such cases
        compress (boolean, optional): if True (default False), compresses the
            '.mrc' files with gzip.

    Returns:
        None
    """
    array_names = ["kappa_1", "kappa_2", "curvedness_VV"]
    names = ["max_curvature", "min_curvature", "curvedness"]

    if mean:
        voxel_value_str = "voxel_mean"
    else:
        voxel_value_str = "voxel_max"
    mrcfilenames = []
    logfilenames = []
    for name in names:
        mrcfilename = "{}.{}.{}.mrc".format(outfile_base, name, voxel_value_str)
        mrcfilenames.append(mrcfilename)
        if log_files:
            logfilename = "{}.{}.{}.log".format(
                outfile_base, name, voxel_value_str)
        else:
            logfilename = None
        logfilenames.append(logfilename)

    # Load the vtkPolyData object from the '.vtp' file, calculate the volumes
    # from arrays, write '.log' files, and save the volumes as '.mrc' files:
    poly = io.load_poly(surf_vtp_file)
    for i, array_name in enumerate(array_names):
        volume = io.poly_array_to_volume(
            poly, array_name, size, logfilename=logfilenames[i], mean=mean)
        io.save_numpy(volume, mrcfilenames[i])

    if compress:
        # Gunzip the '.mrc' files and delete the uncompressed files:
        for mrcfilename in mrcfilenames:
            with open(mrcfilename) as f_in, \
                    gzip.open(mrcfilename + '.gz', 'wb') as f_out:
                f_out.writelines(f_in)
            remove(mrcfilename)
            print('Archive {}.gz was written'.format(mrcfilename))


def new_workflow(
        base_filename, seg_file, fold, pixel_size, radius_hit,
        methods=['VV'], page_curvature_formula=False, area2=True,
        label=1, filled_label=None, unfilled_mask=None, holes=0,
        remove_wrong_borders=True, min_component=100, only_normals=False,
        cores=6, runtimes=''):
    """
    A script for running all processing steps to estimate membrane curvature.

    The three steps are: 1. signed surface generation, 2. surface cleaning using
    a graph, 3. curvature calculation using a graph generated from the clean
    surface.

    It was written for Javier's data. Segmentation is not split into regions.
    Second pass, consisting of normals and curvature calculations, can run in
    parallel on multiple cores (for RVV and AVV, but not for SSVV).

    Args:
        base_filename (str): base file name for saving the output files
        seg_file (str): membrane segmentation mask, pass '' if surface exists
        fold (str): path where the input membrane segmentation is and where the
            output will be written
        pixel_size (float): pixel size in nanometer of the segmentation
        radius_hit (float): radius in length unit of the graph, e.g. nanometers;
            it should be chosen to correspond to radius of smallest features of
            interest on the surface
        methods (list, optional): all methods to run in the second pass ('VV'
            and 'SSVV' are possible, default is 'VV')
        page_curvature_formula (boolean, optional): if True (default False),
            normal curvature formula from Page et al. is used in VV (see
            collect_curvature_votes)
        area2 (boolean, optional): if True (default), votes are weighted by
            triangle area also in the second step (principle directions and
            curvatures estimation)
        label (int, optional): label to be considered in the membrane mask
            (default 1)
        filled_label (int, optional): if the membrane mask was filled with this
            label (default None), complementing it to a compartment
            segmentation, a better surface generation will be used (with a
            slight smoothing; holes are closed automatically by the filling)
        unfilled_mask (numpy.ndarray, optional): if given (default None), apply
            this mask on the extracted surface using a membrane segmentation,
            instead of the segmentation itself; not used if filled_label is
            given
        holes (int, optional): if > 0, small holes in the membrane segmentation
            are closed with a cube of that size in pixels before curvature
            estimation (default 0); not used if filled_label is given
        remove_wrong_borders (boolean, optional): if True (default), wrong
            artefact surface borders will be removed
        min_component (int, optional): if > 0 (default 100), small
            disconnected surface components having triangles within this number
            will be removed
        only_normals (boolean, optional): if True (default False), only normals
            are estimated, without principal directions and curvatures, only the
            graph with the orientations class, normals or tangents is returned.
        cores (int, optional): number of cores to run VV in parallel (default 6)
        runtimes (str, optional): if given, runtimes and some parameters are
            added to this file (default '')

    Returns:
        None
    """
    log_file = '{}{}.{}_rh{}.log'.format(
                fold, base_filename, methods[0], radius_hit)
    sys.stdout = open(log_file, 'w')

    t_begin = time.time()

    surf_file = base_filename + ".surface.vtp"
    if not isfile(fold + surf_file):
        if seg_file == '' or not isfile(fold + seg_file):
            raise pexceptions.PySegInputError(
                expr="new_workflow",
                msg="The segmentation file not given or not found")

        seg = io.load_tomo(fold + seg_file)
        assert(isinstance(seg, np.ndarray))
        data_type = seg.dtype

        if filled_label is not None:  # if lumen segmentation given:
            # Surface generation with compartment segmentation using Marching
            # Cubes algorithm and applying the mask of membrane segmentation.
            print("\nMaking membrane and compartment binary segmentations...")
            binary_seg = (seg == label).astype(data_type)
            if not np.any(binary_seg):
                raise pexceptions.PySegInputError(
                    expr="new_workflow",
                    msg="Label not found in the segmentation!")
            # Combine the membrane and lumen segmentations into the compartment
            # (filled) segmentation:
            filled_binary_seg = np.logical_or(
                seg == label, seg == filled_label).astype(data_type)
            print("\nGenerating a surface...")
            surf = run_gen_surface(
                filled_binary_seg, fold + base_filename, lbl=1,
                other_mask=binary_seg, isosurface=True, sg=1, thr=THRESH_SIGMA1)
            # Write the resulting binary segmentations into a file:
            filled_binary_seg_file = "{}{}.filled_binary_seg.mrc".format(
                fold, base_filename)
            io.save_numpy(filled_binary_seg, filled_binary_seg_file)
            binary_seg_file = "{}{}.binary_seg.mrc".format(fold, base_filename)
            io.save_numpy(binary_seg, binary_seg_file)

        else:  # Surface generation with Hoppe's algorithm and applying the mask
            # of membrane segmentation.
            print("Making the segmentation binary...")
            binary_seg = (seg == label).astype(data_type)
            if not np.any(binary_seg):
                raise pexceptions.PySegInputError(
                    expr="new_workflow",
                    msg="Label not found in the segmentation!")
            if holes > 0:  # close (reduce) holes in the segmentation
                cube_size = abs(holes)
                cube = np.ones((cube_size, cube_size, cube_size))
                print("\nReducing holes in the segmentation...")
                binary_seg = ndimage.binary_closing(
                    binary_seg, structure=cube, iterations=1).astype(data_type)
            # Write the resulting binary segmentation into a file:
            binary_seg_file = "{}{}.binary_seg.mrc".format(
                fold, base_filename)
            # io.save_numpy(binary_seg, binary_seg_file)
            print("\nGenerating a surface from the binary segmentation...")
            surf = run_gen_surface(binary_seg, fold + base_filename, lbl=1,
                                   other_mask=unfilled_mask)
    else:
        print('\nReading in the surface from file...')
        surf = io.load_poly(fold + surf_file)

    clean_graph_file = '{}.scaled_cleaned.gt'.format(base_filename)
    clean_surf_file = '{}.scaled_cleaned.vtp'.format(base_filename)
    if not isfile(fold + clean_graph_file) or not isfile(fold + clean_surf_file):
        print('\nBuilding a triangle graph from the surface...')
        tg = TriangleGraph()
        scale = (pixel_size, pixel_size, pixel_size)
        tg.build_graph_from_vtk_surface(surf, scale)
        if tg.graph.num_vertices() == 0:
            raise pexceptions.PySegInputError(
                expr="new_workflow", msg="The graph is empty!")
        print('The graph has {} vertices and {} edges'.format(
            tg.graph.num_vertices(), tg.graph.num_edges()))

        # Remove the wrong borders (surface generation artefact)
        if remove_wrong_borders:
            b = MAX_DIST_SURF  # "padding" from masking in surface generation
            print('\nFinding triangles that are {} pixels to surface borders...'
                  .format(b))
            tg.find_vertices_near_border(b * pixel_size, purge=True)
            print('The graph has {} vertices and {} edges'.format(
                tg.graph.num_vertices(), tg.graph.num_edges()))

        # Filter out possibly occurring small disconnected fragments
        if min_component > 0:
            print('\nFinding small connected components of the graph...')
            tg.find_small_connected_components(
                threshold=min_component, purge=True, verbose=True)
            print('The graph has {} vertices and {} edges'.format(
                tg.graph.num_vertices(), tg.graph.num_edges()))

        # Saving the scaled (and cleaned) graph and surface:
        tg.graph.save(fold + clean_graph_file)
        surf_clean = tg.graph_to_triangle_poly()
        io.save_vtp(surf_clean, fold + clean_surf_file)
    else:
        print('\nReading in the cleaned graph and surface from files...')
        surf_clean = io.load_poly(fold + clean_surf_file)
        tg = TriangleGraph()
        tg.graph = load_graph(fold + clean_graph_file)

    t_end = time.time()
    duration = t_end - t_begin
    minutes, seconds = divmod(duration, 60)
    print('Surface and graph generation (and cleaning) took: {} min {} s'
          .format(minutes, seconds))

    # Running the modified Normal Vector Voting algorithms:
    gt_file1 = '{}{}.NVV_rh{}.gt'.format(fold, base_filename, radius_hit)
    method_tg_surf_dict = {}
    if not isfile(gt_file1):
        if runtimes != '':
            with open(runtimes, 'w') as f:
                f.write("num_v;radius_hit;g_max;avg_num_neighbors;cores;"
                        "duration1;method;duration2\n")
        method_tg_surf_dict = normals_directions_and_curvature_estimation(
            tg, radius_hit, methods=methods, full_dist_map=False,
            graph_file=gt_file1, page_curvature_formula=page_curvature_formula,
            area2=area2, only_normals=only_normals, poly_surf=surf_clean,
            cores=cores, runtimes=runtimes)
    elif only_normals is False:
        if runtimes != '':
            with open(runtimes, 'w') as f:
                f.write("method;duration2\n")
        for method in methods:
            tg_curv, surface_curv = curvature_estimation(
                radius_hit, graph_file=gt_file1, method=method,
                page_curvature_formula=page_curvature_formula, area2=area2,
                poly_surf=surf_clean, cores=cores, runtimes=runtimes)
            method_tg_surf_dict[method] = (tg_curv, surface_curv)

    if only_normals is False:  # Saving the output (graph and surface objects)
        # for later filtering or inspection in ParaView:
        for method in list(method_tg_surf_dict.keys()):
            (tg, surf) = method_tg_surf_dict[method]
            if method == 'VV':
                if page_curvature_formula and (area2 is False):
                    method = 'NVV'
                elif page_curvature_formula is False:
                    if area2 is False:
                        method = 'RVV'
                    else:
                        method = 'AVV'

            gt_file = '{}{}.{}_rh{}.gt'.format(
                fold, base_filename, method, radius_hit)
            tg.graph.save(gt_file)
            surf_file = '{}{}.{}_rh{}.vtp'.format(
                fold, base_filename, method, radius_hit)
            io.save_vtp(surf, surf_file)


def calculate_PM_curvatures(fold, base_filename, radius_hit, cores=6):
    """
    Calculates plasma membrane curvatures with AVV using a pre-calculated
    estimated normals file.

    Args:
        fold (str): path where the input membrane segmentation is and where the
            output will be written
        base_filename (str): base file name for saving the output files
        radius_hit (float): radius in length unit of the graph, e.g. nanometers;
            it should be chosen to correspond to radius of smallest features of
            interest on the surface
        cores (int, optional): number of cores to run VV in parallel (default 6)

    Returns:
        None
    """
    gt_file_normals = "{}{}.NVV_rh{}.gt".format(fold, base_filename, radius_hit)
    tg = TriangleGraph()
    tg.graph = load_graph(gt_file_normals)

    tg_curv, surf_curv = curvature_estimation(
        radius_hit, graph_file=gt_file_normals, method='VV', cores=cores, sg=tg)

    gt_file_curv = "{}{}.AVV_rh{}.gt".format(fold, base_filename, radius_hit)
    tg_curv.graph.save(gt_file_curv)
    surf_file_curv = "{}{}.AVV_rh{}.vtp".format(fold, base_filename, radius_hit)
    io.save_vtp(surf_curv, surf_file_curv)


def extract_curvatures_after_new_workflow(
        fold, base_filename, radius_hit, methods=['VV'],
        page_curvature_formula=False, area2=True,
        exclude_borders=0, categorize_shape_index=False, regions=1):
    """
    Extracts curvature information from a .gt file generated by new_workflow
    into a .csv file. Optionally, values near surface borders can be excluded
    and shape index can be categorized.

    Args:
        fold (str): path where the input membrane segmentation is and where the
            output will be written
        base_filename (str): base file name for saving the output files
        radius_hit (float): radius in length unit of the graph, here nanometers;
            it should be chosen to correspond to radius of smallest features of
            interest on the surface
        methods (list, optional): all methods to run in the second pass ('VV'
            and 'SSVV' are possible, default is 'VV')
        page_curvature_formula (boolean, optional): if True (default False),
            normal curvature formula from Page et al. is used in VV (see
            collect_curvature_votes)
        area2 (boolean, optional): if True (default), votes are weighted by
            triangle area also in the second step (principle directions and
            curvatures estimation)
        exclude_borders (int, optional): if > 0, triangles within this distance
            from borders in dist and corresponding values will be excluded from
            the output files (graph .gt, surface.vtp file and .csv), iteratively
            starting from 0 until maximally this distance (integer by integer)
        categorize_shape_index (boolean, optional): if True (default False),
            shape index categories will be added to the input graph .gt and
            surface .vtp files as well as the output .csv file
        regions (int, optional): if > 1, extracts from all region files
            (numerated from 1 until this number before the extension) to one CSV
            without creating VTP and GT files, if exclude_borders > 0.

    Returns:
        None
    """
    log_file = '{}{}.{}_rh{}.log'.format(
                fold, base_filename, methods[0], radius_hit)
    sys.stdout = open(log_file, 'a')

    for method in methods:
        if method == 'VV':
            if page_curvature_formula and (area2 is False):
                method = 'NVV'
            elif page_curvature_formula is False:
                if area2 is False:
                    method = 'RVV'
                else:
                    method = 'AVV'
        print("Method: {}".format(method))
        # input graph and surface files
        gt_infile = '{}{}.{}_rh{}.gt'.format(
            fold, base_filename, method, radius_hit)
        vtp_infile = '{}{}.{}_rh{}.vtp'.format(
            fold, base_filename, method, radius_hit)
        # output csv, gt and vtp files (without excluding borders)
        csv_outfile = '{}{}.{}_rh{}.csv'.format(
            fold, base_filename, method, radius_hit)
        if categorize_shape_index:  # overwrite the input files
            gt_outfile = gt_infile
            vtp_outfile = vtp_infile
        else:
            gt_outfile = None
            vtp_outfile = None
        for dist in range(exclude_borders + 1):
            print("\nExtracting curvatures without {} nm from border".format(
                dist))
            if dist > 0:
                eb = "_excluding{}borders".format(dist)
                csv_outfile = '{}{}.{}_rh{}{}.csv'.format(
                    fold, base_filename, method, radius_hit, eb)
                if regions == 1:  # not for multiple regions
                    gt_outfile = '{}{}.{}_rh{}{}.gt'.format(
                        fold, base_filename, method, radius_hit, eb)
                    vtp_outfile = '{}{}.{}_rh{}{}.vtp'.format(
                        fold, base_filename, method, radius_hit, eb)

            if regions == 1:  # normal case
                # Create TriangleGraph object and load the graph file
                tg = TriangleGraph()
                tg.graph = load_graph(gt_infile)

                _extract_curvatures_from_graph(
                    tg, csv_outfile, dist, gt_outfile, vtp_outfile,
                    categorize_shape_index=categorize_shape_index)
            else:  # if multiple regions
                csv_region_outfiles = []
                for i in range(1, regions + 1):
                    print("\nRegion {}:".format(i))
                    # correct in/output files for multiple regions
                    gt_region_infile = gt_infile.replace(
                        base_filename, base_filename+str(i))
                    csv_region_outfile = csv_outfile.replace(
                        base_filename, base_filename+str(i))
                    csv_region_outfiles.append(csv_region_outfile)
                    gt_region_outfile = gt_outfile.replace(
                        base_filename, base_filename + str(i))
                    vtp_region_outfile = vtp_outfile.replace(
                        base_filename, base_filename + str(i))

                    # Create TriangleGraph object and load the graph file
                    tg = TriangleGraph()
                    tg.graph = load_graph(gt_region_infile)

                    _extract_curvatures_from_graph(
                        tg, csv_region_outfile, dist,
                        gt_region_outfile, vtp_region_outfile,
                        categorize_shape_index=categorize_shape_index, region=i)

                # join the region CSV files to one
                combined_df = pd.concat(
                    [pd.read_csv(f) for f in csv_region_outfiles])
                combined_df.to_csv(csv_outfile, index=False)
                # remove the region CSVs
                # for f in csv_region_outfiles:
                #     os.remove(f)


def _extract_curvatures_from_graph(
        sg, csv_file, exclude_borders, gt_file, vtp_file,
        categorize_shape_index, region=0):
    # If don't want to include curvatures near borders, filter out those
    if exclude_borders > 0 and sg.__class__.__name__ == "TriangleGraph":
        sg.find_vertices_near_border(exclude_borders, purge=True)

    # List of shape class labels of all vertices for the csv file:
    shape_index_class = []
    if categorize_shape_index:
        # Add a new property: categorical shape index (one value for class)
        sg.graph.vp.shape_index_cat = sg.graph.new_vertex_property("float")
        for v in sg.graph.vertices():
            si_v = sg.graph.vp.shape_index_VV[v]
            si_cat_v, si_class_v = _shape_index_classifier(si_v)
            sg.graph.vp.shape_index_cat[v] = si_cat_v
            shape_index_class.append(si_class_v)

    # Saving the changes into graph and surface files, if specified:
    if gt_file is not None:
        sg.graph.save(gt_file)
    if vtp_file is not None:
        # Transforming the resulting graph to a surface with triangles:
        surf = sg.graph_to_triangle_poly()
        io.save_vtp(surf, vtp_file)

    # Getting estimated principal curvatures from the output graph:
    kappa_1 = sg.get_vertex_property_array("kappa_1")
    kappa_2 = sg.get_vertex_property_array("kappa_2")
    gauss_curvature = sg.get_vertex_property_array("gauss_curvature_VV")
    mean_curvature = sg.get_vertex_property_array("mean_curvature_VV")
    shape_index = sg.get_vertex_property_array("shape_index_VV")
    curvedness = sg.get_vertex_property_array("curvedness_VV")

    # Writing all the curvature values and errors into a csv file:
    df = pd.DataFrame()
    if region > 0:  # add a column with region number
        df["region"] = [region for i in range(len(kappa_1))]
    df["kappa1"] = kappa_1
    df["kappa2"] = kappa_2
    df["gauss_curvature"] = gauss_curvature
    df["mean_curvature"] = mean_curvature
    df["shape_index"] = shape_index
    if categorize_shape_index:  # want the shape class labels
        df["shape_index_class"] = shape_index_class
    df["curvedness"] = curvedness
    if sg.__class__.__name__ == "TriangleGraph":
        triangle_areas = sg.get_vertex_property_array("area")
        df["triangleAreas"] = triangle_areas
    df.to_csv(csv_file, sep=';')


def _shape_index_classifier(x):
    """
    Maps shape index value to the representative (middle) value of each shape
    class and the class label.

    Args:
        x (float): shape index value, should be in range [-1, 1]

    Returns:
        A tuple of the representative (middle) value of each shape class and
        the class label, e.g. 0, 'Saddle' for values in range [-1/8, +1/8)
    """
    if x < -1:
        return None, None
    elif -1 <= x < -7/8.0:
        return -1, 'Spherical cup'
    elif -7/8.0 <= x < -5/8.0:
        return -0.75, 'Trough'
    elif -5/8.0 <= x < -3/8.0:
        return -0.5, 'Rut'
    elif -3/8.0 <= x < -1/8.0:
        return -0.25, 'Saddle rut'
    elif -1/8.0 <= x < 1/8.0:
        return 0, 'Saddle'
    elif 1/8.0 <= x < 3/8.0:
        return 0.25, 'Saddle ridge'
    elif 3/8.0 <= x < 5/8.0:
        return 0.5, 'Ridge'
    elif 5/8.0 <= x < 7/8.0:
        return 0.75, 'Dome'
    elif 7/8.0 <= x <= 1:
        return 1, 'Spherical cap'
    else:  # x > 1
        return None, None


def from_ply_workflow(
        ply_file, radius_hit, scale=(1, 1, 1), page_curvature_formula=False,
        methods=["VV"], area2=True, cores=6):
    """
    Estimates curvature for each triangle in a triangle mesh in PLY format.

    Args:
        mask_file (str): MRC file for
        radius_hit (float): radius in length unit of the graph, e.g. nanometers;
            it should be chosen to correspond to radius of smallest features of
            interest on the surface
        scale (tuple, optional): pixel size (X, Y, Z) in given units for
            scaling the surface if it is not scaled (default (1, 1, 1))
        page_curvature_formula (boolean, optional): if True (default False),
            normal curvature formula from Page et al. is used in VV (see
            collect_curvature_votes)
        methods (list, optional): all methods to run in the second pass ('VV'
            and 'SSVV' are possible, default is 'VV')
        area2 (boolean, optional): if True (default), votes are weighted by
            triangle area also in the second step (principle directions and
            curvatures estimation)
        cores (int, optional): number of cores to run VV in parallel (default 6)

    Returns:
        None
    """
    base_filename = os.path.splitext(ply_file)[0]
    log_file = '{}.{}_rh{}.log'.format(
        base_filename, methods[0], radius_hit)
    sys.stdout = open(log_file, 'a')

    # Transforming PLY to VTP surface format
    surf_file = base_filename + ".vtp"
    io.ply_file_to_vtp_file(ply_file, surf_file)

    # Reading in the surface and transforming it into a triangle graph
    print('\nReading in the surface file to get a vtkPolyData surface...')
    surf = io.load_poly(surf_file)
    if surf.GetNumberOfCells() == 0:
        print('The surface is empty, exiting.')
        return None
    print('\nBuilding a triangle graph from the surface...')
    tg = TriangleGraph()
    tg.build_graph_from_vtk_surface(surf, scale)
    if tg.graph.num_vertices() == 0:
        raise pexceptions.PySegInputError(
            expr="new_workflow", msg="The graph is empty!")
    print('The graph has {} vertices and {} edges'.format(
        tg.graph.num_vertices(), tg.graph.num_edges()))

    # Running the modified Normal Vector Voting algorithm:
    temp_normals_graph_file = '{}.VV_rh{}_normals.gt'.format(
        base_filename, radius_hit)
    method_tg_surf_dict = normals_directions_and_curvature_estimation(
        tg, radius_hit, methods=methods,
        page_curvature_formula=page_curvature_formula, area2=area2,
        poly_surf=surf, cores=cores, graph_file=temp_normals_graph_file)

    for method in list(method_tg_surf_dict.keys()):
        # Saving the output (TriangleGraph object) for later inspection in
        # ParaView:
        (tg, surf) = method_tg_surf_dict[method]
        if method == 'VV':
            if page_curvature_formula:
                method = 'NVV'
            elif area2:
                method = 'AVV'
            else:
                method = 'RVV'
        surf_file = '{}.{}_rh{}.vtp'.format(base_filename, method, radius_hit)
        io.save_vtp(surf, surf_file)
        gt_file = '{}.{}_rh{}.gt'.format(base_filename, method, radius_hit)
        tg.graph.save(gt_file)
        csv_file = '{}.{}_rh{}.csv'.format(base_filename, method, radius_hit)
        _extract_curvatures_from_graph(tg, csv_file)


def from_vtk_workflow(
        vtk_file, radius_hit, vertex_based, epsilon, eta, scale=(1, 1, 1),
        page_curvature_formula=False, methods=["VV"], area2=True, cores=6,
        reverse_normals=False):
    """
    Estimates curvature for each triangle in a triangle mesh in VTK format.

    Args:
        vtk_file (str): path to the VTK file with the surface
        radius_hit (float): radius in length unit of the graph, e.g. nanometers;
            it should be chosen to correspond to radius of smallest features of
            interest on the surface
        vertex_based (boolean): if True, curvature is calculated per triangle
            vertex instead of triangle center
        epsilon (float): parameter of Normal Vector Voting algorithm influencing
            the number of triangles classified as "crease junction" (class 2)
        eta (float): parameter of Normal Vector Voting algorithm influencing the
            number of triangles classified as "crease junction" (class 2) and
            "no preferred orientation" (class 3)
        scale (tuple, optional): pixel size (X, Y, Z) in given units for
            scaling the surface if it is not scaled (default (1, 1, 1))
        page_curvature_formula (boolean, optional): if True (default False),
            normal curvature formula from Page et al. is used in VV (see
            collect_curvature_votes)
        methods (list, optional): all methods to run in the second pass ('VV'
            and 'SSVV' are possible, default is 'VV')
        area2 (boolean, optional): if True (default), votes are weighted by
            triangle area also in the second step (principle directions and
            curvatures estimation; not possible if vertex_based is True)
        cores (int, optional): number of cores to run VV in parallel (default 6)
        reverse_normals (boolean, optional): if True (default False), original
            surface normals will be reversed

        Returns:
            None
        """
    if reverse_normals:
        reverse_normals_str = "_reversed_normals"
    else:
        reverse_normals_str = ""
    vtk_filename = os.path.basename(vtk_file)
    base_filename = os.path.splitext(vtk_filename)[0] + reverse_normals_str
    log_file = '{}.{}_rh{}_epsilon{}_eta{}.log'.format(
        base_filename, methods[0], radius_hit, epsilon, eta)
    sys.stdout = open(log_file, 'a')

    print('\nReading in the surface file to get a vtkPolyData surface...')
    surf = io.load_poly_from_vtk(vtk_file)

    # Running the modified Normal Vector Voting algorithm:
    normals_graph_file = '{}.VV_rh{}_epsilon{}_eta{}_normals.gt'.format(
        base_filename, radius_hit, epsilon, eta)
    method_tg_surf_dict = {}
    if not isfile(normals_graph_file):
        # Make or read in the graph first:
        if not vertex_based:
            triangle_graph_file = base_filename + ".gt"
            if not isfile(triangle_graph_file):
                # uses TriangleGraph's point_in_cells and triangle_cell_ids
                print('\nBuilding a triangle graph from the surface...')
                tg = TriangleGraph()
                tg.build_graph_from_vtk_surface(
                    surf, scale, reverse_normals=reverse_normals)
                if tg.graph.num_vertices() == 0:
                    raise pexceptions.PySegInputError(
                        expr="new_workflow", msg="The graph is empty!")
                print('The graph has {} vertices and {} edges'.format(
                    tg.graph.num_vertices(), tg.graph.num_edges()))
                tg.graph.save(triangle_graph_file)
            else:
                print('\nReading in the triangle graph from file...')
                tg = TriangleGraph()
                tg.graph = load_graph(triangle_graph_file)
            sg = tg
        else:  # vertex_based
            area2 = False
            point_graph_file = base_filename + "_point.gt"
            if not isfile(point_graph_file):
                print('\nBuilding a point graph from the surface...')
                pg = PointGraph()
                pg.build_graph_from_vtk_surface(
                    surf, scale, reverse_normals=reverse_normals)
                if pg.graph.num_vertices() == 0:
                    raise pexceptions.PySegInputError(
                        expr="new_workflow", msg="The graph is empty!")
                print('The graph has {} vertices and {} edges'.format(
                    pg.graph.num_vertices(), pg.graph.num_edges()))
                pg.graph.save(point_graph_file)
            else:
                print('\nReading in the point graph from file...')
                pg = PointGraph()
                pg.graph = load_graph(point_graph_file)
            sg = pg
        # Estimate normals, directions and curvatures:
        method_tg_surf_dict = normals_directions_and_curvature_estimation(
            sg, radius_hit, epsilon, eta, methods=methods,
            page_curvature_formula=page_curvature_formula,
            area2=area2, poly_surf=surf, cores=cores,
            graph_file=normals_graph_file)
    else:
        # Estimate directions and curvatures using the graph file with normals:
        for method in methods:
            sg_curv, surface_curv = curvature_estimation(
                radius_hit, graph_file=normals_graph_file, method=method,
                page_curvature_formula=page_curvature_formula, area2=area2,
                poly_surf=surf, cores=cores, vertex_based=vertex_based)
            method_tg_surf_dict[method] = (sg_curv, surface_curv)

    for method in list(method_tg_surf_dict.keys()):
        # Saving the output (TriangleGraph object) for later inspection in
        # ParaView:
        (sg_curv, surface_curv) = method_tg_surf_dict[method]
        if method == 'VV':
            if page_curvature_formula:
                method = 'NVV'
            elif area2:
                method = 'AVV'
            else:
                method = 'RVV'
        surf_file = '{}.{}_rh{}_epsilon{}_eta{}.vtp'.format(
            base_filename, method, radius_hit, epsilon, eta)
        io.save_vtp(surface_curv, surf_file)
        gt_file = '{}.{}_rh{}_epsilon{}_eta{}.gt'.format(
            base_filename, method, radius_hit, epsilon, eta)
        sg_curv.graph.save(gt_file)
        csv_file = '{}.{}_rh{}_epsilon{}_eta{}.csv'.format(
            base_filename, method, radius_hit, epsilon, eta)
        _extract_curvatures_from_graph(sg_curv, csv_file)


def from_nii_workflow(
        nii_file, outfold, radius_hit, page_curvature_formula=False,
        methods=["VV"], area2=True, cores=6):
    """
    Extracts surface for every label > 0 in the segmentation in NII format,
    after applying a Gaussian filter with sigma of 1.
    For each surface, estimates curvature for each triangle in a triangle mesh.

    Args:
        nii_file (str): NII file with the segmentation
        outfold (str): output folder
        radius_hit (float): radius in length unit of the graph, e.g. nanometers;
            it should be chosen to correspond to radius of smallest features of
            interest on the surface
        page_curvature_formula (boolean, optional): if True (default False),
            normal curvature formula from Page et al. is used in VV (see
            collect_curvature_votes)
        methods (list, optional): all methods to run in the second pass ('VV'
            and 'SSVV' are possible, default is 'VV')
        area2 (boolean, optional): if True (default), votes are weighted by
            triangle area also in the second step (principle directions and
            curvatures estimation)
        cores (int, optional): number of cores to run VV in parallel (default 6)

        Returns:
            None
        """
    base_filename = os.path.splitext(
        os.path.splitext(os.path.basename(nii_file))[0]
    )[0]  # without the path and without ".nii.gz" extensions
    log_file = '{}.{}_rh{}.log'.format(
        base_filename, methods[0], radius_hit)
    sys.stdout = open(log_file, 'a')

    # Reading in the data and getting the data type and average scaling in mm:
    seg, _, header = io.load_nii(nii_file)
    assert (isinstance(seg, np.ndarray))
    data_type = seg.dtype
    scale = header.get_zooms()
    print("pixel size in mm (x, y, z) = {}".format(scale))

    # Save as MRC file:
    mrc_file = str(os.path.join(outfold, base_filename+".mrc"))
    if not isfile(mrc_file):
        io.save_numpy(seg, mrc_file)

    for label in range(1, np.max(seg)+1):
        print("Label {}".format(label))
        # output base file name with the path and with the label:
        outfile_base = str(os.path.join(outfold, base_filename+str(label)))

        # Surface generation around the filled segmentation using
        # vtkMarchingCubes
        surf_file = outfile_base + ".surface.vtp"
        if not isfile(surf_file):
            filled_binary_seg = (seg == label).astype(data_type)
            if not np.any(filled_binary_seg):
                raise pexceptions.PySegInputError(
                    expr="from_nii_workflow",
                    msg="Label not found in the segmentation!")
            print("\nGenerating a surface...")
            surf = run_gen_surface(
                filled_binary_seg, outfile_base, lbl=1,
                other_mask=None, isosurface=True, sg=1, thr=THRESH_SIGMA1)
        else:
            print('\nReading in the surface from file...')
            surf = io.load_poly(surf_file)

        # Transforming the surface into a triangle graph
        print('\nBuilding a triangle graph from the surface...')
        tg = TriangleGraph()
        tg.build_graph_from_vtk_surface(surf, scale)
        if tg.graph.num_vertices() == 0:
            raise pexceptions.PySegInputError(
                expr="new_workflow", msg="The graph is empty!")
        print('The graph has {} vertices and {} edges'.format(
            tg.graph.num_vertices(), tg.graph.num_edges()))

        # Running the modified Normal Vector Voting algorithm:
        temp_normals_graph_file = '{}.VV_rh{}_normals.gt'.format(
            outfile_base, radius_hit)
        method_tg_surf_dict = normals_directions_and_curvature_estimation(
            tg, radius_hit, methods=methods,
            page_curvature_formula=page_curvature_formula, area2=area2,
            poly_surf=surf, cores=cores, graph_file=temp_normals_graph_file)

        for method in list(method_tg_surf_dict.keys()):
            # Saving the output (TriangleGraph object) for later inspection in
            # ParaView:
            (tg, surf) = method_tg_surf_dict[method]
            if method == 'VV':
                if page_curvature_formula:
                    method = 'NVV'
                elif area2:
                    method = 'AVV'
                else:
                    method = 'RVV'
            surf_file = '{}.{}_rh{}.vtp'.format(
                outfile_base, method, radius_hit)
            io.save_vtp(surf, surf_file)
            gt_file = '{}.{}_rh{}.gt'.format(outfile_base, method, radius_hit)
            tg.graph.save(gt_file)
            csv_file = '{}.{}_rh{}.csv'.format(outfile_base, method, radius_hit)
            _extract_curvatures_from_graph(tg, csv_file)


def main_javier(membrane, radius_hit):
    """
    Main function for running the new_workflow function for Javier's ER or PM.

    Args:
        membrane (string): what membrane segmentation to use 'ER' or 'PM'
        radius_hit (int): neighborhood parameter (in nm)

    Returns:
        None
    """
    t_begin = time.time()

    fold = "../experimental_data_sets/ER/AVV/"
    seg_file = "t2_ny01_lbl.labels_FILLED_half.mrc"
    base_filename = "TCB_180830_l2_t2half.{}".format(membrane)
    pixel_size = 1.368
    holes = 3
    min_component = 100
    runtimes_file = "{}{}_runtimes.csv".format(fold, base_filename)

    if membrane == "PM":
        lbl = 1
        print("\nEstimating normals for {}".format(base_filename))
        new_workflow(
            base_filename, seg_file, fold, pixel_size, radius_hit,
            methods=['VV'], label=lbl, holes=holes,
            min_component=min_component, only_normals=True,
            runtimes=runtimes_file)
    elif membrane == "ER":
        lbl = 2
        filled_lbl = 3  # ER lumen
        print("\nCalculating curvatures for {}".format(base_filename))
        new_workflow(
            base_filename, seg_file, fold, pixel_size, radius_hit,
            methods=['VV'], label=lbl, filled_label=filled_lbl,
            min_component=min_component, runtimes=runtimes_file)

        print("\nExtracting curvatures for {}".format(membrane))
        extract_curvatures_after_new_workflow(
            fold, base_filename, radius_hit, methods=['VV'],
            exclude_borders=1, categorize_shape_index=True)

        # surf_vtp_file = '{}{}.{}_rh{}.vtp'.format(
        #     fold, base_filename, 'AVV', radius_hit)
        # outfile_base = '{}{}.{}_rh{}'.format(
        #     fold, base_filename, 'AVV', radius_hit)
        # scale = (pixel_size, pixel_size, pixel_size)
        # seg = io.load_tomo(fold + seg_file)
        # size = seg.shape
        # convert_vtp_to_stl_surface_and_mrc_curvatures(
        #     surf_vtp_file, outfile_base, scale, size)
    else:
        print("Membrane not known.")
        exit(0)

    t_end = time.time()
    duration = t_end - t_begin
    minutes, seconds = divmod(duration, 60)
    print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))


def main_felix():
    """
    Main function for running the new_workflow function for Felix' data.

    Returns:
        None
    """
    t_begin = time.time()

    # Felix's vesicle:
    base_filename = "t74_vesicle3"
    pixel_size = 2.526
    radius_hit = 10  # nm
    fold = '../experimental_data_sets/vesicle/'
    seg_file = "t74_vesicle3_bin6.Labels.mrc"
    lbl = 1
    min_component = 100
    runtimes_file = "{}{}_runtimes.csv".format(fold, base_filename)

    print("\nCalculating curvatures for {}".format(base_filename))
    new_workflow(
        base_filename, seg_file, fold, pixel_size, radius_hit, methods=['VV'],
        label=lbl, holes=0, min_component=min_component, runtimes=runtimes_file)

    print("\nExtracting curvatures for vesicle")
    extract_curvatures_after_new_workflow(
        fold, base_filename, radius_hit, methods=['VV'],
        exclude_borders=1, categorize_shape_index=True)

    t_end = time.time()
    duration = t_end - t_begin
    minutes, seconds = divmod(duration, 60)
    print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))


def main_till(organelle, regions=True):
    """
    Main function for running the new_workflow function for Till' data.

    Args:
        organelle (string): run this organelle in so named subfolder.
        regions (boolean, optional): if True, split the segmentation into
        commented regions (minimal size of 1000 voxels) and run the workflow per
        region, each on 1 core. If False (default), use a ready surface.

    Returns:
        None
    """
    t_begin = time.time()

    # Input parameters:
    fold = '../experimental_data_sets/Golgi_and_vesicles/{}/'.format(organelle)
    base_filename = "l2_t6_{}".format(organelle)
    lbl = 1
    filled_lbl = 2
    min_comp = 1000  # to remove possible segmentation "noise"
    pixel_size = 1.684  # nm
    radius_hit = 10  # nm

    if regions:
        # Load the segmentation numpy array from the file and join both labels
        # as 1:
        seg_file = '20170217_FIB112_G1_l2_t6_dimifilt_lbl.{}Filled.mrc'.format(
            organelle)
        seg = io.load_tomo(fold + seg_file)
        assert (isinstance(seg, np.ndarray))
        data_type = seg.dtype
        binary_seg = ((seg == lbl) | (seg == filled_lbl)).astype(data_type)
        binary_seg_path = os.path.splitext(fold + seg_file)[0] + 'Binary.mrc'
        io.save_numpy(binary_seg, binary_seg_path)

        # Split the segmentation into binary regions, with 1 instead of both
        # labels:
        binary_seg_regions, _ = split_segmentation(
            infile=binary_seg_path, lbl=1, close=False,
            min_region_size=min_comp)

        cores = len(binary_seg_regions)
        seg_region_files = []
        base_region_files = []
        region_surf_files = []
        for i, binary_seg_region in enumerate(binary_seg_regions):
            # Restore the original labels at place of each binary region:
            seg_region = seg * binary_seg_region

            # Prepare files and filenames for the runs i/o:
            seg_region_file = os.path.splitext(seg_file)[0] + str(i+1) + '.mrc'
            io.save_numpy(seg_region, fold + seg_region_file)
            seg_region_files.append(seg_region_file)

            base_region_file = "{}{}".format(base_filename, str(i + 1))
            base_region_files. append(base_region_file)

            region_surf_file = '{}{}.AVV_rh{}.vtp'.format(
                fold, base_region_file, radius_hit)
            region_surf_files.append(region_surf_file)

        # Run each region on a different core:
        print("\nCalculating curvatures for {} regions".format(cores))
        p = pp.ProcessPool(cores)
        print('Opened a pool with {} processes'.format(cores))
        p.map(partial(new_workflow,
                      fold=fold, pixel_size=pixel_size, radius_hit=radius_hit,
                      methods=['VV'], page_curvature_formula=False, area2=True,
                      label=lbl, filled_label=filled_lbl, unfilled_mask=None,
                      holes=0, min_component=100, remove_wrong_borders=True,
                      only_normals=False, cores=1, runtimes=''),
              base_region_files, seg_region_files)
        p.close()
        p.clear()

        # Extract curvatures and join region VTP files to one VTP file:
        extract_curvatures_after_new_workflow(
            fold, base_filename, radius_hit, methods=['VV'],
            exclude_borders=1, categorize_shape_index=True, regions=cores)
        surf_file = '{}{}.AVV_rh{}.vtp'.format(fold, base_filename, radius_hit)
        io.merge_vtp_files(region_surf_files, surf_file)

    else:  # No regions, use the ready surface
        new_workflow(
            base_filename, "", fold=fold,
            pixel_size=pixel_size, radius_hit=radius_hit, methods=['VV'],
            page_curvature_formula=False, area2=True, label=lbl,
            filled_label=filled_lbl, unfilled_mask=None, holes=0,
            min_component=min_comp, remove_wrong_borders=True,
            only_normals=False, cores=6, runtimes='')

        # Extract curvatures:
        extract_curvatures_after_new_workflow(
            fold, base_filename, radius_hit, methods=['VV'],
            exclude_borders=2, categorize_shape_index=True)

    t_end = time.time()
    duration = t_end - t_begin
    minutes, seconds = divmod(duration, 60)
    print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))


def main_light_microscopy_cells(timepoint=22):
    """
    Main function for running the ply workflow for the whole C.Elegans embryo
    segmentation (by LimeSeg, of data coming from light microscopy) at a given
    time point. The input data structure of PLY files was downloaded from here:
    https://raw.githubusercontent.com/NicoKiaru/TestImages/master/LimeSegOutput/
    DubSeg.zip

    Args:
        timepoint (int): time point, 1-22 are possible for these data set
            (default 22)

    Returns:
        None
    """
    # Input / output settings:
    fold = "/fs/pool/pool-ruben/Maria/curvature/TestImages-LimeSeg/" \
           "LimeSegOutput/DubSeg/"
    cell_seg_file = "T_{}.ply".format(timepoint)
    scale_microns = 0.133
    radius_hit_microns = 3
    cell_vtp_file = "T_{}.AVV_rh{}.vtp".format(timepoint, radius_hit_microns)
    cell_vtp_file_paths_list = []
    embryo_vtp = os.path.join(fold, "embryo_T_{}.AVV_rh{}.vtp".format(
        timepoint, radius_hit_microns))

    # Run the curvature workflow:
    for subfold_p in [x for x in Path(fold).iterdir() if x.is_dir()]:
        cell_seg_file_path = os.path.join(fold, subfold_p.name, cell_seg_file)
        # As not all cells include every time point, check if the file exists
        if os.path.isfile(cell_seg_file_path):
            from_ply_workflow(
                ply_file=cell_seg_file_path, radius_hit=radius_hit_microns,
                scale=(scale_microns, scale_microns, scale_microns))
            # Add the output file path to the list, if it was created
            # (some surfaces were empty and skipped by the workflow):
            out_file_path = os.path.join(fold, subfold_p.name, cell_vtp_file)
            if os.path.isfile(out_file_path):
                cell_vtp_file_paths_list.append(out_file_path)

    # Merge the output VTP files to one:
    merge_vtp_files(cell_vtp_file_paths_list, embryo_vtp)


def main_heart():
    """
    Runs "from_nii_workflow" on the MRI heart data from HVSMR2016 challenge
    (http://segchd.csail.mit.edu/data.html).

    Returns:
     None
    """

    from_nii_workflow(
        nii_file="/fs/pool/pool-ruben/Maria/curvature/HVSMR2016_training_data/"
                 "GroundTruth/training_axial_full_pat0-label.nii.gz",
        outfold="/fs/pool/pool-ruben/Maria/curvature/HVSMR2016_training_data/"
                "GroundTruthOutput",
        radius_hit=5)


def main_brain():
    """
    Runs "from_vtk_workflow" with user-input parameters from a VTK surface file,
    like MRI brain provided by Mindboggle (https://osf.io/36gdy/).
    Comment it in below in the main and run from terminal like this:
    python curvature_calculation.py vtk_file radius_hit [epsilon eta]
    Either of the latter two parameters is optional and will be set to 0 if not
    given.

    Returns:
        None

    Note:
        The surface normals will be reversed to point inside the brain surface,
        according to our convention.
    """
    vtk_file = sys.argv[1]
    radius_hit = float(sys.argv[2])  # Mindboggle's default "radius disk"=2 mm
    if len(sys.argv) > 3:
        epsilon = float(sys.argv[3])
    else:
        epsilon = 0
    if len(sys.argv) > 4:
        eta = float(sys.argv[4])
    else:
        eta = 0
    from_vtk_workflow(
        vtk_file, radius_hit, vertex_based=False, epsilon=epsilon, eta=eta,
        scale=(1, 1, 1), methods=["VV"], cores=6, reverse_normals=True)


def main_pore(isosurface=False, radius_hit=2):
    """
    Main function for running the new_workflow function for a membrane pore.

    Args:
        isosurface (boolean): whether to generate isosurface or single-layered
            surface (default)
        radius_hit (int): neighborhood parameter (in nm, default 2)

    Returns:
        None
    """
    t_begin = time.time()

    fold = "../../../../curvature/4Antonio/"
    seg_file = "binarymembraneAmira.mrc"
    pixel_size = 0.26  # nm
    holes = 0
    min_component = 0
    lbl = 1

    if isosurface:
        base_filename = "binarymembrane_isosurface"
        runtimes_file = "{}{}_runtimes.csv".format(fold, base_filename)
        filled_lbl = 1
        print("\nCalculating curvatures for an isosurface")
        new_workflow(
            base_filename, seg_file, fold, pixel_size, radius_hit,
            methods=['VV'], label=lbl, filled_label=filled_lbl,
            min_component=min_component, runtimes=runtimes_file)

    else:
        base_filename = "binarymembrane_single_layer_surface"
        runtimes_file = "{}{}_runtimes.csv".format(fold, base_filename)
        print("\nCalculating curvatures for a single-layer surface")
        new_workflow(
            base_filename, seg_file, fold, pixel_size, radius_hit,
            methods=['VV'], label=lbl, holes=holes,
            min_component=min_component, runtimes=runtimes_file)

    extract_curvatures_after_new_workflow(
        fold, base_filename, radius_hit, methods=['VV'],
        exclude_borders=1, categorize_shape_index=True)

    t_end = time.time()
    duration = t_end - t_begin
    minutes, seconds = divmod(duration, 60)
    print('\nTotal elapsed time: {} min {} s'.format(minutes, seconds))


if __name__ == "__main__":

    main_javier('ER', 10)

    main_felix()

    # main_till('Golgi')
    # main_till('Vesicles')

    # main_light_microscopy_cells()

    # main_heart()

    # main_brain()

    # main_pore(isosurface=True, radius_hit=4)
