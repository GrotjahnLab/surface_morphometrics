#! /usr/bin/env python
"""Run pycurv on a vtp surface file

Usage: curvarture.py <vtp_file>"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import time
import sys
import click

from curvature_calculation import new_workflow, extract_curvatures_after_new_workflow

@click.command()
@click.argument('filename')
@click.argument('folder')
@click.option('--scale', default=1, help='Pixel size scaling for the vtp file. Our surfaces are generally already scaled, so leave at 1. Default is 1.')
@click.option('--radius_hit', default=10, help='Radius for vector voting in the curvature calculation. Default is 10.')
@click.option('--min_component', default=30, help='Isolated surfaces smaller than this will get removed before calculation. Default is 30.')
@click.option('--exclude_borders', default=0, help='Distance in surface units (angstroms or nm) to exclude from the curvature calculation. Default is 0.')
@click.option('--cores', default=16, help='Number of cores to use for the calculation. Default is 16.')
@click.option('--remove_wrong_borders', default=False, help='Eat back the surface before calculations. If using screened poisson workflow leave this False.')
@click.option('-f', '--force', is_flag=True, default=False, help='Automatically pass all user input prompts.')
def run_pycurv_cli(filename, folder, scale, radius_hit, min_component, exclude_borders, cores, remove_wrong_borders, force):
    logical_cores = os.cpu_count()
    if cores > logical_cores:
        click.echo(f"WARNING: Requested cores ({cores}) exceeds the number of logical cores ({logical_cores}).")
        click.echo("This may cause performance degradation due to oversubscription.")
        if not force:
            if not click.confirm("Continue anyway?", default=False):
                raise SystemExit(1)
    run_pycurv(filename, folder, scale, radius_hit, min_component, exclude_borders, cores, remove_wrong_borders)

def run_pycurv(filename, folder, scale=1, radius_hit=10, min_component=30, exclude_borders=0, cores=16, remove_wrong_borders=False):
    """Run pycurv on a vtp surface file and extract curvatures
    
    filename (str): vtp surface file
    folder (str): folder containing the vtp file
    scale (float): pixel size scaling for the vtp file. Our surfaces are generally already scaled, so leave at 1. Default is 1.
    radius_hit (float): radius for vector voting in the curvature calculation. Default is 10.
    min_component (int): isolated surfaces smaller than this will get removed before calculation. Default is 30.
    exclude_borders (int): distance in surface units (angstroms or nm) to exclude from the curvature calculation. Default is 0.
    cores (int): number of cores to use for the calculation. Default is 16.
    remove_wrong_borders (bool): eat back the surface before calculations. If using screened poisson workflow leave this False.
    """
    assert filename.endswith(".surface.vtp"), "Surface must be a vtp file of a surface, ending with .surface.vtp"
    basename = filename[:-len(".surface.vtp")]
    runtimes_file = "{}{}_runtimes.csv".format(folder, basename)
    seg_file = ""
    t_begin = time.time()
    print("\nCalculating curvatures for {}".format(basename))
    new_workflow(
        basename, seg_file, folder, scale, radius_hit, methods=['VV'],  min_component=min_component, runtimes=runtimes_file, cores=cores, remove_wrong_borders=remove_wrong_borders)
    print("\nExtracting curvatures for all surfaces")
    extract_curvatures_after_new_workflow(
        folder, basename, radius_hit, methods=['VV'],
        exclude_borders=exclude_borders, categorize_shape_index=True)
    t_end = time.time()
    duration = t_end - t_begin
    minutes, seconds = divmod(duration, 60)
    if not folder.endswith("/"):
        folder += "/"
    output_vtp = folder+basename+f'.AVV_rh{radius_hit}.vtp'
    output_csv = folder+basename+f'.AVV_rh{radius_hit}.csv'
    output_gt = folder+basename+f'.AVV_rh{radius_hit}.gt'
    output_log = folder+basename+f'.VV_rh{radius_hit}.log'

    print('\nTotal pycurv time: {} min {} s'.format(minutes, seconds))
    sys.stdout = sys.__stdout__
    print("Final outputs written:")
    print("VTP file for paraview: "+output_vtp)
    print("CSV file for pandas based quantification: "+output_csv)
    print("GT file for further morphometrics quantification: "+output_gt)
    print("Log for troubleshooting: "+output_log)
    
    return

if __name__=="__main__":
    run_pycurv_cli()