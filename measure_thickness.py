"""
Membrane thickness measurement for surface morphometrics pipeline.

Analyze density sampling CSV files to measure membrane thickness by fitting
dual gaussians to density profiles. Generates per-triangle thickness measurements,
refined surface files, and summary statistics.

Pipeline Citation: Barad BA*, Medina M*, Fuentes D, Wiseman RL, Grotjahn DA.
Quantifying organellar ultrastructure in cryo-electron tomography using a surface morphometrics pipeline.
J Cell Biol 2023.

Thickness Citation: Medina M, Chang Y-T, Rahmani H, Frank M, Khan Z, Fuentes D, Heberle FA, Waxham MN,
Barad BA, Grotjahn DA. Surface Morphometrics reveals local membrane thickness variation in organellar
subcompartments. J Cell Biol 2025.
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as opt
import scipy.signal as signal
import scipy.stats as stats
from scipy import spatial
from glob import glob
from pathlib import Path
import os
from sys import argv
import yaml
import click

from pycurv import TriangleGraph, io
from graph_tool import load_graph
from tqdm import tqdm
from morphometrics_stats import histogram
from intradistance_verticality import export_csv

def find_mins(y):
    """Find the indices of minimum values on left and right halves of the curve."""
    mid = int(np.round(len(y)/2))
    left_side = np.argmin(y[:mid])
    right_side = np.argmin(y[mid:]) + mid
    return left_side, right_side, y

def sinc(x, A, mu, sigma):
    return A * (np.sin(np.pi*(x-mu)*sigma) / (np.pi*(x-mu)*sigma))**2

def dual_sinc(x,p): # p = a1, mu1, sigma1, a2, mu2, sigma2, offset
    return sinc(x,*p[0:3])+sinc(x,*p[3:6])+p[6]

def gauss(x, p): # p[0]==mean, p[1]==stdev
    return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))

def monogaussian(x, h, c, w):
    return h*np.exp(-(x-c)**2/(2*w**2))
def dual_gaussian(x, h1, c1, w1, h2, c2, w2, o): # p = h1, c1, w1, h2, c2, w2, offset
    return monogaussian(x,h1,c1,w1)+monogaussian(x,h2,c2,w2)+o

# Fit a gaussian to a series of 21 points and return the thickness
def fit_gaussian(x, thickness_set, skipedge=3):
    x = x[skipedge:-1*skipedge]
    thickness_set = thickness_set[skipedge:-1*skipedge]
    p0 = [0, 4]
    errfunc = lambda p,a,b: gauss(a,p)-b
    p1, success = opt.leastsq(errfunc, p0[:], args=(x,thickness_set))
    fwhm = 2*np.sqrt(2*np.log(2))*p1[1]
    return p1, fwhm


def func(x, *args):
    x = x.reshape(-1, 1)
    a = np.array(args[0::2]).reshape(1, -1)
    b = np.array(args[1::2]).reshape(1, -1)
    return np.sum(a * np.exp(-b * x), axis=1)



def peak_fit(x, y):
    """Use scipy peak width to return a FWHM"""
    # peak = np.argmax(y) 
    # peaks = [peak]
    # peaks, _ = signal.find_peaks(y, 0.13)
    peaks = [np.argmax(y)]
    results_half = signal.peak_widths(y, peaks, rel_height=0.5)
    width = results_half[0][0]
    if width == 0:
        return -1,0,0,0
    height = results_half[1][0]
    h0 = results_half[2][0]+x[0]
    h1=  results_half[3][0]+x[0]
    return width, height, h0, h1
    

def find_two_peaks(x,y):
    peaks, _ = signal.find_peaks(y)
    if len(peaks)<2:
        return 0, 0, 0
    peaks = peaks[::-1]
    pos = np.argsort(y[peaks])
    peaks = np.take_along_axis(peaks, pos, axis=0)
    peak1 = x[peaks[-1]]
    peak2 = x[peaks[-2]]
    width = np.abs(peak1-peak2)

    return width, peak1, peak2


def process_single_surface(filename, average_radius, output_dir):
    """
    Process a single thickness sampling file and generate plots and refined surfaces.

    Parameters
    ----------
    filename : str
        Path to the thickness sampling CSV file
    average_radius : float
        Radius for local averaging in thickness calculations
    output_dir : str
        Directory for output files

    Returns
    -------
    tuple
        (component_info_dict, per_surface_thickness, areas, width, x, fig2, ax2)
    """
    tsname = Path(filename).stem.split(".")[0]
    comp_num = Path(filename).stem.split("_")[-5].split(".")[0]
    print(f"Processing {tsname}")

    # Load thickness csv file (header contains relative positions)
    thickness_set = pd.read_csv(filename, header=0)

    # Extract x positions from column headers
    x = np.array([float(col) for col in thickness_set.columns])
    graph_file = filename[:-13] + ".gt"
    csv_outfile = filename[:-13] + ".csv"
    graph_file_final = filename[:-13] + "_refined.gt"

    tg = TriangleGraph()
    tg.graph = load_graph(graph_file)
    print(tg.graph.vp.points[0], tg.graph.vp.xyz[0])

    areas = tg.graph.vp.area.get_array()
    xx, yy, zz = tg.graph.vp.xyz.get_2d_array([0, 1, 2])
    nvx, nvy, nvz = tg.graph.vp.n_v.get_2d_array([0, 1, 2])
    xyz = tg.graph.vp.xyz.get_2d_array([0, 1, 2]).transpose()
    xyztree = spatial.cKDTree(xyz)

    avg_x = np.average(xx, weights=areas)
    avg_y = np.average(yy, weights=areas)
    avg_z = np.average(zz, weights=areas)
    total_area = np.sum(areas)

    curvedness = tg.graph.vp.curvedness_VV.get_array()
    rad_curv = [1/i for i in curvedness]
    rad_avg = np.average(rad_curv, weights=areas)
    rad_std = np.sqrt(np.cov(rad_curv, aweights=areas))

    surface_file = filename[:-13] + "_refined.vtp"
    per_surface_thickness = []
    per_triangle_offset = []

    fig2, ax2 = plt.subplots()

    # Per-triangle thickness calculation
    for i in tqdm(range(len(rad_curv))):
        l, neighbors = xyztree.query(xyz[i], k=500, distance_upper_bound=average_radius, workers=-1)
        neighbors = neighbors[np.where(l != np.inf)]
        l = l[np.where(l != np.inf)]
        weights = [1/(1+j) for j in l]
        dat = np.asarray(np.average(thickness_set.iloc[neighbors], weights=weights, axis=0)) * -1
        dat = dat - min(dat)
        dat = dat / (80/81 * sum(dat))

        ipk = x[np.argmax(dat)]
        p0 = [0.02, ipk-0.5, 1.5, 0.02, ipk+0.5, 1.5, 0]
        bounds = ([0.005, ipk-6, 0.8, 0.005, ipk-1.5, 0.8, -1],
                  [0.04, ipk+1.5, 2.2, 0.04, ipk+6, 2.2, 1])
        mins = find_mins(dat)

        a = x[mins[0]+2:mins[1]-2]
        b = dat[mins[0]+2:mins[1]-2]

        try:
            p3, _ = opt.curve_fit(dual_gaussian, a, b, p0, bounds=bounds)
            if i % 5000 == 0:
                ax2.plot(x, dat)
            per_surface_thickness.append(np.abs(p3[4]-p3[1]))
            per_triangle_offset.append((p3[4]+p3[1])/2)
        except Exception:
            per_surface_thickness.append(np.nan)
            per_triangle_offset.append(0)

    # Average thickness calculation
    avg = thickness_set.mean(axis=0) * -1
    avg = avg - min(avg)
    avg = avg / (80/81 * sum(avg))
    mins = find_mins(avg)

    ipk = x[np.argmax(avg)]
    p0 = [0.02, ipk-0.5, 1.5, 0.02, ipk+0.5, 1.5, 0]
    bounds = ([0.005, ipk-6, 0.8, 0.005, ipk-2, 0.8, -1],
              [0.04, ipk+2, 2.2, 0.04, ipk+6, 2.2, 1])

    a = x[mins[0]+2:mins[1]-2]
    b = avg[mins[0]+2:mins[1]-2]

    p3, _ = opt.curve_fit(dual_gaussian, a, b, p0, bounds=bounds)
    width = np.abs(p3[1]-p3[4])

    # Update graph with thickness properties
    average_width_prop = tg.graph.new_vertex_property("float")
    average_width_prop.a = [width] * len(thickness_set)
    thick = tg.graph.new_vertex_property("float")
    thick.a = per_surface_thickness
    offset = tg.graph.new_vertex_property("float")
    offset.a = per_triangle_offset

    tg.graph.vp.average_width = average_width_prop
    tg.graph.vp.thickness = thick
    tg.graph.vp.offset = offset
    tg.graph.save(graph_file_final)

    surf = tg.graph_to_triangle_poly()
    io.save_vtp(surf, surface_file)
    export_csv(tg, csv_outfile)

    del tg

    # Component info for CSV output
    component_info = {
        'tsname': tsname,
        'comp_num': comp_num,
        'avg_x': avg_x,
        'avg_y': avg_y,
        'avg_z': avg_z,
        'total_area': total_area,
        'rad_avg': rad_avg,
        'rad_std': rad_std,
        'width': width,
        'p3': p3,
        'avg': avg
    }

    return component_info, per_surface_thickness, areas/np.sum(areas), width, x, fig2, ax2


def run_measure_thickness(config, output_dir=None):
    """
    Main function to run thickness plotting and analysis.

    Parameters
    ----------
    config : dict
        Configuration dictionary loaded from config.yml
    output_dir : str or None
        Output directory for plots (defaults to work_dir)
    """
    # Get settings from config
    work_dir = config.get("work_dir", config.get("data_dir", "./"))
    if not work_dir.endswith("/"):
        work_dir += "/"

    if output_dir is None:
        output_dir = work_dir
    elif not output_dir.endswith("/"):
        output_dir += "/"

    # Get thickness settings
    thickness_config = config.get("thickness_measurements", {})
    components = thickness_config.get("components", [])
    average_radius = thickness_config.get("average_radius", 12)
    radius_hit = config.get("curvature_measurements", {}).get("radius_hit", 9)
    # Note: nsamples and scan_range are now read from CSV headers (set by sample_density.py)

    if not components:
        print("No components specified in config.yml thickness_measurements section")
        return

    print(f"Thickness plots settings:")
    print(f"  Work directory: {work_dir}")
    print(f"  Output directory: {output_dir}")
    print(f"  Components: {components}")
    print(f"  Average radius: {average_radius}")
    print(f"  Radius hit: {radius_hit}")

    # Find files for each component
    filenames = {component: [] for component in components}
    for component in components:
        basename = work_dir + f"*{component}.AVV_rh{radius_hit}_sampling.csv"
        fileset = glob(basename)
        filenames[component].extend(fileset)
        print(f"  Found {len(fileset)} files for {component}")

    # x positions are now read from CSV headers in process_single_surface

    thickness_measurements = []
    area_measurements = []
    widths = {component: [] for component in components}

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    component_list_file = os.path.join(output_dir, "component_list.csv")
    with open(component_list_file, "w") as compfile:
        compfile.write("TS,Component Type,Component Number,Centroid X,Centroid Y,Centroid Z,"
                       "Total Area,Radius of Curvature,Rad_Curv STD,Thickness,"
                       "Peak1 Position,Peak1 Sigma,Peak2 Position,Peak2 Sigma\n")

        for component in components:
            fig, ax = plt.subplots()
            print(f"\nProcessing component: {component}")

            for index, filename in enumerate(filenames[component]):
                info, per_surface_thickness, norm_areas, width, x, fig2, ax2 = \
                    process_single_surface(filename, average_radius, output_dir)

                thickness_measurements.append(per_surface_thickness)
                area_measurements.append(norm_areas)
                widths[component].append(width)

                # Write component info
                p3 = info['p3']
                compfile.write(f"{info['tsname']},{component},{info['comp_num']},"
                               f"{info['avg_x']:.1f},{info['avg_y']:.1f},{info['avg_z']:.1f},"
                               f"{info['total_area']:.1f},{info['rad_avg']:.2f},{info['rad_std']:.2f},"
                               f"{width:.2f},{p3[1]:.2f},{p3[2]:.2f},{p3[4]:.2f},{p3[5]:.2f}\n")

                # Generate individual surface plot
                fig3, ax3 = plt.subplots()
                scan_range = (x[-1] - x[0]) / 2  # Derive scan_range from x positions
                ax.plot(x, info['avg'], label=index)
                ax3.plot(x, info['avg'], label="data")
                ax3.plot(x, dual_gaussian(x, *p3), "-.", label="Dual Gaussian Fit")
                ax3.plot(x, monogaussian(x, *p3[0:3]) + p3[6], "--", label="Gaussian 1")
                ax3.plot(x, monogaussian(x, *p3[3:6]) + p3[6], "--", label="Gaussian 2")
                ax3.axvspan(p3[1], p3[4], facecolor='g', alpha=0.1, label="Dual Gauss Span")
                ax3.legend()
                ax3.set_xlabel("Distance (nm)")
                ax3.set_xlim(x[0], x[-1])
                ax3.set_ylabel("Density")
                ax3.set_title(f"{info['tsname']} - {component} {info['comp_num']}")

                fig2.savefig(os.path.join(output_dir, f"thickness_{component}.png"))
                fig3.savefig(os.path.join(output_dir,
                             f"thickness_average_{info['tsname']}_{component}_{info['comp_num']}_fit.svg"))
                plt.close(fig=fig2)
                plt.close(fig=fig3)

            # Save component summary plot
            ax.set_xlabel("Distance (nm)")
            if 'x' in dir():  # Use x from last processed file
                ax.set_xlim(x[0], x[-1])
            ax.set_ylabel("Density")
            ax.legend()
            ax.set_title(f"{component} - All Curves")
            fig.savefig(os.path.join(output_dir, f"{component}_Averages.png"))
            fig.clear()
            plt.close()

    # Print summary statistics
    thicknesses = []
    for component in components:
        if widths[component]:
            thicknesses.append(widths[component])
            print(f"{component} - mean: {np.mean(widths[component]):.3f}, "
                  f"stdev: {np.std(widths[component]):.3f}")

    # Statistical tests (if we have at least 2 components with data)
    if len(thicknesses) >= 2:
        res = stats.ttest_ind(thicknesses[0], thicknesses[1])
        print(f"Student's T Test pval: {res.pvalue}, df: {res.df}")

    # Generate violin plot
    if thicknesses:
        fig4, ax4 = plt.subplots()
        ax4.violinplot(thicknesses, showmedians=True)
        ax4.set_xticks(range(1, len(components)+1))
        ax4.set_xticklabels(components)
        fig4.savefig(os.path.join(output_dir, "violin.svg"))

        with open(os.path.join(output_dir, "violin.csv"), "w") as violin:
            for component in components:
                if widths[component]:
                    vstring = ",".join([str(i) for i in widths[component]])
                    violin.write(f"{component},{vstring}\n")

    # Generate histogram
    if thickness_measurements and area_measurements:
        histogram(data=thickness_measurements, areas=area_measurements,
                  labels=components, title="Thickness Comparison", xlabel="Thickness (nm)")

    print(f"\nOutput files written to {output_dir}")


@click.command()
@click.argument('configfile', type=click.Path(exists=True))
@click.option('--output', '-o', type=str, default=None,
              help='Output directory for plots (defaults to work_dir from config)')
@click.option('--average_radius', type=float, default=None,
              help='Radius for local averaging (overrides config)')
def measure_thickness_cli(configfile, output, average_radius):
    """
    Measure membrane thickness from density sampling data.

    Fits dual gaussians to density profiles to estimate membrane thickness
    for each triangle in the surface mesh. Generates refined surface files
    with thickness properties and summary statistics.

    CONFIGFILE: Path to the config.yml file
    """
    with open(configfile) as f:
        config = yaml.safe_load(f)

    # Override config settings if specified on command line
    if average_radius is not None:
        if "thickness_measurements" not in config:
            config["thickness_measurements"] = {}
        config["thickness_measurements"]["average_radius"] = average_radius

    run_measure_thickness(config, output_dir=output)


if __name__ == "__main__":
    measure_thickness_cli()
