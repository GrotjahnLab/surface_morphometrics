## Load a thickness csv file and make a bunch of plots of the samples, and fit a gaussian to estimate FWHM

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as opt
import scipy.signal as signal
import scipy.stats as stats
from scipy import spatial
from glob import glob
from pathlib import Path
from pycurv import  TriangleGraph, io
from graph_tool import load_graph
from tqdm import tqdm
from morphometrics_stats import histogram
from intradistance_verticality import export_csv

base_folder = "/Users/bbarad/Downloads/TE/morphometrics/"
components = ["OMM", "IMM", "ER"]
# filenames = [base_folder + f"tomo_1_{component}_sampling_isonet.csv" for component in components]
AVERAGE_RAD = 12
def find_mins(y):
    mid = int(np.round(len(y)/2))
    left_side = np.argmin(y[:mid])
    # for i in range(left_side):
    #     y[i] = y[left_side]
    right_side = np.argmin(y[mid:])+mid
    # for i in range(right_side, len(y)):    
    #     y[i] = y[right_side]
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


filenames = {component:[] for component in components}
outfiles = {component: [] for component in components}
for component in components:
    basename = base_folder + f"*{component}.AVV_rh8_sampling.csv"
    fileset = glob(basename)
    # fileset.sort(key=lambda x: int(Path(x).stem.split("_")[-4].split(".")[0]))
    filenames[component].extend(fileset)


# voxsize = 1
x = np.linspace(-10,10,81)
# x = x/voxsize
thickness_measurements = []
area_measurements = []
widths = {component: [] for component in components}
with open("component_list.csv", "w") as compfile:
    compfile.write("TS,Component Type,Component Number,Centroid X,Centroid Y,Centroid Z,Total Area,Radius of Curvature,Rad_Curv STD,Thickness,Peak1 Position,Peak1 Sigma,Peak2 Position,Peak2 Sigma\n")
    for index2, component in enumerate(components):
        fig, ax = plt.subplots()
        print(component)
        for index, filename in enumerate(filenames[component]):
            # print(index, filename)
            fig2,ax2 = plt.subplots()
            fig3,ax3 = plt.subplots()
            tsname = Path(filename).stem.split(".")[0]
            print(tsname)
            comp_num = Path(filename).stem.split("_")[-5].split(".")[0]
            print(comp_num)
            # Load thickness csv file (headerless)
            thickness_set = pd.read_csv(filename, header=None)
            graph_file = filename[:-13]+"_refined.gt"
            csv_outfile = filename[:-13]+".csv"
            graph_file_final = filename[:-13]+"_refined.gt"
            tg = TriangleGraph()
            tg.graph=load_graph(graph_file)
            print(tg.graph.vp.points[0], tg.graph.vp.xyz[0])
            areas = tg.graph.vp.area.get_array()
            area_measurements.append(areas/np.sum(areas))
            xx,yy,zz = tg.graph.vp.xyz.get_2d_array([0,1,2])
            nvx,nvy,nvz = tg.graph.vp.n_v.get_2d_array([0,1,2])
            xyz = tg.graph.vp.xyz.get_2d_array([0,1,2]).transpose()
            xyztree = spatial.cKDTree(xyz)
            avg_x = np.average(xx, weights=areas)
            avg_y = np.average(yy, weights=areas)
            avg_z = np.average(zz, weights=areas)
            total_area = np.sum(areas)
            curvedness = tg.graph.vp.curvedness_VV.get_array()
            rad_curv = [1/i for i in curvedness]
            rad_avg = np.average(rad_curv, weights=areas)
            rad_std = np.sqrt(np.cov(rad_curv, aweights=areas))
            # print(np.average(rad_curv, weights=areas), np.sqrt(np.cov(rad_curv, aweights=areas)))
            surface_file = filename[:-13]+"_refined.vtp"
            per_surface_thickness = []
            per_triangle_offset = []
            for i in tqdm(range(len(rad_curv))):
                # if i - AVERAGER/2 < 0:
                #     imin = 0
                #     imax = AVERAGER
                # elif i+AVERAGER/2 >= len(rad_curv):
                #     imin=-AVERAGER
                #     imax=-1
                # else:
                #     imin = int(i - AVERAGER/2)
                #     imax = int(i+AVERAGER/2)
                l,neighbors = xyztree.query(xyz[i],k=500, distance_upper_bound=AVERAGE_RAD, workers=-1)
                neighbors = neighbors[np.where(l != np.inf)]
                l = l[np.where(l != np.inf)]
                weights = [1/(1+j) for j in l]
                dat = np.asarray(np.average(thickness_set.iloc[neighbors], weights=weights, axis=0))*-1
                dat = dat-min(dat)
                dat = dat/(80/81*sum(dat))

                ipk = x[np.argmax(dat)]
                p0 = [0.02,ipk-0.5,1.5,0.02,ipk+0.5,1.5,0] # initial guess - positions 1 and 4 are the critical centers!
                bounds = ([0.005,ipk-6,0.8,0.005,ipk-2,0.8,-1],[0.04,ipk+2,2.2,0.04,ipk+6,2.2,1]) # range of 8 for each center (up to 12 apart total), 0.7 to 2.2 thickness 
                mins = find_mins(dat)

                errfunc = lambda p,a,b: (dual_gaussian(a,*p)-b)**2
                a = x[mins[0]+2:mins[1]-2]
                b = dat[mins[0]+2:mins[1]-2]
                
                # p1,success = opt.leastsq(errfunc, p0, args=(x,y2))
                # p2,success2 = opt.leastsq(errfunc, p0, args=(x,y2))
                try:
                    p3, _ = opt.curve_fit(dual_gaussian, a, b, p0, bounds=bounds)
                    if i % 5000 == 0:
                        # print(p3[4]-p3[1])
                        ax2.plot(x, dat)

                    per_surface_thickness.append(np.abs(p3[4]-p3[1]))
                    per_triangle_offset.append((p3[4]+p3[1])/2)
                    # per_surface_thickness.append(1)
                    # per_triangle_offset.append(-1)
                except:
                    per_surface_thickness.append(np.nan)
                    per_triangle_offset.append(0)
            # print(per_triangle_offset)
            thickness_measurements.append(per_surface_thickness)
            # temp_thicknesses = []
            # for _, row in thickness_set.iterrows():
            #     y = np.asarray(row)*-1
            #     # y = y-min(y)
            #     y = y/(80/81*sum(y))
            #     width, _,_,_ = peak_fit(x,y)
            #     # print(width)
            #     if width != -1:
            #         temp_thicknesses.append(width)
            # print(index, np.median(temp_thicknesses),np.std(temp_thicknesses))
            # thicknesses.append(temp_thicknesses)
            avg = thickness_set.mean(axis=0)*-1
            avg = avg-min(avg)
            avg = avg/(80/81*sum(avg))
            mins = find_mins(avg)
            y2 = mins[2]
            # string = ",".join([str(i) for i in avg])+"\n"
            # compcsv.write(string)
            # p1, fwhm = fit_gaussian(x, avg)
            # fit_width, height, h0,h1 = peak_fit(x, avg)
            # print(index, fwhm, fit_width)
            # width, peak1, peak2 = find_two_peaks(x, avg)
            # if width > 0:
            #     widths[component].append(width)
            # print(width, peak1, peak2)
            ipk = x[np.argmax(avg)]
            p0 = [0.02,ipk-0.5,1.5,0.02,ipk+0.5,1.5,0] # initial guess - positions 1 and 4 are the critical centers!
            bounds = ([0.005,ipk-6,0.8,0.005,ipk-2,0.8,-1],[0.04,ipk+2,2.2,0.04,ipk+6,2.2,1]) # range of 8 for each center (up to 12 apart total), 0.7 to 2.2 thickness 

            errfunc = lambda p,a,b: (dual_gaussian(a,*p)-b)**2
            a = x[mins[0]+2:mins[1]-2]
            b = avg[mins[0]+2:mins[1]-2]
            
            # p1,success = opt.leastsq(errfunc, p0, args=(x,y2))
            # p2,success2 = opt.leastsq(errfunc, p0, args=(x,y2))
            p3, _ = opt.curve_fit(dual_gaussian, a, b, p0, bounds=bounds)
            # print(p3)
            # print(np.abs(p1[1]-p1[4]), p1[1], p1[4])
            width = np.abs(p3[1]-p3[4])
            # TS,Component Type,Component Number,Centroid X,Centroid Y,Centroid Z,Total Area,Radius of Curvature,Rad_Curv STD,Thickness,Peak1 Position,Peak1 Sigma,Peak2 Position,Peak2 Sigma
            compfile.write(f"{tsname},{component},{comp_num},{avg_x:.1f},{avg_y:.1f},{avg_z:.1f},{total_area:.1f},{rad_avg:.2f},{rad_std:.2f},{width:.2f},{p3[1]:.2f},{p3[2]:.2f},{p3[4]:.2f},{p3[5]:.2f}\n")
            average_width_prop = tg.graph.new_vertex_property("float")
            average_width_prop.a = [width]*len(thickness_set)
            thick = tg.graph.new_vertex_property("float")
            thick.a = per_surface_thickness
            # Change xx, yy, and zz by the normal vector x,y,z times the offset
            offset = tg.graph.new_vertex_property("float")
            offset.a = per_triangle_offset
            newxyz = np.array([xx-(nvx * per_triangle_offset), yy-(nvy*per_triangle_offset), zz-(nvz * per_triangle_offset)]).transpose()
            # xyzdict = {xyz[i]:j for i,j in enumerate(newxyz)}
            # triangles = tg.graph.vp.vertices.get_2d_array()
            # print(triangles[0])
            # for i in range(len(triangles)):
            #     for j in range(3):
            #         triangles[i][j] = xyzdict[triangles[i][j]]
            # tg.graph.gp.triangle_points.set_2d_array(triangles)
            
            # print(newxyz.shape, tg.graph.vp.xyz.get_2d_array().shape)
            # tg.graph.vp.xyz.set_2d_array(newxyz.transpose())
            tg.graph.vp.average_width = average_width_prop
            tg.graph.vp.thickness = thick
            tg.graph.vp.offset = offset
            tg.graph.save(graph_file_final)
            surf = tg.graph_to_triangle_poly()
            io.save_vtp(surf, surface_file)
            export_csv(tg, csv_outfile)
            
            del tg

            widths[component].append(np.abs(p3[1]-p3[4]))

        

            ax.plot(x, avg, label=index)
            ax3.plot(x, avg, label="data")
            # fig3.savefig(f"thickness_average{component}{index}.png")
            # ax3.plot(x, gauss(x,p1), '-.', label="Gaussian fit")
            ax3.plot(x, dual_gaussian(x,*p3), "-.", label="Dual Gaussian Fit")
            ax3.plot(x, monogaussian(x, *p3[0:3])+p3[6], "--", label="Gaussian 1")
            ax3.plot(x, monogaussian(x, *p3[3:6])+p3[6], "--", label="Gaussian 2")

            ax3.axvspan(p3[1], p3[4], facecolor='g', alpha=0.1, label="Dual Gauss Span")

            ax3.legend()
            # ax2.set_xlabel("Distance (nm)")
            # ax2.set_xlim(-10,10)
            # ax2.set_ylabel("Density")
            ax3.set_xlabel("Distance (nm)")
            ax3.set_xlim(-10,10)
            ax3.set_ylabel("Density")
            # ax2.set_title(f"Component {component} {index}")

            ax3.set_title(f"{tsname} - {component} {comp_num}")
            # thickness = fit_gaussian(x, thickness_set)
            fig2.savefig(f"thickness_{component}.png")
            fig3.savefig(f"thickness_average_{tsname}_{component}_{comp_num}_fit.svg")
            ax3.legend()
            # ax3.hlines(y=np.max(avg),xmin=peak1,xmax=peak2, linewidth=2, linestyle='-', color="r", alpha=0.5, label="Peakfinder FWHM")
            ax3.legend()
            # fig3.savefig(f"thickness_average_{component}_{index}_fit_twice.png")
            # plt.close(fig=fig2)
            plt.close(fig=fig3)
        ax.set_title("All component averages")
        ax.set_xlabel("Distance (nm)")
        ax.set_xlim(-10,10)
        ax.set_ylabel("Density")
        ax.legend()
        ax.set_title(f"{component} - All Curves")
        fig.savefig(f"{component}_Averages.png")
        fig.clear()
        plt.close()
    
thicknesses = []
for component in components:
    thicknesses.append(widths[component])
    print(f"{component} - mean: {np.mean(widths[component])}, stdev: {np.std(widths[component])}")

res = stats.ttest_ind(thicknesses[0], thicknesses[1])
print(f"Student's T Test pval: {res.pvalue}, df: {res.df}")


fig4,ax4 = plt.subplots()
ax4.violinplot(thicknesses, showmedians=True)
ax4.set_xticks(range(1,len(components)+1))
ax4.set_xticklabels(components)
fig4.savefig("violin.svg")
with open("violin.csv", "w") as violin:
    for component in components:
        vstring = ",".join([str(i) for i in widths[component]])
        violin.write(f"{component},{vstring}\n")

histogram(data=thickness_measurements, areas=area_measurements, labels=components,title="Thickness Comparison", xlabel="Thickness (nm)")