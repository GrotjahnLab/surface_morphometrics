import copy
import glob
import os
import pickle

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from textwrap import wrap
SMALL_SIZE = 8
MEDIUM_SIZE = 9
BIGGER_SIZE = 10.5

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

import numpy as np
import pandas as pd
import scipy.stats as stats

purple = [.5, 0, .5, 0.9]
purple_light = [.5,0,.5, 0.3]
green = [0, 1, 0, 0.9]
green_light = [0,1,0,.3]
blue = [0,.3,1, .9]
blue_light = [0,.3,1,.3]
junction_minimum = 19 # nm 
junction_maximum = 40 # nm

my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap.colors[0])


colors = [purple,green, blue]
colors_light = [purple_light, green_light, blue_light]

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

class Experiment():
    """Experiments are containers for a series of tomograms and their pandas dataframes"""
    def __init__(self, name):
        self.name = name
        self.tomograms = {}
        self.tomogram_names = set()
    
    # def __init__(self, names):
    #     """Initialize with a list of tomograms by basename (aka TT6, not TT6_labels.mrc). Must be the basename used by pycurv."""
    #     self.tomograms = {}
    #     self.tomogram_names = set(names)

    #     for name in self.tomogram_names:
    #         self.tomograms[name] = Tomogram(name)

    def add_tomograms(self, names, surface_bases, folder, file_extension):
        """Add tomograms to the experiment. Names must be the basenames used by pycurv."""
        for name in names:
            # print(name)
            surfs = []
            fns = []
            for surface_base in surface_bases:
                filename = folder+name+"_"+surface_base+file_extension
                if os.path.isfile(filename):
                    surfs.append(surface_base)
                    fns.append(filename)
            # print(surface_bases, fns)
            self.tomograms[name] = Tomogram(name, surface_bases, fns)

    def __getitem__(self, name):
        return self.tomograms[name]
    
    def __setitem__(self, name, tomo):
        self.tomogram_names.add(name)
        self.tomograms[name] = tomo

class Tomogram():
    """A tomogram class, which can have a number of different dataframes associated with it."""
    # def __init__(self, name):
    #     self.name = name
    #     self.dataframes = {}
    #     self.dataframe_names = set()
    
    def __init__(self, name, groupnames, csvnames):
        self.name = name
        print(self.name)
        self.dataframes = {}
        self.dataframe_names = set()
        for index, csvname in enumerate(csvnames):
            print(groupnames[index])
            self.add_dataframe(groupnames[index], pd.read_csv(csvname))
    
    def add_dataframe(self, name, dataframe):
        self.dataframes[name] = dataframe
        self.dataframe_names.add(name)
    
    def has_key(self, key):
        return key in self.dataframe_names
    
    def __getitem__(self, key):
        if not key in self.dataframe_names:
            raise KeyError("Tomogram does not have dataframe {}".format(key))
        return self.dataframes[key]
    
    def __setitem__(self, key, value):
        self.dataframe_names.add(key)
        self.dataframes[key] = value

def weighted_median(values, weights):
    ''' compute the weighted median of values list. The 
weighted median is computed as follows:
    1- sort both lists (values and weights) based on values.
    2- select the 0.5 point from the weights and return the corresponding values as results
    e.g. values = [1, 3, 0] and weights=[0.1, 0.3, 0.6] assuming weights are probabilities.
    sorted values = [0, 1, 3] and corresponding sorted weights = [0.6,     0.1, 0.3] the 0.5 point on
    weight corresponds to the first item which is 0. so the weighted     median is 0.
    
    Function provided by Max Ghenis on stack overflow, CC-BY-SA'''

    #convert the weights into probabilities
    sum_weights = sum(weights)
    weights = np.array([(w*1.0)/sum_weights for w in weights])
    #sort values and weights based on values
    values = np.array(values)
    sorted_indices = np.argsort(values)
    values_sorted  = values[sorted_indices]
    weights_sorted = weights[sorted_indices]
    #select the median point
    it = np.nditer(weights_sorted, flags=['f_index'])
    accumulative_probability = 0
    median_index = -1
    while not it.finished:
        accumulative_probability += it[0]
        if accumulative_probability > 0.5:
            median_index = it.index
            return values_sorted[median_index]
        elif accumulative_probability == 0.5:
            median_index = it.index
            it.iternext()
            next_median_index = it.index
            return np.mean(values_sorted[[median_index, next_median_index]])
        it.iternext()

    return values_sorted[median_index]

def weighted_histogram_peak(values, weights, bins, bin_range):
    '''compute the peak of a histogram for a list of values.'''
    hist, bin_edges = np.histogram(values, bins=bins, range=bin_range, weights=weights, density=True)
    i = np.argmax(hist)
    bin_vals = (bin_edges[i], bin_edges[i+1])
    return(np.mean(bin_vals))


def statistics(datasets, basename,  condition_names, morph_names, test_type="median", filename="test.csv", figsize=(5,4), ylabel="Peak Value"):
    raw_filename = filename[:-10]+"rawstats.csv"
    with open(raw_filename, 'w') as raw_file:
        raw_file.write(basename+f" - {test_type}\n")
        for index, val in enumerate(datasets):
            raw_file.write(f"{condition_names[index]} {morph_names[index]},"+",".join(map(str, val))+"\n")
    with open(filename, 'w') as file:
        file.write(f"Base Experiment,Stat Type,Sample A Condition,Sample A Morph,Sample A Mean,Sample B Condition,Sample B Morph,Sample B Mean,U,P_U,T,P_T\n")
        for i, set_a in enumerate(datasets):
            for j, set_b in enumerate(datasets):
                if j<=i: 
                    continue
                try:
                    u,p_u = stats.mannwhitneyu(set_a, set_b)
                    t,p_t = stats.ttest_ind(set_a, set_b, equal_var=False)
                except e:
                    print(e)
                    u,p_u,t,p_t = -1,-1,-1,-1
                print(p_u, p_t)
                file.write(f"{basename},{test_type},{condition_names[i]},{morph_names[i]},{np.mean(set_a)},{condition_names[j]},{morph_names[j]},{np.mean(set_b)},{u},{p_u},{t},{p_t}\n")
    figure_filename = filename[:-3]+"svg"
    fig,ax=plt.subplots(figsize=figsize)
    ax.set_title(basename)
    ax.violinplot(datasets, showmeans=True)
    ax.set_xticks(range(1,len(datasets)+1))
    ax.set_xticklabels([condition_names[i]+"\n"+morph_names[i] for i in range(len(condition_names))])
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    fig.savefig(figure_filename)
    fig.savefig(figure_filename[:-3]+"png")


def histogram(data, areas, labels, title, xlabel, filename="hist.svg", bins=50, range=[0,50], figsize=(2.5,1.8), show=False, logx=False, vlines = True):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    if logx:
        bins = np.logspace(np.log10(range[0]),np.log10(range[1]), bins)
        ax.set_xscale("log")
    for index, value in enumerate(data):
        n,binset,_ = ax.hist(value, bins=bins, weights=areas[index], label=labels[index], ec=colors[index], fc=colors_light[index],histtype="stepfilled", density=False, range = range)
        if vlines:
            # print(binset)
            # print(n)
            # median_val = weighted_median(value, areas[index]) 
            # mean_val = np.average(value, weights=areas[index])
            # print(mean_val)
            # print(median_val)
            idx = n.argmax()
            print(idx)
            ax.axvline(binset[idx], linestyle="--", color=colors[index])
    ax.set_xlim(range)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative Area")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.legend()
    # plt.tight_layout()
    fig.savefig(filename, bbox_inches='tight')
    fig.savefig(filename[:-3]+"notitle.png", bbox_inches='tight', dpi=300)
    ax.set_title(title, loc='left')
    fig.savefig(filename[:-3]+"png", bbox_inches='tight', dpi=300)
    return True

def twod_histogram(data1, data2, areas, data1_label, data2_label, title, bins=None, range=[[0,100],[0,100]], filename="twod_hist.svg", figsize=(3,4), log=True):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title)
    ax.hist2d(data1, data2, bins=bins, range=range, density=True, weights=areas, norm=matplotlib.colors.LogNorm(), cmap=my_cmap)
    ax.set_xlabel(data1_label)
    ax.set_ylabel(data2_label)
    ax.set_xlim(range[0])
    ax.set_ylim(range[1])
    plt.tight_layout()
    fig.savefig(filename, bbox_inches='tight')
    fig.savefig(filename[:-3]+"png", bbox_inches='tight')
    return True

def barchart(bars, errors, labels, title, ylabel, filename="barchart.svg", figsize=(4,3), ymax=None, hline=None):
    print(len(bars))
    x = np.arange(len(bars))
    fig, ax = plt.subplots(figsize=figsize)
    barwidth=0.8
    ax.set_title(title)
    if hline:
        ax.axhline(hline, linestyle="--", alpha=0.5)
    ax.bar(x, bars, yerr=errors, tick_label=labels, width=barwidth, color ="0.7", edgecolor="0.3")
    ax.set_ylabel(ylabel)
    if ymax:
        ax.set_ylim(0,ymax)

    plt.tight_layout()
    fig.savefig(filename, bbox_inches='tight')
    fig.savefig(filename[:-3]+"png", bbox_inches='tight')
    return True



def double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel, legends = ["Untreated", "Treated"], filename="double_barchart.svg", figsize=(4,3)):
    x = np.arange(len(bars1))
    barwidth = 0.35
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title)
    ax.bar(x, bars1, yerr=errors1, width=barwidth, color = colors_light[0], edgecolor=colors[0], label=legends[0])
    ax.bar(x+barwidth, bars2, yerr=errors2, width=barwidth, color=colors_light[1], edgecolor=colors[1], label=legends[1])
    
    ax.set_ylabel(ylabel)
    ax.set_xticks(x+barwidth/2)
    ax.set_xticklabels(labels)
    ax.legend()
    plt.tight_layout()
    fig.savefig(filename, bbox_inches='tight')
    fig.savefig(filename[:-3]+"png", bbox_inches='tight')

if __name__ == "__main__":
    startdir = os.getcwd()
    output_dir = os.getcwd()+"/output/"
    print("Unpickling Treated")
    with open('treated.pkl', 'rb') as file:
        treated = pickle.load(file)
    print("Unpickling Untreated")
    with open('untreated.pkl', 'rb') as file:
        untreated = pickle.load(file)
    print("Unpickling TGGSK")
    with open('tggsk.pkl', 'rb') as file:
        tggsk = pickle.load(file)
    print(treated.tomograms["TT9"]["IMM"])
    print(tggsk.tomograms["GT1"]["IMM"])

    # # Prepare to write out graphs!
    os.chdir(output_dir)
    redo_charts = True
    if redo_charts:
        # OMM-IMM Distance
        
        labels = ["Veh", "Tg"]
        #Long
        title = "OMM-IMM distance (elongated)"
        distances = []
        areas = []
        distance_peak_sets = []
        conditions = ["Vehicle", "Tg", "Vehicle", "Tg"]
        morphologies = ["elongated", "elongated", "fragmented", "fragmented"]
        for file in untreated, treated:
            distance_set = []
            area_set = []
            peak_set = []
            for key,val in file.tomograms.items():
                # Skip Fragmented tomos!
                if key in ["UF3", "TE1"]:
                    continue
                if not key[1] == "F":
                    peak_set.append(weighted_histogram_peak(val["OMM"]["IMM_dist"], val["OMM"]["area"], 100, [5,25]))
                    distance_set.extend(val["OMM"]["IMM_dist"])
                    area_set.extend(val["OMM"]["area"])
            distance_peak_sets.append(peak_set)
            distances.append(distance_set)
            areas.append(list(np.array(area_set)/sum(area_set)))
        histogram(distances, areas, labels, title, xlabel="OMM-IMM Distance (nm)", filename="immdist_comparison_long.svg", bins=100, range=[5,25])
        #Short
        title = "OMM-IMM distance (fragmented)"
        distances = []
        areas = []
        for file in untreated, treated:
            distance_set = []
            area_set = []
            peak_set = []
            for key,val in file.tomograms.items():
                # Skip Non-fragmented tomos!
                if key in ["UF3", "TE1"]:
                    continue
                if key[1] == "F":
                    peak_set.append(weighted_histogram_peak(val["OMM"]["IMM_dist"], val["OMM"]["area"], 100, [5,25]))
                    distance_set.extend(val["OMM"]["IMM_dist"])
                    area_set.extend(val["OMM"]["area"])
            distance_peak_sets.append(peak_set)
            distances.append(distance_set)
            areas.append(list(np.array(area_set)/sum(area_set)))
        histogram(distances, areas, labels, title, xlabel="OMM-IMM Distance (nm)", filename="immdist_comparison_short.svg", bins=100, range=[5,25])
        statistics(distance_peak_sets, "OMM-IMM Distance", conditions, morphologies, test_type="Peaks", filename="imm_omm_distance_violin.csv", ylabel="Distance (nm)")
        #OMM-ER Distance
        
        # Long
        title = "OMM-ER distance (elongated)"
        distances = []
        areas = []
        for file in untreated, treated:
            distance_set = []
            area_set = []
            for key,val in file.tomograms.items():
                # Skip Non-fragmented tomos!
                if key in ["UF3", "TE1"]:
                    continue
                if not key[1] == "F":
                    if "ER_dist" in val["OMM"].columns.values:
                        distance_set.extend(val["OMM"]["ER_dist"])
                        area_set.extend(val["OMM"]["area"])

            distances.append(distance_set)
            areas.append(list(np.array(area_set)/sum(area_set)))
        histogram(distances, areas, labels, title, xlabel="OMM-ER Distance (nm)", filename="omm_er_comparison_long.svg", bins=100, range=[5,40])
        # Short
        title = "OMM-ER distance (fragmented)"
        distances = []
        areas = []
        for file in untreated, treated:
            distance_set = []
            area_set = []
            for key,val in file.tomograms.items():
                # Skip Non-fragmented tomos!
                if key in ["UF3", "TE1"]:
                    continue
                if key[1] == "F":
                    if "ER_dist" in val["OMM"].columns.values:
                        distance_set.extend(val["OMM"]["ER_dist"])
                        area_set.extend(val["OMM"]["area"])

            distances.append(distance_set)
            areas.append(list(np.array(area_set)/sum(area_set)))
        histogram(distances, areas, labels, title, xlabel="OMM-ER Distance (nm)", filename="omm_er_comparison_short.svg", bins=100, range=[5,40])
        # OMM/ER vs OMM/IMM 
        for file, label in ((untreated, "Vehicle"), (treated, "Tg"), (tggsk, "Tg+GSK")):
            # Short
            filename = "ER_OMM_IMM_2D_"+label+"_short.svg"
            title = f"{label} ER distance vs IMM distance for fragmented cells"
            imm_distances = []
            er_distances = []
            areas = []
            for key, val in file.tomograms.items():
                if key in ["UF3", "TE1"]:
                    continue
                if key[1] == "F":
                    if "ER_dist" in val["OMM"].columns.values:
                        er_distances.extend(val["OMM"]["ER_dist"])
                        imm_distances.extend(val["OMM"]["IMM_dist"])
                        areas.extend(val["OMM"]["area"])
            twod_histogram(imm_distances, er_distances, areas, "OMM-IMM Distance (nm)", "OMM-ER Distance (nm)", title, bins=(30,30), range=[[5,30],[5,50]], filename=filename, figsize=(4,4))
            # Long
            filename = "ER_OMM_IMM_2D_"+label+"_long.svg"
            title = f"{label} ER distance vs IMM distance for elongated cells"
            imm_distances = [] 
            er_distances = []
            areas = []
            for key, val in file.tomograms.items():
                if key in ["UF3", "TE1"]:
                    continue
                if not key[1] == "F":
                    if "ER_dist" in val["OMM"].columns.values:
                        er_distances.extend(val["OMM"]["ER_dist"])
                        imm_distances.extend(val["OMM"]["IMM_dist"])
                        areas.extend(val["OMM"]["area"])
            twod_histogram(imm_distances, er_distances, areas, "OMM-IMM Distance (nm)", "OMM-ER Distance (nm)", title, bins=(30,30), range=[[5,30],[5,50]], filename=filename, figsize=(4,4))

        # Absolute surface areas
        names = ["OMM", "IMM", "ER"]
        untreated_areas = {i: [] for i in names}
        treated_areas = {i: [] for i in names}
        tggsk_areas = {i: [] for i in names}
        for file, label, areas in ((untreated, "Vehicle", untreated_areas), (treated, "Tg", treated_areas)):
            for key, val in file.tomograms.items():
                if key in ["UF3", "TE1"]:
                    continue
                for tomoname in names:
                    if tomoname in val.dataframes:
                        tomo = val[tomoname]
                        areas[tomoname].append(sum(tomo["area"]))
                    else:
                        areas[tomoname].append(0)
            filename = f"absolute_areas_{label}.svg"
            title = f"Surface area per tomogram - {label}"
            labels = names
            bars = [np.mean(areas[i])/1000000 for i in names]
            print(bars)
            errors = [np.std(areas[i])/1000000 for i in names]
            print(errors)
            barchart(bars, errors, labels, title, ylabel="Area (um^2)", filename=filename)
        filename = "absolute_areas_comparison.svg"
        title = "Surface area per tomogram"
        labels = names
        legends = ["Vehicle", "Tg"]
        bars1 = [np.mean(untreated_areas[i])/1000000 for i in names] 
        errors1 = [np.std(untreated_areas[i])/1000000 for i in names]
        bars2 = [np.mean(treated_areas[i])/1000000 for i in names]
        errors2 = [np.std(treated_areas[i])/1000000 for i in names]
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)
        # Short
        names = ["OMM", "IMM", "ER"]
        untreated_areas = {i: [] for i in names}
        treated_areas = {i: [] for i in names}
        tggsk_areas = {i: [] for i in names}
        for file, label, areas in ((untreated, "Vehicle", untreated_areas), (treated, "Tg", treated_areas)):
            for key, val in file.tomograms.items():
                if key in ["UF3", "TE1"]:
                    continue
                if not key[1] == "F":
                    continue
                for tomoname in names:
                    if tomoname in val.dataframes:
                        tomo = val[tomoname]
                        areas[tomoname].append(sum(tomo["area"]))
                    else:
                        areas[tomoname].append(0)
            filename = f"absolute_areas_{label}_short.svg"
            title = f"Surface area per tomogram in fragmented cells - {label}"
            labels = names
            bars = [np.mean(areas[i])/1000000 for i in names]
            print(bars)
            errors = [np.std(areas[i])/1000000 for i in names]
            print(errors)
            barchart(bars, errors, labels, title, ylabel="Area (um^2)", filename=filename)
        filename = "absolute_areas_comparison_short.svg"
        title = "Surface area per tomogram in fragmented cells"
        labels = names
        legends = ["Untreated", "Treated"]
        bars1 = [np.mean(untreated_areas[i])/1000000 for i in names] 
        errors1 = [np.std(untreated_areas[i])/1000000 for i in names]
        bars2 = [np.mean(treated_areas[i])/1000000 for i in names]
        errors2 = [np.std(treated_areas[i])/1000000 for i in names]
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)
        # Long
        names = ["OMM", "IMM", "ER"]
        untreated_areas = {i: [] for i in names}
        treated_areas = {i: [] for i in names}
        tggsk_areas = {i: [] for i in names}
        for file, label, areas in ((untreated, "Vehicle", untreated_areas), (treated, "Tg", treated_areas), (tggsk, "Tg+GSK", tggsk_areas)):
            for key, val in file.tomograms.items():
                if key in ["UF3", "TE1"]:
                    continue
                if key[1] == "F":
                    continue
                for tomoname in names:
                    if tomoname in val.dataframes:
                        tomo = val[tomoname]
                        areas[tomoname].append(sum(tomo["area"]))
                    else:
                        areas[tomoname].append(0)
            filename = f"absolute_areas_{label}_long.svg"
            title = f"Surface area per tomogram in elongated cells- {label}"
            labels = names
            bars = [np.mean(areas[i])/1000000 for i in names]
            print(bars)
            errors = [np.std(areas[i])/1000000 for i in names]
            print(errors)
            barchart(bars, errors, labels, title, ylabel="Area (um^2)", filename=filename)
        filename = "absolute_areas_comparison_long.svg"
        title = "Surface area per tomogram in elongated cells"
        labels = names
        legends = ["Untreated", "Treated", "Cotreatment"]
        bars1 = [np.mean(untreated_areas[i])/1000000 for i in names] 
        errors1 = [np.std(untreated_areas[i])/1000000 for i in names]
        bars2 = [np.mean(treated_areas[i])/1000000 for i in names]
        errors2 = [np.std(treated_areas[i])/1000000 for i in names]
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)

        # ER Contact Site Area
        conditions = ["Vehicle", "Tg", "Vehicle", "Tg"]
        morphologies = ["elongated", "elongated", "fragmented", "fragmented"]
        # Long
        names = ["OMM", "ER"]
        untreated_contact_areas = {i: [] for i in names}
        untreated_contact_fractions = {i: [] for i in names}
        treated_contact_areas = {i: [] for i in names}
        treated_contact_fractions = {i: [] for i in names}
        tggsk_contact_areas = {i: [] for i in names}
        tggsk_contact_fractions = {i: [] for i in names}
        for file, label, contact_areas, contact_fractions in ((untreated, "Vehicle", untreated_contact_areas, untreated_contact_fractions), (treated, "Tg", treated_contact_areas, treated_contact_fractions)):
                for key, val in file.tomograms.items():
                    if key in ["UF3", "TE1", "TE4"]:
                        continue
                    if key[1] == "F":
                        continue
                    if not "ER" in val.dataframes:
                        continue
                    omm = val["OMM"]
                    er = val["ER"]
                    omm_area = sum(omm.area)
                    contact_areas["OMM"].append(omm_area)
                    omm_fraction = sum(omm[omm["ER_dist"] < 30].area)/omm_area
                    contact_fractions["OMM"].append(omm_fraction)
                    er_area = sum(er.area)
                    contact_areas["ER"].append(er_area)
                    er_fraction = sum(er[er["OMM_dist"] < 30].area)/er_area
                    contact_fractions["ER"].append(er_fraction)
        bars1 = [np.mean(untreated_contact_areas[i])/1000000 for i in names] 
        errors1 = [np.std(untreated_contact_areas[i])/1000000 for i in names]
        bars2 = [np.mean(treated_contact_areas[i])/1000000 for i in names]
        errors2 = [np.std(treated_contact_areas[i])/1000000 for i in names]
        labels = names
        legends = ["Veh", "Tg"]
        filename = "er_contact_site_areas_long.svg"
        title = "ER contact site area in elongated cells (<30nm)"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)
        bars1 = [np.mean(untreated_contact_fractions[i]) for i in names] 
        errors1 = [np.std(untreated_contact_fractions[i]) for i in names]
        bars2 = [np.mean(treated_contact_fractions[i]) for i in names]
        errors2 = [np.std(treated_contact_fractions[i]) for i in names]
        filename = "er_contact_site_fractions_long.svg"
        title = "Fraction of membrane in ER contact site in elongated cells"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Fraction of total area in contact site", legends=legends, filename=filename)
        omm_contact_areas = [untreated_contact_areas["OMM"], treated_contact_areas["OMM"]]
        omm_contact_fractions = [untreated_contact_fractions["OMM"], treated_contact_fractions["OMM"]]
        er_contact_areas = [untreated_contact_areas["ER"], treated_contact_areas["ER"]]
        er_contact_fractions = [untreated_contact_fractions["ER"], treated_contact_fractions["ER"]]
        # Short
        names = ["OMM", "ER"]
        untreated_contact_areas = {i: [] for i in names}
        untreated_contact_fractions = {i: [] for i in names}
        treated_contact_areas = {i: [] for i in names}
        treated_contact_fractions = {i: [] for i in names}
        tggsk_contact_areas = {i: [] for i in names}
        tggsk_contact_fractions = {i: [] for i in names}
        for file, label, contact_areas, contact_fractions in ((untreated, "Vehicle", untreated_contact_areas, untreated_contact_fractions), (treated, "Tg", treated_contact_areas, treated_contact_fractions)):
                for key, val in file.tomograms.items():
                    if key in ["UF3", "TE1", "TE4"]:
                        continue
                    if not key[1] == "F":
                        continue
                    if not "ER" in val.dataframes:
                        continue
                    omm = val["OMM"]
                    er = val["ER"]
                    omm_area = sum(omm.area)
                    contact_areas["OMM"].append(omm_area)
                    omm_fraction = sum(omm[omm["ER_dist"] < 30].area)/omm_area
                    contact_fractions["OMM"].append(omm_fraction)
                    er_area = sum(er.area)
                    contact_areas["ER"].append(er_area)
                    er_fraction = sum(er[er["OMM_dist"] < 30].area)/er_area
                    contact_fractions["ER"].append(er_fraction)
        bars1 = [np.mean(untreated_contact_areas[i])/1000000 for i in names] 
        errors1 = [np.std(untreated_contact_areas[i])/1000000 for i in names]
        bars2 = [np.mean(treated_contact_areas[i])/1000000 for i in names]
        errors2 = [np.std(treated_contact_areas[i])/1000000 for i in names]
        labels = names
        legends = ["Veh", "Tg"]
        filename = "er_contact_site_areas_short.svg"
        title = "ER contact site area in fragmented cells(<30nm)"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)
        bars1 = [np.mean(untreated_contact_fractions[i]) for i in names] 
        errors1 = [np.std(untreated_contact_fractions[i]) for i in names]
        bars2 = [np.mean(treated_contact_fractions[i]) for i in names]
        errors2 = [np.std(treated_contact_fractions[i]) for i in names]
        filename = "er_contact_site_fractions_short.svg"
        title = "Fraction of membrane in ER contact site in fragmented cells"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Fraction of total area in contact site", legends=legends, filename=filename)
        omm_contact_areas.extend([untreated_contact_areas["OMM"], treated_contact_areas["OMM"]])
        omm_contact_fractions.extend([untreated_contact_fractions["OMM"], treated_contact_fractions["OMM"]])
        er_contact_areas.extend([untreated_contact_areas["ER"], treated_contact_areas["ER"]])
        er_contact_fractions.extend([untreated_contact_fractions["ER"], treated_contact_fractions["ER"]])
        statistics(omm_contact_areas, "OMM-ER Contact Area", conditions,morphologies, test_type="Peaks", filename="OMM_ER_Contact_OMM_violin.csv", ylabel="OMM Area (nm^2)")
        statistics(omm_contact_fractions, "Fraction of OMM contacting ER", conditions,morphologies, test_type="Peaks", filename="OMM_ER_Fraction_OMM_violin.csv", ylabel="Fraction of OMM Area")
        statistics(er_contact_areas, "OMM-ER Contact Area", conditions,morphologies, test_type="Peaks", filename="OMM_ER_Contact_ER_violin.csv", ylabel="ER Area (nm^2)")
        statistics(er_contact_fractions, "Fraction of ER contacting OMM", conditions,morphologies, test_type="Peaks", filename="OMM_ER_Fraction_ER_violin.csv", ylabel="Fraction of ER Area")
                
        # 2D Curvature
        # Long
        for file, label in ((untreated, "Vehicle"), (treated, "Tg"), (tggsk, "Tg+GSK")):
            curvatures = []
            distances = []
            areas = []
            for key, val in file.tomograms.items():
                if key in ["UF3", "TE1"]:
                    continue
                if key[1] == "F":
                    continue
                imm = val["IMM"]
                curvatures.extend(imm["curvedness_VV"])
                distances.extend(imm["OMM_dist"])
                areas.extend(imm["area"])
            filename = f"curvedness_vs_distance_{label}_long.svg"
            twod_histogram(distances, curvatures, areas, "IMM-OMM Distance (nm)", "Curvedness (1/nm)", title="Curvedness vs Distance in the IMM in elongated cells", filename = filename, bins=(50,100), range=[[0,50],[0,0.1]], figsize=(4,4), log=True)

        #localized curvature
        # Curvature of cristae
        # Long
        for feature in ['IMM', 'IBM', 'crista', 'junction']:
            basename = f"{feature[0].title()+feature[1:]} curvedness"
            filename = f"curvedness_{feature}_violin.csv"
            curve_peaks = []
            conditions = []
            morphologies = []
            for run in ['elongated', 'fragmented']:
                curvature_list = []
                log_curvature_list = []
                area_list = []
                labels = ["Veh", "Tg"]
                counts = []
                for file, label in ((untreated, "Vehicle"), (treated, "Tg")):
                    peaks = []
                    morphologies.append(run)
                    conditions.append(label)
                    curvatures = []
                    # distances = []
                    areas = []
                    count = 0
                    for key, val in file.tomograms.items():
                        if key in ["UF3", "TE1"]:
                            continue
                        if run == 'elongated':
                            if key[1] == "F":
                                continue
                        elif run == 'fragmented':
                            if not key[1] == 'F':
                                continue

                        count += 1
                        imm = val["IMM"]
                        # imm = imm[imm["OMM_dist"]>23]
                        if feature == 'IBM':
                            imm = imm[imm["OMM_dist"]<junction_minimum]
                        elif feature == 'crista':
                            imm = imm[imm["OMM_dist"]>junction_maximum]
                        elif feature == 'junction':
                            imm = imm[imm["OMM_dist"]<junction_maximum]
                            imm = imm[imm["OMM_dist"]>junction_minimum]
                        curvatures.extend(imm["curvedness_VV"])
                        peaks.append(weighted_histogram_peak(imm["curvedness_VV"], imm["area"], 100, [0,0.1]))
                        # distances.extend(imm["OMM_dist"])
                        areas.extend(imm["area"])
                    # filename = f"curvedness_{label}_long.svg"
                    # twod_histogram(distances, curvatures, areas, "IMM-OMM Distance (nm)", "Curvedness (1/nm)", title="Curvedness vs Distance in the IMM", filename = filename, bins=(50,100), range=[[0,50],[0,0.1]], figsize=(4,4), log=True)
                    curve_peaks.append(peaks)
                    areas = np.array(areas)
                    area_list.append(areas/sum(areas))
                    counts.append(count)
                    curvature_list.append(curvatures)
                    log_curvature_list.append(np.log(np.array(curvatures)))
                histogram(curvature_list, area_list, labels, title=f"{feature[0].title()+feature[1:]} curvedness ({run})", xlabel="Curvedness (1/nm)", bins=50, range=[5e-4,2e-1], filename=f"{feature}_curvedness_xlog_{run}.svg",logx=True)
                histogram(curvature_list, area_list, labels, title=f"{feature[0].title()+feature[1:]} curvedness ({run})", xlabel="Curvedness (1/nm)", bins=50, range=[0,0.1], filename=f"{feature}_curvedness_{run}.svg",logx=False)
                statistics(curve_peaks, basename, conditions,morphologies, test_type="Peaks", filename=filename, ylabel="Curvedness (1/nm)" )
                # av0, std0 = weighted_avg_and_std(log_curvature_list[0], area_list[0])
                # av1,std1 = weighted_avg_and_std(log_curvature_list[1], area_list[1])
                # print(counts)
                # print(stats.ttest_ind_from_stats(av0, std0, counts[0], av1, std1, counts[1]))
                # print(stats.ks_2samp(curvature_list[0], curvature_list[1]))

        # Relative` surface area
        labels = ["Veh", "Tg"]
        for run in ['elongated', 'fragmented', 'all']:
            area_set = []
            std_set = []
            for file, label in ((untreated, "Vehicle"), (treated, "Tg")): 
                areas = []
                for key, val in file.tomograms.items():
                    if key in ["UF3", "TE1"]:
                        continue
                    if run == 'elongated':
                        if key[1] == "F":
                            continue
                    elif run == 'fragmented':
                        if not key[1] == 'F':
                            continue
                    print(key)
                    omm_area = sum(val["OMM"]["area"])

                    imm = val["IMM"]
                    imm_area = sum(imm['area'])/omm_area
                    ibm = imm[imm["OMM_dist"]<junction_minimum]
                    ibm_area = sum(ibm['area'])/omm_area
                    cristae = imm[imm["OMM_dist"] > junction_maximum]
                    cristae_area = sum(cristae['area'])/omm_area
                    junctions = imm[imm["OMM_dist"] < junction_maximum]
                    junctions = junctions[junctions["OMM_dist"]>junction_minimum]
                    junctions_area = sum(junctions['area'])/omm_area

                    areas.append((ibm_area, junctions_area, cristae_area, imm_area))
                
                areas = np.array(areas)
                mean = np.mean(areas, axis=0)
                std = np.std(areas, axis=0)
                area_set.append(mean)
                std_set.append(std)
                barchart(mean, std, labels=["IBM", "Junctions", "Cristae", "All"], title=f'Relative Surface Areas of {label} in {run} cells', ylabel="Surface Area Relative to OMM", filename=f"{label}_relative_areas_{run}.svg", figsize=(4,3))
            average = [area_set[1][i]/area_set[0][i] for i in range(4)]
            std = [np.sqrt((std_set[0][i]/area_set[0][i])**2+(std_set[1][i]/area_set[1][i])**2) for i in range(4)]
            barchart(average, std, labels=["IBM", "Junctions", "Cristae", "All"], title = f"Treated Surface Area Relative to Untreated in {run} cells", ylabel="Relative Surface Area", filename=f"ratio_relative_areas_{run}.svg", figsize=(4,3), ymax = 2.5, hline = 1)
            # average = [area_set[2][i]/area_set[0][i] for i in range(4)]
            # std = [np.sqrt((std_set[0][i]/area_set[0][i])**2+(std_set[2][i]/area_set[2][i])**2) for i in range(4)]
            # barchart(average, std, labels=["IBM", "Junctions", "Cristae", "All"], title = f"Cotreatment Surface Area Relative to Untreated in {run} cells", ylabel="Relative Surface Area", filename=f"ratio_relative_areas_cotreatment_{run}.svg", figsize=(4,3), ymax = 2.5, hline = 1)

            
            
            double_barchart(area_set[0], area_set[1], std_set[0], std_set[1], labels=["IBM", "Junctions", "Cristae", "All"], title="Relative Surface Area", ylabel="Surface Area Relative to OMM", legends = labels, filename=f"double_barchart_relative_areas_{run}.svg", figsize=(4,3))

        # Intra-IMM Spacing
        #  2D
        labels = ["Veh", "Tg"]
        for feature in ["crista", "junction", "IBM", "IMM"]:
            basename_close = f"Intra-{feature} Spacing"
            basename_far = f"Inter-{feature} Spacing"
            filename_close = f"intra{feature}_spacing_violin.csv"
            filename_far = f"inter{feature}_spacing_violin.csv"
            close_peaks = []
            far_medians = []
            far_peaks = []
            conditions = []
            morphs = []
            for run in ['elongated', 'fragmented']: #, 'all'
                close_list = []
                far_list = []
                area_list = []
                for file, label in ((untreated, "Vehicle"), (treated, "Tg")):
                    conditions.append(label)
                    morphs.append(run)
                    distances_short = []
                    areas = []
                    distances_long = []
                    peak_list_close = []
                    peak_list_far = []
                    for key, val in file.tomograms.items():
                        if key in ["UF3", "TE1"]:
                            continue
                        if run == 'elongated':
                            if key[1] == "F":
                                continue
                        elif run == 'fragmented':
                            if not key[1] == 'F':
                                continue
                        imm = val["IMM"]
                        # imm = imm[imm["OMM_dist"]>23]
                        if feature == 'IBM':
                            imm = imm[imm["OMM_dist"]<junction_minimum]
                        elif feature == 'crista':
                            imm = imm[imm["OMM_dist"]>junction_maximum]
                        elif feature == 'junction':
                            imm = imm[imm["OMM_dist"]<junction_maximum]
                            imm = imm[imm["OMM_dist"]>junction_minimum]
                        imm = imm[imm["self_id_min"] != -1]
                        peak_list_close.append(weighted_histogram_peak(imm["self_dist_min"], imm["area"], 100, [10,40]))
                        peak_list_far.append(weighted_histogram_peak(imm["self_dist_far"], imm["area"], 100, [20,200]))
                        distances_long.extend(imm["self_dist_far"])
                        distances_short.extend(imm["self_dist_min"])
                        areas.extend(imm["area"])
                    close_peaks.append(peak_list_close)
                    far_peaks.append(peak_list_far)
                    close_list.append(distances_short)
                    far_list.append(distances_long)
                    areas = np.array(areas)
                    area_list.append(areas/sum(areas))
                
                histogram(close_list, area_list, labels, title=f"Intra-{feature} spacing ({run})", xlabel="Distance (nm)", bins=50, range=[10,40], filename=f"{feature}_close_intra_{run}.svg", vlines=True)
                histogram(far_list, area_list, labels, title=f"Inter-{feature} spacing ({run})", xlabel="Distance (nm)", bins=50, range=[20,100], filename=f"{feature}_far_intra_{run}.svg", vlines=True)
            statistics(close_peaks, basename_close, conditions, morphs, test_type="Peaks", ylabel="Distance (nm)", filename=filename_close)
            statistics(far_peaks, basename_far, conditions, morphs, test_type="Peaks", ylabel="Distance (nm)", filename=filename_far)

        for run in ['elongated', 'fragmented', 'all']:
            curvature_list = []
            log_curvature_list = []
            area_list = []
            labels = ["Veh", "Tg"]
            counts = []
            for file, label in ((untreated, "Vehicle"), (treated, "Tg")):
                curvatures = []
                # distances = []
                areas = []
                count = 0
                for key, val in file.tomograms.items():
                    if key in ["UF3", "TE1"]:
                        continue
                    if run == 'elongated':
                        if key[1] == "F":
                            continue
                    elif run == 'fragmented':
                        if not key[1] == 'F':
                            continue

                    count += 1
                    omm = val["OMM"]

                    curvatures.extend(omm["curvedness_VV"])
                    # distances.extend(imm["OMM_dist"])
                    areas.extend(omm["area"])
                # filename = f"curvedness_{label}_long.svg"
                # twod_histogram(distances, curvatures, areas, "IMM-OMM Distance (nm)", "Curvedness (1/nm)", title="Curvedness vs Distance in the IMM", filename = filename, bins=(50,100), range=[[0,50],[0,0.1]], figsize=(4,4), log=True)
                areas = np.array(areas)
                area_list.append(areas/sum(areas))
                counts.append(count)
                curvature_list.append(curvatures)
                log_curvature_list.append(np.log(np.array(curvatures)))
            histogram(curvature_list, area_list, labels, title=f"OMM Curvedness ({run})", xlabel="Curvedness (1/nm)", bins=50, range=[5e-4,2e-1], filename=f"OMM_curvedness_xlog_{run}.svg",logx=True)
            histogram(curvature_list, area_list, labels, title=f"OMM Curvedness ({run})", xlabel="Curvedness (1/nm)", bins=50, range=[0,0.1], filename=f"OMM_curvedness_{run}.svg",logx=False)


    ## Verticality
    print("Verticality")
    for feature in ['IMM', 'IBM', 'cristae', 'junctions']:
        vert_medians = []
        vert_means = []
        vert_stds = []
        conditions = []
        morphs = []
        # filename_means = f"{feature}_verticality_means_violin.csv"
        # filename_medians = f"{feature}_verticality_medians_violin.csv"
        filename_stds = f"{feature}_verticality_stds_violin.csv"
        basename = f"Orientation of {feature} relative to the growth plane"
        for run in ['elongated', 'fragmented']:
            angle_list = []
            area_list = []
            labels = ["Veh", "Tg"]
            counts = []
            for file, label in ((untreated, "Vehicle"), (treated, "Tg")):
                conditions.append(label)
                morphs.append(run)
                angles = []
                # distances = []
                # medians = []
                # means = []
                stds = []
                areas = []
                count = 0
                for key, val in file.tomograms.items():
                    if key in ["UF3", "TE1"]:
                        continue
                    if run == 'elongated':
                        if key[1] == "F":
                            continue
                    elif run == 'fragmented':
                        if not key[1] == 'F':
                            continue

                    count += 1
                    imm = val["IMM"]
                    # imm = imm[imm["OMM_dist"]>23]
                    if feature == 'IBM':
                        imm = imm[imm["OMM_dist"]<junction_minimum]
                    elif feature == 'cristae':
                        imm = imm[imm["OMM_dist"]>junction_maximum]
                    elif feature == 'junctions':
                        imm = imm[imm["OMM_dist"]<junction_maximum]
                        imm = imm[imm["OMM_dist"]>junction_minimum]
                    # curvatures.extend(imm["curvedness_VV"])
                    z = imm["n_v_z"]
                    vert = 90-np.abs(np.arccos(z)*180/np.pi-90)
                    angles.extend(vert)
                    areas.extend(imm["area"])
                    # medians.append(weighted_median(vert, imm["area"]))
                    # means.append(np.average(vert, weights=imm["area"]))
                    stds.append(np.std(vert))
                # filename = f"curvedness_{label}_long.svg"
                # twod_histogram(distances, curvatures, areas, "IMM-OMM Distance (nm)", "Curvedness (1/nm)", title="Curvedness vs Distance in the IMM", filename = filename, bins=(50,100), range=[[0,50],[0,0.1]], figsize=(4,4), log=True)
                areas = np.array(areas)
                area_list.append(areas/sum(areas))
                counts.append(count)
                angle_list.append(angles)
                # vert_means.append(means)
                # vert_medians.append(medians)
                vert_stds.append(stds)
            histogram(angle_list, area_list, labels, title=f"{feature[0].title()+feature[1:]}-growth plane angle ({run})", xlabel="Angle (º)", bins=60, range=[0,90], filename=f"{feature}_verticality_{run}.svg")
        statistics(vert_stds, basename, conditions, morphs, test_type="STD", ylabel="Standard Deviation of Distance (nm)", filename=filename_stds)
        # statistics(vert_means, basename, conditions, morphs, test_type="Means", ylabel="Distance (nm)", filename=filename_means)
        # statistics(vert_medians, basename, conditions, morphs, test_type="Medians", ylabel="Distance (nm)", filename=filename_medians)
        

    ## Relative angle  
    print("Relative Angle")
    for feature in ['IMM', 'IBM', 'cristae', 'junctions']:
        ang_stds = []
        morphologies = []
        conditions = []
        filename=f"rel_angle_{feature}_violin.csv"
        basename=f"Angle between {feature} and OMM"
        for run in ['elongated', 'fragmented']:
            angle_list = []
            area_list = []
            labels = ["Veh", "Tg"]
            counts = []
            for file, label in ((untreated, "Vehicle"), (treated, "Tg")):
                angles = []
                # distances = []
                areas = []
                count = 0
                stds = []
                morphologies.append(run)
                conditions.append(label)
                for key, val in file.tomograms.items():
                    if key in ["UF3", "TE1"]:
                        continue
                    if run == 'elongated':
                        if key[1] == "F":
                            continue
                    elif run == 'fragmented':
                        if not key[1] == 'F':
                            continue

                    count += 1
                    imm = val["IMM"]

                    omm = val["OMM"]
                    # imm = imm[imm["OMM_dist"]>23]
                    if feature == 'IBM':
                        imm = imm[imm["OMM_dist"]<junction_minimum]
                    elif feature == 'cristae':
                        imm = imm[imm["OMM_dist"]>junction_maximum]
                    elif feature == 'junctions':
                        imm = imm[imm["OMM_dist"]<junction_maximum]
                        imm = imm[imm["OMM_dist"]>junction_minimum]
                    # curvatures.extend(imm["curvedness_VV"])
                    a_x = imm["n_v_x"]
                    a_y = imm["n_v_y"]
                    a_z = imm["n_v_z"]
                    b_x = omm["n_v_x"][imm["OMM_neighbor_index"]]
                    b_y = omm["n_v_y"][imm["OMM_neighbor_index"]]
                    b_z = omm["n_v_z"][imm["OMM_neighbor_index"]]
                    a_vec = np.array([a_x,a_y,a_z]).transpose()
                    b_vec = np.array([b_x,b_y,b_z]).transpose()
                    ang = [np.arccos(np.abs(np.dot(a_vec[i], b_vec[i])))*180/np.pi for i in range(len(imm))]
                    angles.extend(ang)
                    areas.extend(imm["area"])
                    stds.append(np.std(ang))
                # filename = f"curvedness_{label}_long.svg"
                # twod_histogram(distances, curvatures, areas, "IMM-OMM Distance (nm)", "Curvedness (1/nm)", title="Curvedness vs Distance in the IMM", filename = filename, bins=(50,100), range=[[0,50],[0,0.1]], figsize=(4,4), log=True)
                areas = np.array(areas)
                area_list.append(areas/sum(areas))
                counts.append(count)
                angle_list.append(angles)
                ang_stds.append(stds)
            histogram(angle_list, area_list, labels, title=f"{feature[0].title()+feature[1:]}-OMM angle ({run})", xlabel="Angle (º)", bins=90, range=[0,90], filename=f"{feature}_angle_from_OMM_{run}.svg")
        statistics(ang_stds, basename, conditions, morphs, test_type="STD", ylabel="Standard Deviation of Angle (º)", filename=filename)


    # if not os.path.isdir(output_dir):
    #     os.mkdir(output_dir)
    # experiments = ["TgGSK", "Treated", "Untreated"]
    # # # # # # experiments = ["Treated"]
    # folders = ["/gpfs/group/grotjahn/bbarad/Final_Dataset/"+i+"/" for i in experiments]
    # basenames = ["OMM", "IMM", "ER"]
    # file_extension = ".AVV_rh12.csv"   
    # experiment_classes = {}
    # for index, experiment in enumerate(experiments):
    #     os.chdir(folders[index])
    #     print(experiment)
    #     tomo_names = glob.glob("*_labels.mrc")
    #     tomo_names = [i.split("_")[0] for i in tomo_names]
    #     experiment_classes[experiment] = Experiment(experiment)
    #     experiment_classes[experiment].add_tomograms(tomo_names, basenames, folders[index], file_extension)
    # os.chdir(startdir)
    # print("Treated")
    # treated = experiment_classes["Treated"]
    # with open("treated.pkl", 'wb') as file: 
    #     pickle.dump(treated, file)
    # print("Untreated")
    # untreated = experiment_classes["Untreated"]
    # with open("untreated.pkl", 'wb') as file: 
    #     pickle.dump(untreated, file)
    # print("TgGSK")
    # tggsk = experiment_classes["TgGSK"]
    # with open("tggsk.pkl", 'wb') as file: 
    #     pickle.dump(tggsk, file)
    # # print(tggsk)
    # print(tggsk.tomogram_names)
    # print(tggsk["GF1"])
    # print(tggsk["GF1"]["OMM"])

    # Histograms
    # Distance comparison
    

    
 




