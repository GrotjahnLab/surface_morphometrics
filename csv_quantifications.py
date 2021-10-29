import copy
import glob
import os
import pickle

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats

purple = [.5, 0, .5, 0.9]
purple_light = [.5,0,.5, 0.3]
green = [0, 1, 0, 0.9]
green_light = [0,1,0,.3]
blue = [0,.3,1, .9]
blue_light = [0,.3,1,.3]

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

def histogram(data, areas, labels, title, xlabel, filename="hist.svg", bins=50, range=[0,50], figsize=(4,3), show=False, logx=False):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title)
    if logx:
        bins = np.logspace(np.log10(range[0]),np.log10(range[1]), bins)
        ax.set_xscale("log")
    for index, value in enumerate(data):
        ax.hist(value, bins=bins, weights=areas[index], label=labels[index], ec=colors[index], fc=colors_light[index],histtype="stepfilled", density=False, range = range)
    ax.set_xlim(range)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative Area")
    ax.legend()
    plt.tight_layout()
    fig.savefig(filename, bbox_inches='tight')
    fig.savefig(filename[:-3]+"png", bbox_inches='tight')
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

def barchart(bars, errors, labels, title, ylabel, filename="barchart.svg", figsize=(4,3)):
    print(len(bars))
    x = np.arange(len(bars))
    fig, ax = plt.subplots(figsize=figsize)
    barwidth=0.8
    ax.set_title(title)
    ax.bar(x, bars, yerr=errors, tick_label=labels, width=barwidth, color ="0.7", edgecolor="0.3")
    ax.set_ylabel(ylabel)
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
    # Prepare to write out graphs!
    os.chdir(output_dir)
    redo_charts = False
    if redo_charts:
        # OMM-IMM Distance
        title = "OMM-IMM Distances"
        labels = ["Untreated", "Treated"]
        #Long
        distances = []
        areas = []
        for file in untreated, treated:
            distance_set = []
            area_set = []
            for key,val in file.tomograms.items():
                # Skip Fragmented tomos!
                if not key[1] == "F":
                    distance_set.extend(val["OMM"]["IMM_dist"])
                    area_set.extend(val["OMM"]["area"])

            distances.append(distance_set)
            areas.append(area_set)
        histogram(distances, areas, labels, title, xlabel="OMM-IMM Distance (nm)", filename="immdist_comparison_long.svg", bins=100, range=[5,25])
        #Short
        distances = []
        areas = []
        for file in untreated, treated:
            distance_set = []
            area_set = []
            for key,val in file.tomograms.items():
                # Skip Non-fragmented tomos!
                if key[1] == "F":
                    distance_set.extend(val["OMM"]["IMM_dist"])
                    area_set.extend(val["OMM"]["area"])

            distances.append(distance_set)
            areas.append(area_set)
        histogram(distances, areas, labels, title, xlabel="OMM-IMM Distance (nm)", filename="immdist_comparison_short.svg", bins=100, range=[5,25])
        
        #OMM-ER Distance
        title = "OMM-ER Distances"
        # Short
        distances = []
        areas = []
        for file in untreated, treated:
            distance_set = []
            area_set = []
            for key,val in file.tomograms.items():
                # Skip Non-fragmented tomos!
                if key[1] == "F":
                    if "ER_dist" in val["OMM"].columns.values:
                        distance_set.extend(val["OMM"]["ER_dist"])
                        area_set.extend(val["OMM"]["area"])

            distances.append(distance_set)
            areas.append(area_set)
        histogram(distances, areas, labels, title, xlabel="OMM-ER Distance (nm)", filename="omm_er_comparison_short.svg", bins=100, range=[5,40])
        # Long
        distances = []
        areas = []
        for file in untreated, treated:
            distance_set = []
            area_set = []
            for key,val in file.tomograms.items():
                # Skip Non-fragmented tomos!
                if not key[1] == "F":
                    if "ER_dist" in val["OMM"].columns.values:
                        distance_set.extend(val["OMM"]["ER_dist"])
                        area_set.extend(val["OMM"]["area"])

            distances.append(distance_set)
            areas.append(area_set)
        histogram(distances, areas, labels, title, xlabel="OMM-ER Distance (nm)", filename="omm_er_comparison_long.svg", bins=100, range=[5,40])
        
        # OMM/ER vs OMM/IMM 
        for file, label in ((untreated, "Untreated"), (treated, "Treated")):
            # Short
            filename = "ER_OMM_IMM_2D_"+label+"_short.svg"
            title = f"{label} ER distance vs IMM distance"
            imm_distances = []
            er_distances = []
            areas = []
            for key, val in file.tomograms.items():
                if key[1] == "F":
                    if "ER_dist" in val["OMM"].columns.values:
                        er_distances.extend(val["OMM"]["ER_dist"])
                        imm_distances.extend(val["OMM"]["IMM_dist"])
                        areas.extend(val["OMM"]["area"])
            twod_histogram(imm_distances, er_distances, areas, "OMM-IMM Distance (nm)", "OMM-ER Distance (nm)", title, bins=(30,30), range=[[5,30],[5,50]], filename=filename, figsize=(4,4))
            # Long
            filename = "ER_OMM_IMM_2D_"+label+"_long.svg"
            title = f"{label} ER distance vs IMM distance"
            imm_distances = [] 
            er_distances = []
            areas = []
            for key, val in file.tomograms.items():
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
        for file, label, areas in ((untreated, "Untreated", untreated_areas), (treated, "Treated", treated_areas)):
            for key, val in file.tomograms.items():
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
        legends = ["Untreated", "Treated"]
        bars1 = [np.mean(untreated_areas[i])/1000000 for i in names] 
        errors1 = [np.std(untreated_areas[i])/1000000 for i in names]
        bars2 = [np.mean(treated_areas[i])/1000000 for i in names]
        errors2 = [np.std(treated_areas[i])/1000000 for i in names]
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)
        # Short
        names = ["OMM", "IMM", "ER"]
        untreated_areas = {i: [] for i in names}
        treated_areas = {i: [] for i in names}
        for file, label, areas in ((untreated, "Untreated", untreated_areas), (treated, "Treated", treated_areas)):
            for key, val in file.tomograms.items():
                if not key[1] == "F":
                    continue
                for tomoname in names:
                    if tomoname in val.dataframes:
                        tomo = val[tomoname]
                        areas[tomoname].append(sum(tomo["area"]))
                    else:
                        areas[tomoname].append(0)
            filename = f"absolute_areas_{label}_short.svg"
            title = f"Surface area per tomogram - {label}"
            labels = names
            bars = [np.mean(areas[i])/1000000 for i in names]
            print(bars)
            errors = [np.std(areas[i])/1000000 for i in names]
            print(errors)
            barchart(bars, errors, labels, title, ylabel="Area (um^2)", filename=filename)
        filename = "absolute_areas_comparison_short.svg"
        title = "Surface area per tomogram"
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
        for file, label, areas in ((untreated, "Untreated", untreated_areas), (treated, "Treated", treated_areas)):
            for key, val in file.tomograms.items():
                if key[1] == "F":
                    continue
                for tomoname in names:
                    if tomoname in val.dataframes:
                        tomo = val[tomoname]
                        areas[tomoname].append(sum(tomo["area"]))
                    else:
                        areas[tomoname].append(0)
            filename = f"absolute_areas_{label}_long.svg"
            title = f"Surface area per tomogram - {label}"
            labels = names
            bars = [np.mean(areas[i])/1000000 for i in names]
            print(bars)
            errors = [np.std(areas[i])/1000000 for i in names]
            print(errors)
            barchart(bars, errors, labels, title, ylabel="Area (um^2)", filename=filename)
        filename = "absolute_areas_comparison_long.svg"
        title = "Surface area per tomogram"
        labels = names
        legends = ["Untreated", "Treated"]
        bars1 = [np.mean(untreated_areas[i])/1000000 for i in names] 
        errors1 = [np.std(untreated_areas[i])/1000000 for i in names]
        bars2 = [np.mean(treated_areas[i])/1000000 for i in names]
        errors2 = [np.std(treated_areas[i])/1000000 for i in names]
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)

        # ER Contact Site Area
        # Long
        names = ["OMM", "ER"]
        untreated_contact_areas = {i: [] for i in names}
        untreated_contact_fractions = {i: [] for i in names}
        treated_contact_areas = {i: [] for i in names}
        treated_contact_fractions = {i: [] for i in names}
        for file, label, contact_areas, contact_fractions in ((untreated, "Untreated", untreated_contact_areas, untreated_contact_fractions), (treated, "Treated", treated_contact_areas, treated_contact_fractions)):
                for key, val in file.tomograms.items():
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
        legends = ["Untreated", "Treated"]
        filename = "er_contact_site_areas_long.svg"
        title = "ER contact site area (<30nm)"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)
        bars1 = [np.mean(untreated_contact_fractions[i]) for i in names] 
        errors1 = [np.std(untreated_contact_fractions[i]) for i in names]
        bars2 = [np.mean(treated_contact_fractions[i]) for i in names]
        errors2 = [np.std(treated_contact_fractions[i]) for i in names]
        filename = "er_contact_site_fractions_long.svg"
        title = "Fraction of membrane in ER contact site"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Fraction of total area in contact site", legends=legends, filename=filename)

        # Short
        names = ["OMM", "ER"]
        untreated_contact_areas = {i: [] for i in names}
        untreated_contact_fractions = {i: [] for i in names}
        treated_contact_areas = {i: [] for i in names}
        treated_contact_fractions = {i: [] for i in names}
        for file, label, contact_areas, contact_fractions in ((untreated, "Untreated", untreated_contact_areas, untreated_contact_fractions), (treated, "Treated", treated_contact_areas, treated_contact_fractions)):
                for key, val in file.tomograms.items():
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
        legends = ["Untreated", "Treated"]
        filename = "er_contact_site_areas_short.svg"
        title = "ER contact site area (<30nm)"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Area (um^2)", legends=legends, filename=filename)
        bars1 = [np.mean(untreated_contact_fractions[i]) for i in names] 
        errors1 = [np.std(untreated_contact_fractions[i]) for i in names]
        bars2 = [np.mean(treated_contact_fractions[i]) for i in names]
        errors2 = [np.std(treated_contact_fractions[i]) for i in names]
        filename = "er_contact_site_fractions_short.svg"
        title = "Fraction of membrane in ER contact site"
        double_barchart(bars1, bars2, errors1, errors2, labels, title, ylabel="Fraction of total area in contact site", legends=legends, filename=filename)

                
        # 2D Curvature
        # Long
        for file, label in ((untreated, "Untreated"), (treated, "Treated")):
            curvatures = []
            distances = []
            areas = []
            for key, val in file.tomograms.items():
                if key[1] == "F":
                    continue
                imm = val["IMM"]
                curvatures.extend(imm["curvedness_VV"])
                distances.extend(imm["OMM_dist"])
                areas.extend(imm["area"])
            filename = f"curvedness_vs_distance_{label}_long.svg"
            twod_histogram(distances, curvatures, areas, "IMM-OMM Distance (nm)", "Curvedness (1/nm)", title="Curvedness vs Distance in the IMM", filename = filename, bins=(50,100), range=[[0,50],[0,0.1]], figsize=(4,4), log=True)

    #localized curvature
    # Curvature of cristae
    # Long
    for feature in ['imm', 'ibm', 'cristae', 'junctions']:
        for run in ['long', 'short', 'all']:
            curvature_list = []
            log_curvature_list = []
            area_list = []
            labels = ["Untreated", "Treated", "Cotreated"]
            counts = []
            for file, label in ((untreated, "Untreated"), (treated, "Treated"), (tggsk,"Cotreated")):
                curvatures = []
                # distances = []
                areas = []
                count = 0
                for key, val in file.tomograms.items():
                    if run == 'long':
                        if key[1] == "F":
                            continue
                    elif run == 'short':
                        if not key[1] == 'F':
                            continue

                    count += 1
                    imm = val["IMM"]
                    # imm = imm[imm["OMM_dist"]>23]
                    if feature == 'ibm':
                        imm = imm[imm["OMM_dist"]<16]
                    elif feature == 'cristae':
                        imm = imm[imm["OMM_dist"]>23]
                    elif feature == 'junctions':
                        imm = imm[imm["OMM_dist"]<23]
                        imm = imm[imm["OMM_dist"]>16]
                    curvatures.extend(imm["curvedness_VV"])
                    # distances.extend(imm["OMM_dist"])
                    areas.extend(imm["area"])
                # filename = f"curvedness_{label}_long.svg"
                # twod_histogram(distances, curvatures, areas, "IMM-OMM Distance (nm)", "Curvedness (1/nm)", title="Curvedness vs Distance in the IMM", filename = filename, bins=(50,100), range=[[0,50],[0,0.1]], figsize=(4,4), log=True)
                areas = np.array(areas)
                area_list.append(areas/sum(areas))
                counts.append(count)
                curvature_list.append(curvatures)
                log_curvature_list.append(np.log(np.array(curvatures)))
            histogram(curvature_list, area_list, labels, title=f"Curvedness of {feature}", xlabel="Curvedness (1/nm)", bins=50, range=[5e-4,2e-1], filename=f"{feature}_curvedness_triple_xlog_{run}.svg",logx=True)
            histogram(curvature_list, area_list, labels, title=f"Curvedness of {feature}", xlabel="Curvedness (1/nm)", bins=50, range=[0,0.1], filename=f"{feature}_curvedness_triple_{run}.svg",logx=False)
            av0, std0 = weighted_avg_and_std(log_curvature_list[0], area_list[0])
            av1,std1 = weighted_avg_and_std(log_curvature_list[1], area_list[1])
            print(counts)
            print(stats.ttest_ind_from_stats(av0, std0, counts[0], av1, std1, counts[1]))
            print(stats.ks_2samp(curvature_list[0], curvature_list[1]))

    # Relative surface area
    labels = ["Untreated", "Treated"]
    for run in ['long', 'short', 'all']:
        area_set = []
        std_set = []
        for file, label in ((untreated, "Untreated"), (treated, "Treated"), (tggsk, "Cotreated")): 
            areas = []
            for key, val in file.tomograms.items():
                omm_area = sum(val["OMM"]["area"])

                imm = val["IMM"]
                imm_area = sum(imm['area'])/omm_area
                ibm = imm[imm["OMM_dist"]<16]
                ibm_area = sum(ibm['area'])/omm_area
                cristae = imm[imm["OMM_dist"] > 23]
                cristae_area = sum(cristae['area'])/omm_area
                junctions = imm[imm["OMM_dist"] < 23]
                junctions = junctions[junctions["OMM_dist"]>16]
                junctions_area = sum(junctions['area'])/omm_area

                areas.append((ibm_area, junctions_area, cristae_area, imm_area))
            
            areas = np.array(areas)
            mean = np.mean(areas, axis=0)
            std = np.std(areas, axis=0)
            area_set.append(mean)
            std_set.append(std)
            barchart(mean, std, labels=["IBM", "Junctions", "Cristae", "All"], title=f'Relative Surface Areas of {label}', ylabel="Surface Area Relative to OMM", filename=f"{label}_relative_areas_{run}.svg", figsize=(4,3))
        # double_barchart(area_set[0], area_set[1], std_set[0], std_set[1], labels=["IBM", "Junctions", "Cristae", "All"], title="Relative Surface Area", ylabel="Surface Area Relative to OMM", legends = labels, filename=f"double_barchart_relative_areas_{run}.svg", figsize=(4,3))
                


    # if not os.path.isdir(output_dir):
    #     os.mkdir(output_dir)
    # experiments = ["TgGSK"]
    # # # experiments = ["Treated"]
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
    # # pickle.dump(experiment_classes, "experiments.pkl")
    # # os.chdir(output_dir)
    # # untreated = experiment_classes["Untreated"]
    # # pickle.dump(untreated, "untreated.pkl")
    # untreated = experiment_classes["Untreated"]
    # with open("untreated.pkl", 'wb') as file: 
    #     pickle.dump(untreated, file)
    # tggsk = experiment_classes["TgGSK"]
    # with open("tggsk.pkl", 'wb') as file: 
    #     pickle.dump(tggsk, file)
    # # pickle.dump(tggsk, "tggsk.pkl")
    # print(tggsk)
    # print(tggsk.tomogram_names)
    # print(tggsk["GF1"])
    # print(tggsk["GF1"]["OMM"])

    # Histograms
    # Distance comparison
    

    
 




