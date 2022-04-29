#! /usr/bin/env python
"""Script to generate all the plots and statistics for the mitochondria ultrastructure 
analysis in the surface morphometrics paper."""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import pickle
import numpy as np

from morphometrics_stats import Experiment, Tomogram, weighted_median, weighted_histogram_peak, statistics, histogram, twod_histogram, barchart, double_barchart


junction_minimum = 19 # nm 
junction_maximum = 40 # nm

startdir = os.getcwd()
output_dir = os.getcwd()+"/output/"
print("Unpickling Treated")
with open('treated.pkl', 'rb') as file:
    treated = pickle.load(file)
print("Unpickling Untreated")
with open('untreated.pkl', 'rb') as file:
    untreated = pickle.load(file)
print("Unpickling Cotreatment")
with open('tggsk.pkl', 'rb') as file:
    tggsk = pickle.load(file)


# # Prepare to write out graphs!
os.chdir(output_dir)

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
