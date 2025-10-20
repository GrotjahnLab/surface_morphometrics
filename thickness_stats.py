"""Script to generate all the plots and statistics for the joint 
analysis in the bilayer thickness paper."""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import pickle
import numpy as np
import time
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


from morphometrics_stats import Experiment, Tomogram, weighted_median, weighted_histogram_peak, statistics, histogram, twod_histogram, barchart, double_barchart, ks_statistics, bootstrap


startdir = os.getcwd()+"/pkl/2_pkl/"
output_dir = os.getcwd()+"/stats_analysis/"

print("Unpickling UE")
with open(startdir+'UE.pkl', 'rb') as file:
    ue = pickle.load(file)
print("Unpickling UF")
with open(startdir+'UF.pkl', 'rb') as file:
    uf = pickle.load(file)
print("Unpickling Other")
with open(startdir+'other_org_4.pkl', 'rb') as file:
    other = pickle.load(file)

# Remove cases where the thickness failed to move from the exact initialization value - these are artifactual
for groups in ue,uf,other:
    for key in groups.tomograms:
        tomo = groups.tomograms[key]
        for name in tomo.dataframe_names:
            tomo[name] = tomo[name][tomo[name].thickness != 4.00]
    
# Test
for key, value in enumerate(ue.tomograms):
    tomo = ue.tomograms[value]
    print(weighted_histogram_peak(tomo["IMM"].thickness, tomo["IMM"].area, 50, [2,6]), weighted_median(tomo["IMM"].thickness, tomo["IMM"].area))

# Calculate thickness statistics
# Violin plot
thicknesses = {}
total_thicknesses = {}
total_areas = {}
for toms in ue,uf,other:
    for key in toms.tomograms:
        tomo = toms.tomograms[key]
        for name in tomo.dataframe_names:
            name = str(name)
            if name == "COP":
                continue
            if name == "ER":
                continue
            if name == "Golgi":
                continue
            if name not in thicknesses:
                thicknesses[name] = []
                total_thicknesses[name] = []
                total_areas[name] = []
            if "unique_component" not in tomo[name].columns:
                thicknesses[name].append(weighted_median(tomo[name].thickness, tomo[name].area))
                continue

            for component in tomo[name].unique_component.unique():
                df = tomo[name][tomo[name].unique_component == component]
                if sum(df.area) < 1000:
                    continue
                thicknesses[name].append(weighted_median(df.thickness, df.area))
            # thicknesses[name].append(weighted_histogram_peak(tomo[name].thickness, tomo[name].area, 100, [1.5,7.5]))
            total_thicknesses[name].extend(df.thickness)
            total_areas[name].extend(df.area)


print(thicknesses)
print(thicknesses.keys())
conditions = thicknesses.keys()
# conditions = list(conditions)
conditions=['IMM', 'OMM', 'Smooth_ER', 'Ribo_ER', 'Vesicles'] # , 'Golgi', 'COP'
morphologies = ["" for i in conditions]
thicknesses_all = [total_thicknesses[name] for name in conditions]
print(len(thicknesses_all))
areas_all = [total_areas[name] for name in conditions]
histogram(thicknesses_all, areas_all, conditions, "Thickness (nm)", filename="thick_hist.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)

# conditions=['IMM', 'OMM', 'Smooth_ER', 'Ribo_ER', 'Vesicles'] # 
morphologies = ["" for i in conditions]
thickness_sets = [thicknesses[i] for i in conditions]
statistics(thickness_sets, "Thickness", conditions, morphologies, test_type="medians", filename="thickness_violin_main.csv", ylabel="Distance (nm)")

# # Curvature statistics
# # Calculate thickness statistics
# # Violin plot
# # thicknesses = {}
# # total_thicknesses = {}
# # total_areas = {}
# # for list in ue,uf,other:
# #     for key in list.tomograms:
# #         tomo = list.tomograms[key]
# #         for name in tomo.dataframe_names:
# #             name = str(name)
# #             if name not in thicknesses:
# #                 thicknesses[name] = []
# #                 total_thicknesses[name] = []
# #                 total_areas[name] = []
# #             thicknesses[name].append(weighted_histogram_peak(tomo[name].thickness, tomo[name].area), 100, [0.001,0.02])
# #             # thicknesses[name].append(weighted_histogram_peak(tomo[name].thickness, tomo[name].area, 100, [1.5,7.5]))
# #             total_thicknesses[name].extend(tomo[name].thickness)
# #             total_areas[name].extend(tomo[name].area)


# # print(thicknesses)
# # print(thicknesses.keys())
# # conditions=['IMM', 'OMM', 'ER', 'Smooth_ER', 'Ribo_ER', 'Vesicles'] # , 'Golgi', 'COP'
# # morphologies = ["" for i in conditions]
# # thickness_sets = [thicknesses[i] for i in conditions]
# # thicknesses_all = [total_thicknesses[name] for name in conditions]
# areas_all = [total_areas[name] for name in conditions]

# # histogram(thicknesses_all, areas_all, conditions, "Thickness (nm)", filename="thick_hist.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)
# # statistics(thickness_sets, "Thickness", conditions, morphologies, test_type="medians", filename="thickness_violin.csv", ylabel="Distance (nm)")


# # # Make a matplotlib violin plor with conditions as x-axis
# # plt.figure(figsize=(10, 6))
# # plt.violinplot(thickness_sets, showmeans=True, showmedians=True)
# # plt.xticks(range(1, len(conditions) + 1), conditions)
# # plt.xlabel("Conditions")
# # plt.ylabel("Thickness (nm)")
# # plt.title("Violin plot of thicknesses")
# # plt.savefig(output_dir+"violin_plot.png")
# # plt.savefig(output_dir+"violin_plot.svg")


# thicknesses = []
# total_thicknesses = []

# conditions= ["Extended"]*3+["Fragmented"]*3
# compartments = ["IBM", "Junction", "Crista Body"]*2
# names = [conditions[i]+" "+compartments[i] for i in range(len(conditions))]
# j = 0
# thicknesses = []
# total_thicknesses = []
# total_areas = []
# for list in ue,uf:
#     for i, sub in enumerate(["IBM", "Junction", "Crista Body"]):
#         total_thicknesses.append([])
#         total_areas.append([])
#         thicknesses.append([])
#         for key in list.tomograms:
#             tomo = list.tomograms[key]
#             df = tomo["IMM"]
#             if sub == "IBM":
#                 df = df[df.subcompartment == 0]
#             elif sub == "Junction":
#                 df = df[df.subcompartment == 1]
#             elif sub == "Crista Body":
#                 df = df[df.subcompartment == 2]
#             else:
#                 raise ValueError("Invalid compartment name")
#             thicknesses[j].append(weighted_median(df.thickness, df.area))
#             total_thicknesses[j].extend(df.thickness)
#             total_areas[j].extend(df.area)
#         j += 1

# histogram(total_thicknesses, total_areas, names, "Thickness (nm)", filename="thick_hist_subcompartments.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)
# statistics(thicknesses, "Thickness", conditions, compartments, test_type="medians", filename="thickness_violin_subcompartments.csv", ylabel="Distance (nm)")


# Combined UE and UF
compartments = ["IBM", "Junction", "Crista Body"]
thicknesses = []
total_thicknesses = []
total_areas = []
for i, sub in enumerate(["IBM", "Junction", "Crista Body"]):
    total_thicknesses.append([])
    total_areas.append([])
    thicknesses.append([])
    for list in ue,uf:
        for key in list.tomograms:
            tomo = list.tomograms[key]
            df = tomo["IMM"]
            if sub == "IBM":
                df = df[df.subcompartment == 0]
            elif sub == "Junction":
                df = df[df.subcompartment == 1]
            elif sub == "Crista Body":
                df = df[df.subcompartment == 2]
            else:
                raise ValueError("Invalid compartment name")
            thicknesses[i].append(weighted_median(df.thickness, df.area))
            total_thicknesses[i].extend(df.thickness)
            total_areas[i].extend(df.area)

histogram(total_thicknesses, total_areas, compartments, "Thickness (nm)", filename="thick_hist_subcompartments_all.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)
statistics(thicknesses, "Thickness", [""]*3, compartments, test_type="medians", filename="thickness_violin_subcompartments_all.csv", ylabel="Distance (nm)")

# Thickness for Fragmented and Extended

conditions = ["Fragmented"]*2+["Extended"]*2
compartments = ["IMM", "OMM"]*2

thicknesses = []
total_thicknesses = []
total_areas = []
i=0
for list in ue,uf:
    for compartment in ["IMM", "OMM"]:
        total_thicknesses.append([])
        total_areas.append([])
        thicknesses.append([])
        for key in list.tomograms:
            tomo = list.tomograms[key]
            # Not always ER
            if not compartment in tomo.dataframe_names:
                continue
            df = tomo[compartment]
            thicknesses[i].append(weighted_median(df.thickness, df.area))
            total_thicknesses[i].extend(df.thickness)
            total_areas[i].extend(df.area)
        i += 1

names = [conditions[j]+" "+compartments[j] for j in range(len(conditions))]
histogram(total_thicknesses, total_areas, names, "Thickness (nm)", filename="thick_hist_treatment.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)
statistics(thicknesses, "Thickness", conditions, compartments, test_type="medians", filename="thick_violin_treatment.csv", ylabel="Distance (nm)")

# # Curvature vs thickness
# # generate 2d histogram of curvature vs thickness
# # for IMM
# curvature = []
# thickness = []
# curvature_areas = []
# for list in ue,uf:
#     for key in list.tomograms:
#         tomo = list.tomograms[key]
#         if "IMM" in tomo.dataframe_names:
#             df = tomo["IMM"]
#             df=df[df.thickness > 0]
#             curvature.extend(df.curvedness_VV)
#             thickness.extend(df.thickness)
#             curvature_areas.extend(df.area)


# curvature = np.log10(curvature)
# twod_histogram(curvature, thickness, curvature_areas, "Log Curvedness (Log 1/nm)", "Thickness (nm)", "Thickness vs Curvedness", filename="curvature_thickness.svg", bins=[100,100], range=[[-3, -1], [2.5, 4.5]], log=True)
# reg = LinearRegression().fit(np.reshape(curvature, (-1,1)), np.reshape(thickness,(-1,1)))
# print(reg.score(np.reshape(curvature, (-1,1)), np.reshape(thickness,(-1,1))), reg.intercept_, reg.coef_)

# verticality vs thickness
# generate 2d histogram of curvature vs thickness
# for IMM
verticality = []
thickness = []
curvature_areas = []
for list in ue,uf:
    for key in list.tomograms:
        tomo = list.tomograms[key]
        if "IMM" in tomo.dataframe_names:
            df = tomo["IMM"]
            df=df[df.thickness > 0]
            verticality.extend(df.verticality)
            thickness.extend(df.thickness)
            curvature_areas.extend(df.area)


twod_histogram(verticality, thickness, curvature_areas, "Verticality (º)", "Thickness (nm)", "Thickness vs Verticality", filename="verticality_thickness.svg", bins=[100,100], range=[[0, 90], [2.5, 4.5]], log=True)
reg = LinearRegression().fit(np.reshape(verticality, (-1,1)), np.reshape(thickness,(-1,1)))
print(reg.score(np.reshape(verticality, (-1,1)), np.reshape(thickness,(-1,1))), reg.intercept_, reg.coef_)

# # Curvature vs thickness cristae
# # generate 2d histogram of curvature vs thickness
# # for IMM
# curvature = []
# thickness = []
# curvature_areas = []
# for list in ue,uf:
#     for key in list.tomograms:
#         tomo = list.tomograms[key]
#         if "IMM" in tomo.dataframe_names:
#             df = tomo["IMM"]
#             df = df[df.subcompartment == 2]
#             df=df[df.thickness > 0]
#             curvature.extend(df.curvedness_VV)
#             thickness.extend(df.thickness)
#             curvature_areas.extend(df.area)

# curvature = np.log10(curvature)
# twod_histogram(curvature, thickness, curvature_areas, "Log Curvedness (Log 1/nm)", "Thickness (nm)", "Thickness vs Curvedness (Cristae)", filename="curvature_thickness_cristae.svg", bins=[100,100], range=[[-3, -1], [2.5, 4.5]], log=True)
# reg = LinearRegression().fit(np.reshape(curvature, (-1,1)), np.reshape(thickness,(-1,1)))
# print(reg.score(np.reshape(curvature, (-1,1)), np.reshape(thickness,(-1,1))), reg.intercept_, reg.coef_)

# # Calculate average thickness for IMM at different curvature quantiles
# quantile = [0,0.5,0.9,0.95,0.99,1]
# names = ["0-0.5", "0.5-0.9", "0.9-0.95", "0.95-0.99", "0.99-1"]
# conditions = [""]*5
# thickness = [[],[],[],[],[]]
# thicknesses_all = [[],[],[],[],[]]
# areas_all = [[],[],[],[],[]]
# for list in ue,uf:
#     for key in list.tomograms:
#         tomo = list.tomograms[key]
#         if "IMM" in tomo.dataframe_names:
#             df = tomo["IMM"]
#             df=df[df.thickness > 0]
#             for i in range(len(quantile)-1):
#                 df_q = df[(df.curvedness_VV > np.quantile(df.curvedness_VV, quantile[i])) & (df.curvedness_VV < np.quantile(df.curvedness_VV, quantile[i+1]))]
#                 thickness[i].append(weighted_median(df_q.thickness, df_q.area))
#                 thicknesses_all[i].extend(df_q.thickness)
#                 areas_all[i].extend(df_q.area)

# statistics(thickness, "Thickness", names, conditions, test_type="medians", filename="thickness_quantiles.csv", ylabel="Thickness (nm)")
# histogram(thicknesses_all, areas_all, names, "Thickness (nm)", filename="thick_hist_quantiles.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)     

# # Calculate average thickness for IMM at different curvature quantiles
# quantile = [0,0.25,0.5,0.75,1]
# names = ["0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"]
# conditions = [""]*4
# thickness = [[],[],[],[]]

# for list in ue,uf:
#     for key in list.tomograms:
#         tomo = list.tomograms[key]
#         if "IMM" in tomo.dataframe_names:
#             df = tomo["IMM"]
#             df=df[df.thickness > 0]
#             for i in range(len(quantile)-1):
#                 df_q = df[(df.curvedness_VV > np.quantile(df.curvedness_VV, quantile[i])) & (df.curvedness_VV < np.quantile(df.curvedness_VV, quantile[i+1]))]
#                 thickness[i].append(weighted_median(df_q.thickness, df_q.area))

# statistics(thickness, "Thickness", names, conditions, test_type="medians", filename="thickness_quantiles_2.csv", ylabel="Thickness (nm)")

# # Calculate average thickness for IMM at different curvature quantiles
# quantile = [0,0.33,0.66,1]
# names = ["0-0.33", "0.33-0.66", "0.66-1"]
# conditions = [""]*3
# thickness = [[],[],[]]
# for list in ue,uf:
#     for key in list.tomograms:
#         tomo = list.tomograms[key]
#         if "IMM" in tomo.dataframe_names:
#             df = tomo["IMM"]
#             df=df[df.thickness > 0]
#             for i in range(len(quantile)-1):
#                 df_q = df[(df.curvedness_VV > np.quantile(df.curvedness_VV, quantile[i])) & (df.curvedness_VV < np.quantile(df.curvedness_VV, quantile[i+1]))]
#                 thickness[i].append(weighted_median(df_q.thickness, df_q.area))

# statistics(thickness, "Thickness", names, conditions, test_type="medians", filename="thickness_quantiles_3.csv", ylabel="Thickness (nm)")
                     

# Measure OMM thickness when <50nm from ER vs when more than 80nm from ER

# thicknesses = [[],[]]
# total_thicknesses = [[],[]]
# total_areas = [[],[]]
# conditions = ["ER Contact", "Non-ER Contact"]
# morphologies = [""]*2
# for list in [ue,uf]:
#     for key in list.tomograms:
#         tomo = list.tomograms[key]
#         if "ER" in tomo.dataframe_names:
#             df = tomo["OMM"]
#             print(df.columns)
#             for i in range(2):
#                 if i == 0:
#                     df_q = df[(df.ER_dist < 40)]
#                     if len(df_q) == 0:
#                         continue
#                 else:
#                     df_q = df[(df.ER_dist >= 40)]

#                 thicknesses[i].append(weighted_median(df_q.thickness, df_q.area))
#                 total_thicknesses[i].extend(df_q.thickness)
#                 total_areas[i].extend(df_q.area)

# histogram(total_thicknesses, total_areas, conditions, "Thickness (nm)", filename="thick_hist_ER_contact.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)
# statistics(thicknesses, "Thickness", conditions, morphologies, test_type="medians", filename="thickness_ER_contact.csv", ylabel="Distance (nm)")


# Calculate thickness statistics
# Violin plot
thicknesses = {}
total_thicknesses = {}
total_areas = {}
for list in [other]:
    for key in list.tomograms:
        tomo = list.tomograms[key]
        for name in tomo.dataframe_names:
            name = str(name)

            if name not in thicknesses:
                thicknesses[name] = []
                total_thicknesses[name] = []
                total_areas[name] = []

            thicknesses[name].append(weighted_median(tomo[name].thickness, tomo[name].area))
                # thicknesses[name].append(weighted_histogram_peak(tomo[name].thickness, tomo[name].area, 100, [1.5,7.5]))
            total_thicknesses[name].extend(tomo[name].thickness)
            total_areas[name].extend(tomo[name].area)


print(thicknesses)
print(thicknesses.keys())
conditions=['Smooth_ER', 'Ribo_ER'] # , 'Golgi', 'COP'
morphologies = ["" for i in conditions]
thickness_sets = [thicknesses[i] for i in conditions]
thicknesses_all = [total_thicknesses[name] for name in conditions]
areas_all = [total_areas[name] for name in conditions]

histogram(thicknesses_all, areas_all, conditions, "Thickness (nm)", filename="thick_hist_smoothrough.svg", bins=100, range=[1.5, 7.5], xlabel="Thickness (nm)", vlines=False)
statistics(thickness_sets, "Thickness", conditions, morphologies, test_type="medians", filename="thickness_violin_ERsub.csv", ylabel="Distance (nm)")
