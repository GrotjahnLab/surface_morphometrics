import os
import pickle
import numpy as np
import time
# from glob import glob

from pycurv import  TriangleGraph, io
from graph_tool import load_graph

from intradistance_verticality import export_csv
from morphometrics_stats import Experiment, Tomogram, weighted_median, weighted_histogram_peak, statistics, histogram, twod_histogram, barchart, double_barchart, ks_statistics, bootstrap

conditions = ["TE", "TF", "UE", "UF"]

# basename = f"MIM029_3_lam5_ts_005_PHB.labels_IMM.AVV_rh8_refined.gt" # IMM ONLY

basename = f"*IMM.AVV_rh8_refined.gt"

junction_minimum_offset = 4 # nm 
junction_maximum_offset = 14 # nm

for condition in conditions:
    startdir = "/scratch1/users/mmedina/PHB_tomos_final_dataset/Amira_seg/"+condition+"/morphometrics/"
    # files = glob(startdir+basename)
    files = [startdir+basename]
    total_values_condition = {"IBM": [], "Junction": [], "Crista Body":[]}
    for file in files:
        ommname = file.replace("IMM", "OMM")
        omm_tg = TriangleGraph()
        omm_tg.graph=load_graph(ommname)
        imm_dist = omm_tg.graph.vp.IMM_dist.get_array()
        omm_areas = omm_tg.graph.vp.area.get_array()
        ibm_dist = weighted_histogram_peak(imm_dist, omm_areas, 100, [5,25])
        print(ibm_dist)
        junction_minimum = ibm_dist + junction_minimum_offset
        junction_maximum = ibm_dist + junction_maximum_offset
        ranges = {"IBM": [0, junction_minimum],
                  "Junction": [junction_minimum, junction_maximum], 
                  "Crista Body": [junction_maximum, 100000]}
        values = {"IBM": [], "Junction": [], "Crista Body":[]}
        
        tg = TriangleGraph()
        tg.graph = load_graph(file)
        areas = tg.graph.vp.area.get_array()
        omm_dist = tg.graph.vp.OMM_dist.get_array()
        thickness = tg.graph.vp.thickness.get_array()
        for key,value in ranges.items():
            print(key, value)
            truth = np.argwhere(np.bitwise_and(omm_dist > value[0], omm_dist<=value[1]))
            print(truth)
            print(len(truth))
            sub_thickness = thickness[truth]
            sub_areas = areas[truth]
            histogram(data=[sub_thickness], areas=[sub_areas], labels = [key], title="Thickness by Subcompartment", xlabel="Thickness (nm)", filename=file+"_"+key+".svg")
            print(key, weighted_histogram_peak(sub_thickness, sub_areas, 100, [2,7]))










        




