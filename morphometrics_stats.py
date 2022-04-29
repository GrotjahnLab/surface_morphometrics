#! /usr/bin/env python
"""Library of useful functions for statistical analysis of morphometrics data

If called as a script, this will assemble an Experiment object based on the config.yml file
Usage: python morphometrics_stats.py config.yml experiment_name"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import copy
import glob
import os
import pickle

import click
import yaml
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
    """Experiments are containers for a series of tomograms and their pandas dataframes
    
    Access pattern: experiment[tomoname][label]."""
    def __init__(self, name):
        self.name = name
        self.tomograms = {}
        self.tomogram_names = set()
    

    def add_tomograms(self, names, surface_bases, folder, file_extension):
        """Add tomograms to the experiment. Names must be the basenames used by pycurv.
        
        Args:
            names (list): list of tomogram prefixes
            surface_bases (list): list of label names
            folder (str): path to the folder containing surfaces
            file_extension (str): file extension of the surface files"""
        for name in names:
            # print(name)
            surfs = []
            fns = []
            for surface_base in surface_bases:
                filename = folder+name+"_"+surface_base+file_extension
                if os.path.isfile(filename):
                    surfs.append(surface_base)
                    fns.append(filename)
            self.tomograms[name] = Tomogram(name, surface_bases, fns)

    def __getitem__(self, name):
        return self.tomograms[name]
    
    def __setitem__(self, name, tomo):
        self.tomogram_names.add(name)
        self.tomograms[name] = tomo

class Tomogram():
    """A tomogram class, which can have a number of different dataframes associated with it.
    Each dataframe is associated with a differen segmented feature, such as ER or IMM.
"""
    
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
    """compute the weighted median of values list. The 
    weighted median is computed as follows:
    1- sort both lists (values and weights) based on values.
    2- select the 0.5 point from the weights and return the corresponding values as results
    e.g. values = [1, 3, 0] and weights=[0.1, 0.3, 0.6] assuming weights are probabilities.
    sorted values = [0, 1, 3] and corresponding sorted weights = [0.6, 0.1, 0.3] the 0.5 point on
    weight corresponds to the first item which is 0. so the weighted median is 0.
    
    Function provided by Max Ghenis on stack overflow, CC-BY-SA
    
    Args:
        values (array-like): list of values
        weights (array-like): list of weights"""

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
    '''compute the peak of a histogram for a list of values.
    
    Args:
        values (array-like): list of values
        weights (array-like): list of weights
        bins (array-like): list of bins
        bin_range (tuple): (min, max) of the histogram
    '''
    hist, bin_edges = np.histogram(values, bins=bins, range=bin_range, weights=weights, density=True)
    i = np.argmax(hist)
    bin_vals = (bin_edges[i], bin_edges[i+1])
    return(np.mean(bin_vals))


def statistics(datasets, basename,  condition_names, morph_names, test_type="median", filename="test.csv", figsize=(5,4), ylabel="Peak Value"):
    """Generate mann-whitney U test and violin plots for a set of datasets.
    Args:
        datasets (list): list of datasets to compare
        basename (str): base name of the output files
        condition_names (list): list of condition names
        morph_names (list): list of morph names
        test_type (str): type of test to perform.
        filename (str): name of the output file
        figsize (tuple): size of the figure
        ylabel (str): label of the y axis
        """
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


def histogram(data, areas, labels, title, xlabel, filename="hist.svg", bins=50, range=None, figsize=(5,3.6), logx=False, vlines = True):
    """Construct an area-weighted histogram of the data.
    Args:
        data (array-like): list of arrays to be independently plotted.
        areas (array-like): list of area values for each data array.
        labels (array-like): list of labels for each data array.
        title (str): title of the plot.
        xlabel (str): x-axis label.
        filename (str): name of the output file.
        bins (int): number of bins.
        range (list): range of the x-axis.
        figsize (tuple): figure size.
        logx (bool): whether to use a logarithmic x-axis.
        vlines (bool): whether to plot vertical lines at the histogram peaks.
        """
    assert len(data)==len(areas)
    assert len(data)==len(labels)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    if logx:
        bins = np.logspace(np.log10(range[0]),np.log10(range[1]), bins)
        ax.set_xscale("log")
    for index, value in enumerate(data):
        n,binset,_ = ax.hist(value, bins=bins, weights=areas[index], label=labels[index], ec=colors[index], fc=colors_light[index],histtype="stepfilled", density=False, range = range)
        if vlines: 
            delta = (binset[1]-binset[0])/2
            idx = n.argmax()
            print(idx)
            ax.axvline(binset[idx]+delta, linestyle="--", color=colors[index])
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

def twod_histogram(data1, data2, areas, data1_label, data2_label, title, bins=(50,50), range=None, filename="twod_hist.svg", figsize=(5,5), log=True):
    """Construct an area-weighted 2D histogram of the data.
    
    Args:
        data1 (list): list of values to for x axis.
        data2 (list): list of values to for y axis.
        areas (list): list of area values for weighting.
        data1_label (str): label for x axis.
        data2_label (str): label for y axis.
        title (str): title of the plot.
        bins (tuple): number of bins for x and y axis.
        range (list): range of the x and y axis. (format is [xmin, xmax, ymin, ymax])
        filename (str): name of the output file.
        figsize (tuple): figure size.
        log (bool): whether to use a logarithmic x-axis.
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title)
    ax.hist2d(data1, data2, bins=bins, range=range, density=True, weights=areas, norm=matplotlib.colors.LogNorm(), cmap=my_cmap)
    ax.set_xlabel(data1_label)
    ax.set_ylabel(data2_label)
    if range:
        ax.set_xlim(range[0])
        ax.set_ylim(range[1])
    plt.tight_layout()
    fig.savefig(filename, bbox_inches='tight')
    fig.savefig(filename[:-3]+"png", bbox_inches='tight')
    return True

def barchart(bars, errors, labels, title, ylabel, filename="barchart.svg", figsize=(4,3), ymax=None, hline=None):
    """Construct a bar chart with error bars.
    
    Args:
        bars (list): list of bar heights.
        errors (list): list of error bars.
        labels (list): list of bar labels.
        title (str): title of the plot.
        ylabel (str): y-axis label.
        filename (str): name of the output file.
        figsize (tuple): figure size.
        ymax (float): maximum value of the y-axis.
        hline (float): Optional value of the horizontal line.
        """
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
    """Construct a side by side bar chart for 2 sets of data.

        Args:
            bars1 (list): list of values for the first set of bars.
            bars2 (list): list of values for the second set of bars.
            errors1 (list): list of errors for the first set of bars.
            errors2 (list): list of errors for the second set of bars.
            labels (list): list of labels for the bars.
            title (str): title of the plot.
            ylabel (str): label for the y-axis.
            legends (list): list of legends for the two datasets
            filename (str): name of the output file.
            figsize (tuple): figure size."""
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

@click.command()
@click.argument('configfile', type=click.Path(exists=True), required=True)
@click.argument('experimentname', required=True)
def assemble_experiment_pickle(configfile, experimentname):
    """Assemble a pickle file for an experiment from all the tomos in the config folder."""
    with open(configfile) as f:
        config = yaml.safe_load(f)
        if not config["data_dir"]:
            print("data_dir not specified in config.yml")
            exit()
        elif not config["data_dir"].endswith("/"):
            config["data_dir"] += "/"
        if not config["work_dir"]:
            print("work_dir not specified in config.yml - data_dir will be used for output")
            config["work_dir"] = config["data_dir"]
        elif not config["work_dir"].endswith("/"):
            config["work_dir"] += "/"
        
        output_file = config["work_dir"] + experimentname+".pkl"

        input_files = glob.glob(config["data_dir"] + "*.mrc")
        input_names = [os.path.basename(x)[:-4] for x in input_files]
        rh = config["curvature_measurements"]["radius_hit"]
        extension = ".AVV_rh{}.csv".format(rh)
        labels = list(config["segmentation_values"].keys())
        experiment = Experiment(experimentname)

        experiment.add_tomograms(input_names, labels, config["work_dir"], extension) 
        with open(output_file, "wb") as f:
            pickle.dump(experiment, f)
        print(output_file)






if __name__ == "__main__":
    assemble_experiment_pickle()