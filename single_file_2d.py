#! /usr/bin/env python
"""Generate an area weighted 2D histogram of a [aor pf features] feature from a single mesh from a single segmentation
Useful to get a sense of how your configuration for meshing and quantification looks. Also nice for responding to reviewer requests!

Usage: python single_file_histogram.py input_csv.csv -n1 featurename1 -n2 featurename2
"""
__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import click
import pandas as pd
from morphometrics_stats import twod_histogram

@click.command()
@click.argument('input', type=str)
@click.option('-n1', '--featurename1', type=str, required=True, default='curvedness_VV', help="Name of feature 1 (x axis) to extract")
@click.option('-n2', '--featurename2', type=str, required=True, default='curvedness_VV', help="Name of feature 2 (y axis) to extract")
def main(input, featurename1, featurename2):
    """Generate an area weighted histogram of a single feature from a single mesh from a single segmentation
    Useful to get a sense of how your configuration for meshing and quantification looks. Also nice for responding to reviewer requests!
    """
    df = pd.read_csv(input)
    data1 = df[featurename1].values
    data2 = df[featurename2].values
    areas = df['area'].values
    output_name = input[:-4] + '_' + featurename1 +"_"+featurename2+'.svg'
    print('Generating 2D histogram for feature {} vs feature {}'.format(featurename2,featurename1))
    twod_histogram(data1, data2,areas, featurename1, featurename2, os.path.basename(input)+ " "+ featurename1+" vs "+featurename2, filename=output_name)

if __name__ == '__main__':
    main()