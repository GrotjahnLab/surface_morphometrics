#! /usr/bin/env python
"""Generate an area weighted histogram of a single feature from a single mesh from a single segmentation
Useful to get a sense of how your configuration for meshing and quantification looks. Also nice for responding to reviewer requests!

Usage: python single_file_histogram.py input_csv.csv -n featurename
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"


import os
import click
import pandas as pd
from morphometrics_stats import histogram, weighted_median

@click.command()
@click.argument('input', type=str)
@click.option('-n', '--featurename', type=str, required=True, default='curvedness_VV', help="Name of feature to extract")
@click.option('-s','--figuresize', type=tuple, required=False, default=(11,8.5),help="Size of figure in inches (x,y)")
def main(input, featurename, figuresize):
    """Generate an area weighted histogram of a single feature from a single mesh from a single segmentation
    Useful to get a sense of how your configuration for meshing and quantification looks. Also nice for responding to reviewer requests!
    """
    df = pd.read_csv(input)
    data = df[featurename].values
    areas = df['area'].values
    output_name = input[:-4] + '_' + featurename + '.svg'
    print('Generating histogram for feature {}'.format(featurename))
    print('Area Weighted Median: {}'.format(weighted_median(data, areas)))
    histogram([data], [areas], [os.path.basename(input)[:-4]], os.path.basename(input)+ " "+ featurename, featurename, filename=output_name, figsize=figuresize)

if __name__ == '__main__':
    main()
