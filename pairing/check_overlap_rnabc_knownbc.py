#!/usr/bin/env python3

import pandas as pd
import altair as alt
from typing import List, Tuple, Dict, Set, Callable
from itertools import product
import argparse
from pathlib import Path
import pysam
from splshared.var_mapping import GeneModelMapper

import numpy as np
import matplotlib.pyplot as plt

def read_starcode_histo(fpath: str) -> pd.DataFrame:
    """Read a starcode histogram file and return a pandas DataFrame."""
    with open(fpath, 'r') as f:
        lines = f.readlines()

    # Extract the header line
    header = lines[0].strip().split('\t')


def main():
    parser = argparse.ArgumentParser(description='Compute overlap of barcodes between a barcode pairing table and a set of MPSA RNA-seq barcode histograms')

    parser.add_argument('--pairing_tbl', type=str, required=True,
                      help='path to pairing table')
    
    parser.add_argument('--rnaseq_samp_tbl', type=str, required=True,
                      help='path to sample table of RNA-seq tables; requied columns = libname, bc_histo (containing path to barcode histogram in starcode format)')

    parser.add_argument('--output_plot_base', type=str, default=None,
                      help='Output plot file path')

    parser.add_argument('--lmin_barcode', type=str, required=True,
                    help='list of minimum barcode count thresholds, e.g., 1,5,10')

    args = parser.parse_args()
    
    # Read input data
    pairtbl = pd.read_table(args.pairing_tbl)
    samp_tbl = pd.read_table(args.rnaseq_samp_tbl)
    
    
if __name__ == '__main__':
    main() 