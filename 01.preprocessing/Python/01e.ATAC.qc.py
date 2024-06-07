#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="preprocess raw snap file")
parser.add_argument("-i", "--input", type=str, dest="input", help="input raw snap file")
parser.add_argument("-p", "--plot", action="store_true", dest="plot", default=False, help="output qc plot")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ #'2.3.2dev'
import numpy as np
import pandas as pd

# Input files
raw_path = args.input
plot = args.plot
out_prefix = args.output

# preprocess
data = snap.read(raw_path, backed='r')
out_file = out_prefix + ".qc.h5ad"
data.copy(out_file)
data.close()
data = snap.read(out_file)

if plot:
	out_plot = out_prefix + ".qc.pdf"
	snap.pl.tsse(data, out_file=out_plot, show=False)

snap.pp.filter_cells(data, min_counts=1000, min_tsse=10)

ATAC_cells_qc = np.array(data.obs_names)
cells_filename = out_prefix + ".cells_qc.txt"
np.savetxt(cells_filename, ATAC_cells_qc, fmt="%s", delimiter=',')
data.close()





