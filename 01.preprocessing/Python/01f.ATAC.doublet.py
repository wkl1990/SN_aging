#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="remove doublet from snap file")
parser.add_argument("-i", "--input", type=str, dest="input", help="input qc snap file")
parser.add_argument("-b", "--blacklist", type=str, default=None, dest="blacklist", help="blacklist file")
#parser.add_argument("-p", "--plot", action="store_true", dest="plot", default=False, help="output qc plot")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ #'2.3.2dev'
import numpy as np
import pandas as pd

# Input files
qc_path = args.input
blacklist = args.blacklist
#plot = args.plot
out_prefix = args.output

# preprocess
data = snap.read(qc_path, backed='r')
out_file = out_prefix + ".doublet.h5ad"
data.copy(out_file)
data.close()
data = snap.read(out_file)

#if plot:
#	out_plot = out_prefix + ".doublet.pdf"
#	snap.pl.scrublet(data, out_file=out_plot, show=False)

snap.pp.add_tile_matrix(data)
snap.pp.select_features(data, n_features=500000, blacklist=blacklist)
snap.pp.scrublet(data)
snap.pp.filter_doublets(data)

ATAC_cells_doublet = np.array(data.obs_names)
cells_filename = out_prefix + ".cells_doublet.txt"
np.savetxt(cells_filename, ATAC_cells_doublet, fmt="%s", delimiter=',')
data.close()





