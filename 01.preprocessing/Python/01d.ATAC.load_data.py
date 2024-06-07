#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="load fragment file to snap")
parser.add_argument("-i", "--input", type=str, dest="input", help="input fragment file")
parser.add_argument("-c", "--cells", type=str, default=None, dest="cells", help="cell barcode file")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ #'2.3.2dev'
from pathlib import Path
import numpy as np
import pandas as pd

# Input files
frag_path = args.input
cells_file = args.cells
out_prefix = args.output
file = Path(frag_path)

# load data
h5ad_filename = out_prefix + ".h5ad"
if cells_file is None:
	data = snap.pp.import_data(file, chrom_sizes=snap.genome.mm10, file=h5ad_filename, sorted_by_barcode=False, min_num_fragments=0)
else:
	data = snap.pp.import_data(file, chrom_sizes=snap.genome.mm10, file=h5ad_filename, sorted_by_barcode=False, min_num_fragments=0, whitelist=cells_file)


ATAC_cells_raw = np.array(data.obs_names)
cells_filename = out_prefix + ".cells_raw.txt"
np.savetxt(cells_filename, ATAC_cells_raw, fmt="%s", delimiter=',')
data.close()

