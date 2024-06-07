#!/bin/python

# example: python 03f.ATAC.pmat.py -i final_cells.h5ad -b /peak_calling/macs2/SN.macs2.union.peaks.bed -o /pmat/final_cells.pmat


import argparse

parser = argparse.ArgumentParser(description="extract cell data from h5ad sample")
parser.add_argument("-i", "--input", type=str, dest="input", required=True, help="input sample h5ad file")
parser.add_argument("-b", "--bed", type=str, dest="bed", required=True, help="input bed file")
#parser.add_argument("-p", "--processes", type=int, dest="processes", help="number of processes to use")
parser.add_argument("-o", "--output", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
print(snap.__version__)
import numpy as np
#from multiprocessing import Pool
import sys


def main():
	sample_file = args.input
	bed_file = args.bed
#	processes = args.processes
	out_prefix = args.output
	print("Step1: read data!")
	hdfile = snap.read(sample_file, backed='r')

	print("Step2: get pmat data!")
	out_file = ".".join([out_prefix, "h5ad"])
	snap.pp.make_peak_matrix(adata = hdfile,
	                        inplace = False,
	                        file = out_file,
	                        backend = "hdf5",
	                        peak_file = bed_file,
	                        chunk_size = 10000,
	                        use_x = False)


	sys.exit(0)

if __name__ == "__main__":
	main()






