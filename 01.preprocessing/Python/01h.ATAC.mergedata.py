#!/bin/python

# example: python 01h.ATAC.mergedata.py -i samle_info_celltypeage.txt -o /peak_calling/subset/


import argparse

parser = argparse.ArgumentParser(description="extract cell data from each sample")
parser.add_argument("-i", "--input", type=str, dest="input", required=True, help="input sample info. file")
#parser.add_argument("-c", "--cell", type=str, dest="cell", required=True, help="input cell info. file")
#parser.add_argument("-p", "--processes", type=int, dest="processes", help="number of processes to use")
parser.add_argument("-o", "--output", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
print(snap.__version__)
import numpy as np
#from multiprocessing import Pool
import sys


#def parse_cluster_cells(cells):
#	dict_data_cells = {}
#	for row in cells:
#		sample = (row[0], row[2])
#		cell = row[1]
#		if sample in dict_data_cells:
#			dict_data_cells[sample].append(cell)
#		else:
#			dict_data_cells[sample] = [cell]
#	return dict_data_cells

def main():
	sample_file = args.input
#	cell_file = args.cell
#	processes = args.processes
	out_prefix = args.output
	print("Step1: read data!")
	sample_info = np.loadtxt(sample_file, dtype=str)
	dataset = {}
	for i in range(len(sample_info)):
		hdfile = sample_info[i][2]
		index = (sample_info[i][0], sample_info[i][1])
		if index in dataset:
			dataset[index].append(hdfile)
		else:
			dataset[index] = [hdfile]

	print("Step2: merge data!")

	dataset_subset = {}
	for sample, hdfiles in dataset.items():
		if len(hdfiles) != 2:
			print("Error: num of hdfile is not 2!!!")
			break
		hdfiles_bk = [(hdfile, snap.read(hdfile, backed='r')) for hdfile in hdfiles] 
		output_path = "".join([out_prefix, ".".join(sample), ".h5ads"])
		merge_dataset = snap.AnnDataSet(
			adatas=hdfiles_bk,
			filename=output_path,
			add_key='sample',
		)
		merge_dataset.obs_names = merge_dataset.obs['sample'].to_numpy() + ":" + merge_dataset.obs_names

	print("Subset dataset written to file:", output_path)

	sys.exit(0)

if __name__ == "__main__":
	main()

