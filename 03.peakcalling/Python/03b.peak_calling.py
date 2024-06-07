#!/bin/python

# example: python 03b.peak_calling.py -i c1.h5ads -p c1 -o /snapatac2/peaks/

import argparse

parser = argparse.ArgumentParser(description="call peaks for pool and replications")
parser.add_argument("-i", "--input", type=str, dest="input", required=True, help="input h5ad file")
parser.add_argument("-r", "--rep", type=str, dest="replication", help="name of the replication info. in h5ad")
parser.add_argument("-m", "--macs", action="store_true", dest="macs", default=False, help="call peaks use masc2 instead of snapatac2")
parser.add_argument("-s", "--split", action="store_true", dest="split", default=False, help="split the bed file instead of shuffling cells")
parser.add_argument("-b", "--bigwig", action="store_true", dest="bigwig", default=False, help="export tn5 bigwig file by snapatac2")
parser.add_argument("-p", "--prefix", type=str, dest="prefix", help="prefix of output file")
#parser.add_argument("-p", "--processes", type=int, dest="processes", help="number of processes to use")
parser.add_argument("-o", "--output", type=str, dest="output", help="output path")

args = parser.parse_args()

import snapatac2 as snap
print(snap.__version__)
import numpy as np
#from multiprocessing import Pool
import sys
import os

def main():
	h5ad_file = args.input
	biorep = args.replication
	out_path = args.output
	prefix = args.prefix
	macs = args.macs
	split = args.split
	bigwig = args.bigwig
	if not out_path.endswith('/'):
		out_path = out_path + '/'
	print("Step1: read data!")
	subset_dataset = snap.read_dataset(h5ad_file)
	print("Step2: check replications!")
	subset_dataset.obs["pool"] = np.repeat("pool", subset_dataset.n_obs)
	if biorep is None:
		ids = subset_dataset.obs["sample"].to_numpy().astype('U')
		sam_rep = np.char.rsplit(ids, sep="_", maxsplit=1)
		rep = np.array([r[1] for r in sam_rep])
		subset_dataset.obs["biorep"] = rep
		biorep = "biorep"
	else:
		rep = subset_dataset.obs[biorep].to_numpy()
	if not split:
		np.random.seed(seed=2023)
		np.random.shuffle(rep)
		replace_dict = {"rep1": "pseudo1", "rep2": "pseudo2"}
		pseudo = np.copy(rep).astype('U7')
		for key, value in replace_dict.items():
			mask = rep == key
			pseudo[mask] = value
		subset_dataset.obs["pseudo"] = pseudo
	print("Step3: call peaks!")
	if macs:
		snap.ex.export_fragments(subset_dataset, groupby="pool", selections=None, ids=None, out_dir=out_path, prefix=prefix, suffix='.bed.gz')
		os.system("mv {0}{1}pool.bed.gz {0}{1}.bed.gz".format(out_path, prefix))
		snap.ex.export_fragments(subset_dataset, groupby=biorep, selections=set(["rep1"]), ids=None, out_dir=out_path, prefix=prefix+'.', suffix='.bed.gz')
		snap.ex.export_fragments(subset_dataset, groupby=biorep, selections=set(["rep2"]), ids=None, out_dir=out_path, prefix=prefix+'.', suffix='.bed.gz')
		if split:
			os.system("python 03a.pseudo_tags.py -i {0}{1}.bed.gz -o {0}{1}".format(out_path, prefix))
		else:
			snap.ex.export_fragments(subset_dataset, groupby="pseudo", selections=set(["pseudo1"]), ids=None, out_dir=out_path, prefix=prefix+'.', suffix='.bed.gz')
			snap.ex.export_fragments(subset_dataset, groupby="pseudo", selections=set(["pseudo2"]), ids=None, out_dir=out_path, prefix=prefix+'.', suffix='.bed.gz')
		os.system("macs2 callpeak -t {0}{1}.bed.gz -f BED -n {1} -g mm -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir {0}".format(out_path, prefix))
		os.system("macs2 callpeak -t {0}{1}.rep1.bed.gz -f BED -n {1}.rep1 -g mm -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir {0}".format(out_path, prefix))
		os.system("macs2 callpeak -t {0}{1}.rep2.bed.gz -f BED -n {1}.rep2 -g mm -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir {0}".format(out_path, prefix))
		os.system("macs2 callpeak -t {0}{1}.pseudo1.bed.gz -f BED -n {1}.pseudo1 -g mm -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir {0}".format(out_path, prefix))
		os.system("macs2 callpeak -t {0}{1}.pseudo2.bed.gz -f BED -n {1}.pseudo2 -g mm -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir {0}".format(out_path, prefix))
	else:
		pool_peaks = snap.tl.call_peaks(subset_dataset, groupby="pool", q_value=0.01, out_dir=out_path, inplace=False)
		os.system("gunzip -c {0}pool.NarrowPeak.gz > {0}{1}_peaks.narrowPeak; mv {0}pool_insertion.bed.gz {0}{1}_insertion.bed.gz; rm {0}pool.NarrowPeak.gz".format(out_path, prefix))
		biorep_peaks = snap.tl.call_peaks(subset_dataset, groupby=biorep, q_value=0.01, out_dir=out_path, inplace=False)
		os.system("gunzip -c {0}rep1.NarrowPeak.gz > {0}{1}.rep1_peaks.narrowPeak; mv {0}rep1_insertion.bed.gz {0}{1}.rep1_insertion.bed.gz; rm {0}rep1.NarrowPeak.gz".format(out_path, prefix))
		os.system("gunzip -c {0}rep2.NarrowPeak.gz > {0}{1}.rep2_peaks.narrowPeak; mv {0}rep2_insertion.bed.gz {0}{1}.rep2_insertion.bed.gz; rm {0}rep2.NarrowPeak.gz".format(out_path, prefix))
		#biorep1_peaks = snap.tl.call_peaks(subset_dataset, groupby="biorep", selections=set(["rep1"]), q_value=0.01, out_dir="/".join([out_path, "biorep1"]), inplace=False)
		#biorep2_peaks = snap.tl.call_peaks(subset_dataset, groupby="biorep", selections=set(["rep2"]), q_value=0.01, out_dir="/".join([out_path, "biorep2"]), inplace=False)
		pseudo_peaks = snap.tl.call_peaks(subset_dataset, groupby="pseudo", q_value=0.01, out_dir=out_path, inplace=False)
		os.system("gunzip -c {0}pseudo1.NarrowPeak.gz > {0}{1}.pseudo1_peaks.narrowPeak; mv {0}pseudo1_insertion.bed.gz {0}{1}.pseudo1_insertion.bed.gz; rm {0}pseudo1.NarrowPeak.gz".format(out_path, prefix))
		os.system("gunzip -c {0}pseudo2.NarrowPeak.gz > {0}{1}.pseudo2_peaks.narrowPeak; mv {0}pseudo2_insertion.bed.gz {0}{1}.pseudo2_insertion.bed.gz; rm {0}pseudo2.NarrowPeak.gz".format(out_path, prefix))
		#pseudo1_peaks = snap.tl.call_peaks(subset_dataset, groupby="pseudo", selections=set(["pseudo1"]), q_value=0.01, out_dir="/".join([out_path, "pseudo1"]), inplace=False)
		#pseudo2_peaks = snap.tl.call_peaks(subset_dataset, groupby="pseudo", selections=set(["pseudo2"]), q_value=0.01, out_dir="/".join([out_path, "pseudo2"]), inplace=False)
		if bigwig:
			snap.ex.export_bigwig(subset_dataset, groupby="pool", selections=None, resolution=1, out_dir=out_path, prefix=prefix, suffix='.bw')
			os.system("mv {0}{1}pool.bw {0}{1}.bw".format(out_path, prefix))
			snap.ex.export_bigwig(subset_dataset, groupby=biorep, selections=set(["rep1"]), resolution=1, out_dir=out_path, prefix=prefix+'.', suffix='.bw')
			snap.ex.export_bigwig(subset_dataset, groupby=biorep, selections=set(["rep2"]), resolution=1, out_dir=out_path, prefix=prefix+'.', suffix='.bw')
			snap.ex.export_bigwig(subset_dataset, groupby="pseudo", selections=set(["pseudo1"]), resolution=1, out_dir=out_path, prefix=prefix+'.', suffix='.bw')
			snap.ex.export_bigwig(subset_dataset, groupby="pseudo", selections=set(["pseudo2"]), resolution=1, out_dir=out_path, prefix=prefix+'.', suffix='.bw')			
	print("Results written to file:", out_path)

	sys.exit(0)

if __name__ == "__main__":
	main()



