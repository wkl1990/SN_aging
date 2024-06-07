#!/bin/python

# example: python /projects/ps-renlab2/kaw033/WK_aging/scripts/03a.pseudo_tags.py -i /projects/ps-renlab2/kaw033/CEMBAv2/CEMBA_data/snapatac2/peaks/c1/pool/c1_pool1_head12.bed.gz -o /projects/ps-renlab2/kaw033/CEMBAv2/CEMBA_data/snapatac2/peaks/c1/pool/c1_head12

import argparse

parser = argparse.ArgumentParser(description="split tn5 file to two pseudo replications")
parser.add_argument("-i", "--input", type=str, dest="input", required=True, help="input h5ad file")
#parser.add_argument("-p", "--processes", type=int, dest="processes", help="number of processes to use")
parser.add_argument("-o", "--output", type=str, dest="output", help="prefix of output")

args = parser.parse_args()

import gzip
import random

def main():
    tn5_file = args.input
    out_prefix = args.output
    #Read file!
    with gzip.open(tn5_file, 'rt') as f:
        lines = f.readlines()
    #Shuffle the lines while keeping each pair together
    lines_index = [i for i in range(0, len(lines), 2)]
    random.shuffle(lines_index)
    #Calculate the split point
    num_lines = len(lines_index)
    split_point = num_lines // 2
    #Split the shuffled lines into two equal parts
    lines_index1 = lines_index[:split_point]
    lines_index2 = lines_index[split_point:]
    shuffled_lines1 = [lines[i] + lines[i+1] for i in lines_index1]
    shuffled_lines2 = [lines[i] + lines[i+1] for i in lines_index2]

    # Write the two parts to gzip compressed files
    output_file1 = ".".join([out_prefix, "pseudo1.bed.gz"])
    output_file2 = ".".join([out_prefix, "pseudo2.bed.gz"])
    with gzip.open(output_file1, 'wt') as f1:
        f1.writelines(shuffled_lines1)
    with gzip.open(output_file2, 'wt') as f2:
        f2.writelines(shuffled_lines2)

if __name__ == "__main__":
    main()
