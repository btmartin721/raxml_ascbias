#!/usr/bin/env python3

import argparse
import time

import pandas as pd
import numpy as np

#from itertools import product
from Bio import AlignIO

def Get_Arguments():

    parser = argparse.ArgumentParser(description="Extracts variant sites for use with RAxML ascertainment-bias correction")

    parser.add_argument("-f", "--file", type=str, required=True, help="Input filename")
    parser.add_argument("-o", "--outfile", type=str, required=False,
                        help="Output filename; Default = out.phy", nargs="?", default="out.phy")

    args = parser.parse_args()

    return args

# Uses AlignIO to read input PHYLIP file
def Read_Alignment(infile):

    my_id_list = []
    with open(infile) as fin:
        alignment = AlignIO.read(fin, "phylip-relaxed")
        for record in alignment:
            id = record.id
            seq = record.seq
            my_id_list.append(id)

        matrix = [[char for char in seq] for seq in alignment]  # 2d list

    return matrix, my_id_list

# Drops invariable columns from pandas DataFrame
def filter_invariants(dframe):

    tot = len(dframe.columns+1)

    bases = ["A","G","C","T"]

    LOG_EVERY_N = 1000

    for i in dframe.columns:

        current = i+1
        total = tot

        progress = (current / tot) * 100

        if (i % LOG_EVERY_N) == 0:
            print("Progress: {0:.2f}%".format(round(progress, 2)))

        if not len(check_intersect(bases, dframe[i])) > 1:
            dframe.drop(i, axis=1, inplace=True)

## Unused code; may implement at some point but needs optimization
        #for j in product(*[ambiguity_codes(j) for j in dframe[i].values]):
            #expanded_seq = "".join(j)
            #if expanded_seq == len(expanded_seq) * expanded_seq[0]:
                #dframe.drop(i, axis=1, inplace=True)
                #break

# Dictionary to phase each column in pandas DataFrame
# def ambiguity_codes(char):
#
#     iupac = {
#         'A': ["A"],
#         'G': ["G"],
#         'C': ["C"],
#         'T': ["T"],
#         'N': [""],
#         '-': [""],
#         'Y': ["C", "T"],
#         'R': ["A", "G"],
#         'W': ["A", "T"],
#         'S': ["G", "C"],
#         'K': ["T", "G"],
#         'M': ["C", "A"],
#         'B': ["C", "G", "T"],
#         'D': ["A", "G", "T"],
#         'H': ["A", "C", "T"],
#         'V': ["A", "C", "G"]
#         }
#
#
#     return iupac[char]

def write_phylip(dframe, outfile, ids):
    df_size = dframe.shape

    header = str(df_size[0]) + " " + str(df_size[1])

    seq_lst = dframe.values.tolist()
    sample_lst = dframe.index.tolist()
    joined_seqs = list()

    [joined_seqs.append("".join(i)) for i in seq_lst]

    with open(outfile, "w") as fout:
        fout.write(header + "\n")

        [fout.write("{}\t{:>15}\n".format(str(sample), str(seq)))
        for sample, seq
        in zip(sample_lst, joined_seqs)]

def check_intersect(nt, col):
    return list(set(nt) & set(col))

######################################MAIN######################################################################

start = time.time()

arguments = Get_Arguments() # argparse library

data, ids = Read_Alignment(arguments.file) # Reads PHYLIP file using biopython

df = pd.DataFrame(data, ids) # Creates pandas DataFrame

# For each column of pandas DataFrame
# Drops column if it is invariant
filter_invariants(df)

write_phylip(df, arguments.outfile, ids) # Write DataFrame to PHYLIP outfile

# Prints execution time for script
end = time.time()
delay = (end - start)
print("\nExecution time: {} seconds\n".format(round(delay, 2)))
