#!/usr/bin/env python3

import argparse
import time
import sys

import pandas as pd
import numpy as np

#from itertools import product
from Bio import AlignIO

def Get_Arguments():

    parser = argparse.ArgumentParser(description="Does Lewis, Felsenstein, and Stamatkis ascertainment bias corrections for RAxML")

    parser.add_argument("-p", "--file", type=str, required=True, help="Input filename")
    parser.add_argument("-o", "--outfile", type=str, required=False,
                        help="Output file prefix; Default = out", nargs="?", default="out")
    parser.add_argument("-l", "--lewis", action="store_true",
                        help="Boolean; Specifies Lewis correction to remove all invariant sites; default=True")

    parser.add_argument("-f", "--felsenstein", action="store_true",
                        help="Boolean; Specifies Felsenstein correction to count invariant sites; default=False")
    parser.add_argument("-s", "--stamatkis", action="store_true",
                        help="Boolean; Specifies Stamatkis correction to count invarant sites with A C G and T; default=False")

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
def filter_invariants(dframe, lew, fels, stam, outfile):

    tot = len(dframe.columns+1)

    bases = ["A","G","C","T"]

    total_inv_sites = 0
    A = 0
    C = 0
    G = 0
    T = 0

    LOG_EVERY_N = 2000

    for i in dframe.columns:

        current = i+1
        total = tot
        progress = (current / tot) * 100

        if (i % LOG_EVERY_N) == 0:
            print("Progress: {0:.2f}%".format(round(progress, 2)))

        if not len(check_intersect(bases, dframe[i])) > 1 and not lew and not fels and not stam:
            dframe.drop(i, axis=1, inplace=True)
        elif not len(check_intersect(bases, dframe[i])) > 1 and lew:
            dframe.drop(i, axis=1, inplace=True)
        elif not len(check_intersect(bases, dframe[i])) > 1 and fels:
            total_inv_sites+=1
        elif not len(check_intersect(bases, dframe[i])) > 1 and stam:
            As,Cs,Gs,Ts = stamatkis_correction(dframe, i)
            if As > 0:
                A += 1
            elif Cs > 0:
                C += 1
            elif Gs > 0:
                G += 1
            elif Ts > 0:
                T +=1

    write_output(dframe, outfile, ids, total_inv_sites, A, C, G, T)

def write_output(dframe, outfile, ids, total, A, C, G, T):

    if not arguments.lewis and not arguments.felsenstein and not arguments.stamatkis:
        filename = (outfile + ".phy")
        write_phylip(dframe, filename, ids) # Write DataFrame to PHYLIP outfile

    elif arguments.lewis:
        filename = (outfile + ".phy")
        write_phylip(dframe, filename, ids) # Write DataFrame to PHYLIP outfile

    elif arguments.felsenstein:
        filename = outfile + ".felsenstein"
        with open(filename, "w") as fout:
            fout.write(str(total)) # Writes number of invariant sites to outfile

    elif arguments.stamatkis:
        filename = (outfile + ".stamatkis")
        with open(filename, "w") as fout:
            # Write number of invariant sites containing A C G T to outfile
            fout.write(str(A) + " " + str(C) + " " + str(G) + " " + str(T) + "\n")

# For Lewis correction. Writes only variant sites to output file
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

def stamatkis_correction(dframe, col):
    A = dframe[col].str.count("A").sum()
    C = dframe[col].str.count("C").sum()
    G = dframe[col].str.count("G").sum()
    T = dframe[col].str.count("T").sum()

    return A,C,G,T

######################################MAIN######################################################################

start = time.time() # time library

arguments = Get_Arguments() # argparse library

# Checks to make sure only one correction opion is used
if arguments.lewis and arguments.stamatkis:
    sys.exit("Error: Cannot use more than one type of correction")
elif arguments.lewis and arguments.felsenstein:
    sys.exit("Error: Cannot use more than one type of correction")
elif arguments.felsenstein and arguments.stamatkis:
    sys.exit("Error: Cannot use more than one type of correction")
elif arguments.lewis and arguments.felsenstein and arguments.stamatkis:
    sys.exit("Error: Cannot use more than one type of correction")

data, ids = Read_Alignment(arguments.file) # Reads PHYLIP file using biopython's AlignIO

df = pd.DataFrame(data, ids) # Creates pandas DataFrame

# For each column of pandas DataFrame
# Drops column if it is invariant (lewis)
# Counts number of invariant sites (Felsenstein)
# Or Count number of invariant sites containing A C G and T (Stamatkis)
filter_invariants(df, arguments.lewis, arguments.felsenstein, arguments.stamatkis, arguments.outfile)

# Prints execution time for script
end = time.time()
delay = (end - start)
print("\nExecution time: {} seconds\n".format(round(delay, 2)))
