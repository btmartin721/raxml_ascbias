#!/usr/bin/env python3

import argparse
import time
import sys

# To count invariant sites for Stamatakis correction
from collections import Counter

# Creates dataframe from PHYLIP input file; need pandas library
import pandas as pd
import numpy as np

# For reading PHYLIP input file; need biopython library
from Bio import AlignIO

def Get_Arguments():

    parser = argparse.ArgumentParser(description="Does ascertainment bias correction for RAxML and does Felsenstein and Stamatakis counts")

    parser.add_argument("-p", "--phylip", type=str, required=True, help="Input PHYLIP filename")
    parser.add_argument("-o", "--outfile", type=str, required=False,
                        help="Output filename; invariant site count filenames will append .felsenstein and .stamatakis; Default = out.phy",
                        nargs="?", default="out.phy")

    args = parser.parse_args()

    return args

# Uses AlignIO to read input PHYLIP file
def Read_Alignment(infile):

    my_id_list = []
    with open(infile) as fin:
        # from Bio import AlignIO
        alignment = AlignIO.read(fin, "phylip-relaxed")
        for record in alignment:
            id = record.id
            seq = record.seq
            my_id_list.append(id)

        matrix = [[char for char in seq] for seq in alignment]  # 2d list

    return matrix, my_id_list

# Identifies and drops invariant columns from pandas DataFrame
def filter_invariants(dframe):

    bases = ["A","G","C","T"]

    # collections::Counter library
    stamatakis_cnt = Counter()
    fels_cnt = 0

    invariant_lst = list()

    # Loop through each dataframe column
    for i in dframe.columns:

        # Gets unique values at each column and saves to list
        column_unique = dframe[i].unique().tolist()

        # Intersects column_unique with bases list
        intersect = [value for value in bases if value in column_unique]

        # If column contains only ambigous or IUPAC characters
        # Save the column index for dropping later
        if not any(value for value in bases if value in column_unique):
            invariant_lst.append(i)

        # If site is invariant (only A, C, G, or T); ignores N's and "-"
        if len(intersect) == 1:
            # Uses collections::Counter to get Stamatakis counts
            stamatakis_cnt[intersect[0]] += 1

            # Counts number of invariant sites for Felsenstein count
            fels_cnt += 1

            # Saves column indexes to list
            invariant_lst.append(i)

    # Drops invariant sites from dataframe
    dframe.drop(invariant_lst, axis=1, inplace=True)

    return stamatakis_cnt, fels_cnt, dframe

# Writes three output files: *.phy, *.phy.stamatakis, *.phy.felsenstein
def write_output(dframe, outfile, ids, st, fel):

    write_phylip(dframe, outfile, ids) # Write DataFrame to PHYLIP outfile

    felsenstein = (outfile + ".felsenstein")
    with open(felsenstein, "w") as fout:
        fout.write(str(fel)) # Writes number of invariant sites to outfile

    stamatakis = (outfile + ".stamatakis")
    with open(stamatakis, "w") as fout:
        # Write number of invariant sites containing A C G T to outfile
        fout.write(str(st["A"]) + " " + str(st["C"]) + " " + str(st["G"]) + " " + str(st["T"]) + "\n")

# Writes only variant sites to output file
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

######################################MAIN######################################################################

start = time.time() # time library

arguments = Get_Arguments() # argparse library

data, ids = Read_Alignment(arguments.phylip) # Reads PHYLIP file using biopython's AlignIO

df = pd.DataFrame(data, ids) # Creates pandas DataFrame

# For each column of pandas DataFrame
# Drops column if it is invariant
# Counts number of invariant sites (Felsenstein)
# And counts number of invariant sites containing A C G and T (Stamatakis)
stam, fels, df = filter_invariants(df)

write_output(df, arguments.outfile, ids, stam, fels)

# Prints execution time for script
end = time.time()
delay = (end - start)
print("\nExecution time: {} seconds\n".format(round(delay, 2)))
