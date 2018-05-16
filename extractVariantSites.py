#!/usr/bin/env python3

import argparse

import pandas as pd
import numpy as np

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
def drop_invariable_cols(dframe, iupac):

    df_copy = dframe.copy()

    for key, val in iupac.items():
        for item in val:
            for col in dframe.columns:

                df_copy[col] = df_copy[col].replace(key, item) # Replace iupac char with all possibilities

                unique_cols = df_copy[col].nunique()

                if unique_cols == 1: # If site is monomorphic after phasing drop column
                    dframe.drop(col, axis=1, inplace=True)

    return dframe # modified DataFrame

# Dictionary to phase each column in pandas DataFrame
def ambiguity_codes():

    iupac = {

            'N': ["A", "G", "C", "T"],
            '-': ["A", "G", "C", "T"],
            'Y': ["C", "T"],
            'R': ["A", "G"],
            'W': ["A", "T"],
            'S': ["G", "C"],
            'K': ["T", "G"],
            'M': ["C", "A"],
            'B': ["C", "G", "T"],
            'D': ["A", "G", "T"],
            'H': ["A", "C", "T"],
            'V': ["A", "C", "G"]
            }


    return iupac

def write_phylip(dframe, outfile, ids):
    df_size = dframe.shape

    header = str(df_size[0]) + " " + str(df_size[1])

    seq_lst = dframe.values.tolist()
    sample_lst = dframe.index.tolist()
    joined_seqs = list()

    for i in seq_lst:
        joined_seqs.append("".join(i))

    with open(outfile, "w") as fout:
        fout.write(header + "\n")
        for sample, seq in zip(sample_lst, joined_seqs):
            fout.write("{}\t{:>15}\n".format(str(sample), str(seq)))



######################################MAIN######################################################################

arguments = Get_Arguments()

data, ids = Read_Alignment(arguments.file) # Reads PHYLIP file

ambig = ambiguity_codes()

df = pd.DataFrame(data, ids) # Creates pandas DataFrame

iupac_dict = ambiguity_codes()

df = drop_invariable_cols(df, iupac_dict)

write_phylip(df, arguments.outfile, ids)
