#!/usr/bin/env python3

import argparse

import pandas as pd
import numpy as np

from Bio import AlignIO
from itertools import product
from itertools import combinations

def Get_Arguments():

    parser = argparse.ArgumentParser(description="extracts variant sites for use with RAxML ascertainment bias correction")

    parser.add_argument("-f", "--file", type=str, required=True, help="Input filename")
    parser.add_argument("-o", "--outfile", type=str, required=False,
                        help="Output filename; Default = out.txt", nargs="?", default="out.txt")

    args = parser.parse_args()

    return args
    
def Read_Alignment(infile):
    
    my_id_list = []
    with open(infile) as fin:
        alignment = AlignIO.read(fin, "phylip-relaxed")
        for record in alignment:
            id = record.id
            seq = record.seq
            my_id_list.append(id)

        matrix = [[char for char in seq] for seq in alignment]

    return matrix, my_id_list
    
def read_phylip(file):
    
    loci = dict()
    for line in file:
        line = line.rstrip()
        id_seq = line.strip().split()
        #print(id_seq)
        #loci[id_seq[0]] = id_seq[1]
        ids = id_seq[0]
        sequences = id_seq[1]

    yield ids, sequences
    
def drop_invariable_cols(dframe, iupac):
    
    #dframe = dframe.replace(['N'], np.nan, regex=True)
    #dframe = dframe.replace(['-'], np.nan, regex=True)    
    
    df_copy = dframe.copy()
    
    for col in dframe:
        for key, val in iupac.items():
            for item in val:
              
                #if dframe[col].str.contains(key).any():
                df_copy[col] = df_copy[col].replace(key, item)
                print(df_copy)
                unique_cols = df_copy[col].nunique()
                if unique_cols == 1:
                    #print(unique_cols)
                    dframe = dframe.drop(col, axis=1)
                    break
            break
                        

    print(dframe)    

def ambiguity_codes():

    iupac = {
            
            'N': ["A", "G", "C", "T"],
            '-': ["-"],
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
    

   
######################################MAIN######################################################################

arguments = Get_Arguments()

#ids = []
#data = []
data, ids = Read_Alignment(arguments.file)

ambig = ambiguity_codes()


#with open(arguments.file, "r") as fin:
    
    #header = fin.readline()
    #seqs = read_phylip(fin)
        

#for key, val in seqs.items():
#    print(key + "\t" + val)
    

#nt = ["A", "G", "C", "T"]

df = pd.DataFrame(data, ids)
#print(df)

iupac_dict = ambiguity_codes()
#with open(arguments.outfile, "w") as fout:

drop_invariable_cols(df, iupac_dict)
        
        #for id, seq in read_phylip(fin):
            #for i in expand_iupac(seq, ambig):
                #fout.write(str(i))
                
#for i in expanded:
    #print(i)
    

