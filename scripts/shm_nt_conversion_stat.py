#!/usr/bin/env python
'''
    This script is for calculating statistic of 
    the frequency of nucleotide transition
    
    Usage: python nt_conversion_stat.py sampleid_file isotype indir outdir
'''
import sys, csv, re, os
import glob
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import linear_model
from itertools import combinations

def is_nonsyn(series):
    '''
        Discern the locus type
    '''
    indexName = ['is_syn'+nt for nt in ['A', 'T', 'C', 'G'] if nt != series['ref']]
    return True if series[indexName].tolist() == ['No', 'No', 'No'] else False
        
    
def main(isotype):
    ################### Initialization module
    # define variable
    total_dict = {"A":0, "C":0, "G":0, "T":0}
    conv_dict = {}
    df_conversion = pd.DataFrame({
        'A':[0,0,0,0],
        'C':[0,0,0,0],
        'G':[0,0,0,0],
        'T':[0,0,0,0],
    }, index=["A","C","G","T"]) # used as input for seaborn
    
    ################### Perform statistic analysis module
    # obtain the path for all allele position mutation statistics file
    path_list = []
    for [sample] in sampleHandle:
        path = glob.glob(r'%s/%s*IGH%s/*.pos.mut.type.stat'%(indir,sample,isotype))
        if path != []:
            path_list = path_list + path
    # iterated over all stat fl to perform conversion statistics
    for path in path_list:
        # read the positional mutation statistics file
        df = pd.read_csv(path, sep="\t", index_col=0)
        nReads = df.loc[0, "nReads"]
        # considering different type of locus independently
        df = df[df.apply(is_nonsyn, axis=1)]
        # add the number of each nucleotide as the background
        # number for later normalization
        nt_list = df["ref"].tolist()
        for nt in nt_list:
            total_dict[nt] += nReads
        # statistics of the conversion number of each pair nt
        for nt, a, c, g, t in zip(nt_list, df["A"].tolist(), df["C"].tolist(), df["G"].tolist(), df["T"].tolist()):
            try:
                conv_dict[nt+'-'+"A"] += a
            except:
                conv_dict[nt+'-'+"A"] = a
            try:
                conv_dict[nt+'-'+"C"] += c
            except:
                conv_dict[nt+'-'+"C"] = c
            try:
                conv_dict[nt+'-'+"G"] += g
            except:
                conv_dict[nt+'-'+"G"] = g
            try:
                conv_dict[nt+'-'+"T"] += t
            except:
                conv_dict[nt+'-'+"T"] = t
    
    ################### Postprocessing module
    # assign the value in the conversion dictionary
    # to dataframe
    for from_nt in ["A", "C", "G", "T"]:
        for to_nt in ["A", "C", "G", "T"]:
            df_conversion.loc[to_nt, from_nt] = conv_dict[from_nt+'-'+to_nt]
    # normalize the dataframe
    df_conversion1 = df_conversion
    for from_nt in ["A", "C", "G", "T"]:
        df_conversion1[from_nt] = df_conversion1[from_nt]/total_dict[from_nt]*100
    df_conversion1.to_csv("%s/nt_conversion_matrix_%s_%s.txt"%(outdir,isotype,type))
    ################### Visualization module
    plt.figure()
    sns.heatmap(df_conversion1, cmap="Blues", square=True)
    plt.savefig("%s/nt_conversion_matrix_%s_%s.png"%(outdir,isotype,type), dpi=600)
    
if __name__ == "__main__":
    sampleidFl = sys.argv[1]
    sampleHandle = csv.reader(open(sampleidFl, "r"), delimiter="\t")
    isotype = sys.argv[2]
    indir = sys.argv[3]
    outdir = sys.argv[4]
    main(isotype)
