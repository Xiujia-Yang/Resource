#!/usr/bin/env python
'''
    This script is for generating positional mutation
    frequency matrix with each row representing an allele
    and each column representing a position.
    
    Usage: python positional_mut_freq_profile.py allelelist indir
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

def mutRateCal(s):
    '''
        This script is to calculate mutation rate using
        merged dF returned by function allele_positional_mutation_summary
    '''
    return (s["A"]+s["T"]+s["C"]+s["G"])/s["nReads"]*100

def allele_positional_mutation_summary(pathList, allele, outfl):
    '''
        This function is to combine individual positional mutation
        statistic files for a single allele and then output a com-
        bined pandas Series object storing percentage of mutation 
        for each positon
    '''
    for i, path in enumerate(pathList):
        df = pd.read_csv(path, sep="\t")
        if i > 0:
            dF[["A", "T", "C", "G", "nReads"]] = dF[["A", "T", "C", "G", "nReads"]] + df[["A", "T", "C", "G", "nReads"]]
        else:
            dF = df
    # output the number of samples containing these allele
    # and the the number of total clones recombined from this
    # allele
    outfl.writerow([allele, i+1, dF["nReads"].tolist()[0]])
    # calculate positional mutation rate using apply function
    pos_mutRate_series = dF.apply(mutRateCal, axis=1, result_type="reduce")
    # give the resulting series a name (allele ID)
    pos_mutRate_series.name = allele
    return pos_mutRate_series

def main(outfl):
    # read the region delimitation file
    df_delim = pd.read_csv("data/v.allele.region.delimitation.txt", sep="\t")
    df_delim.vAllele = df_delim.vAllele.apply(lambda x:x.replace('/', '.').replace('*', '.'))
    # read the allele list, of which we want to explore its possiblity
    # alleles were assumed to be ordered
    allele_fl_handle = csv.reader(open(allele_list, "rU"), delimiter="\t")

    # iterated over each allele, merge mutation statistics
    # and then merge the combined results into a DataFrame
    # with each column represents an allele and each row 
    # represents a position. This dataframe can be directly 
    # applied for visualization
    n = 0
    for [allele] in allele_fl_handle:
        if allele in df_delim.vAllele.tolist():
            # obtain the begining index and length of each functional region
            fr1, cdr1, fr2, cdr2, fr3, cdr3 = df_delim[df_delim.vAllele==allele].ix[:,range(1,7)].T.squeeze().to_list()
            fr1_len = cdr1-fr1; cdr1_len = fr2-cdr1; fr2_len = cdr2-fr2; cdr2_len = fr3-cdr2; fr3_len = cdr3-fr3
            # obtain all allele statistic file
            pathList = glob.glob(r'%s/*-IGHG/%s.pos.mut.type.stat'%(indir, allele))
            if pathList != []:
                print (allele)
                pos_mutRate_series = allele_positional_mutation_summary(pathList, allele, outfl)
                if n > 0:
                    if fr1_len != 0:
                        s_fr1 = pos_mutRate_series[fr1:cdr1].reset_index(drop=True)
                        dF_fr1 = pd.concat([dF_fr1, s_fr1], axis=1)
                    if cdr1_len != 0:
                        s_cdr1 = pos_mutRate_series[cdr1:fr2].reset_index(drop=True)
                        dF_cdr1 = pd.concat([dF_cdr1, s_cdr1], axis=1)
                    if fr2_len != 0:
                        s_fr2 = pos_mutRate_series[fr2:cdr2].reset_index(drop=True)
                        dF_fr2 = pd.concat([dF_fr2, s_fr2], axis=1)
                    if cdr2_len != 0:
                        s_cdr2 = pos_mutRate_series[cdr2:fr3].reset_index(drop=True)
                        dF_cdr2 = pd.concat([dF_cdr2, s_cdr2], axis=1)
                    if fr3_len != 0:
                        s_fr3 = pos_mutRate_series[fr3:cdr3].reset_index(drop=True)
                        dF_fr3 = pd.concat([dF_fr3, s_fr3], axis=1)
                else:
                    df = pos_mutRate_series.to_frame()
                    if fr1_len != 0:
                        dF_fr1 = df[fr1:cdr1].reset_index(drop=True)
                    if cdr1_len != 0:
                        dF_cdr1 = df[cdr1:fr2].reset_index(drop=True)
                    if fr2_len != 0:
                        dF_fr2 = df[fr2:cdr2].reset_index(drop=True)
                    if cdr2_len != 0:
                        dF_cdr2 = df[cdr2:fr3].reset_index(drop=True)
                    if fr3_len != 0:
                        dF_fr3 = df[fr3:cdr3].reset_index(drop=True)
                n = n+1
    # concat different fuctional region
    df = pd.concat([dF_fr1,dF_cdr1,dF_fr2,dF_cdr2,dF_fr3], join='outer')
    df = df.reset_index(drop=True)
    # output the calcualted statistics
    df.T.to_csv("positional_mut_freq_matrix.txt"%allele_list, sep="\t")
    # output the length of each functional region
    outfl2.writerow(["FR1", dF_fr1.shape[0]])
    outfl2.writerow(["CDR1", dF_cdr1.shape[0]])
    outfl2.writerow(["FR2", dF_fr2.shape[0]])
    outfl2.writerow(["CDR2", dF_cdr2.shape[0]])
    outfl2.writerow(["FR3", dF_fr3.shape[0]])

if __name__ == "__main__":
    allele_list = sys.argv[1]
    indir = sys.argv[2]
    # open a file to write the number of sample and clone containing this allele
    outfl = csv.writer(open("allele_sample_clone_number.txt"%allele_list, "w"), delimiter="\t")
    outfl.writerow(["allele", "sample", "clone"])
    #
    outfl2 = csv.writer(open("region_length_for_matrix_anno.txt"%allele_list, "w"), delimiter="\t")
    outfl2.writerow(["region", "length"])
    main(outfl)
