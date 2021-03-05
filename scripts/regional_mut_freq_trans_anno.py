#!/usr/bin/env python
'''
    This script is for transforming the regional mutation
    statistic file into a format can be accepted by seaborn
    and annotating each sample.

    Usage: python regional_mut_freq_trans_anno.py sample statfl outdir
'''
import sys, csv, re, os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import linear_model
from itertools import combinations

def main():
    # Read the metadata information
    meta="data/metadata-20191219-ZZHLAB-original-name.csv"
    df_meta = pd.read_csv(meta, 
                          usecols=["Run","ageGroup","gender-classification","tissue-classification",
                                  "disease-classification"],
                         index_col=0)
    df_meta = df_meta.rename(columns={"ageGroup":"age","gender-classification":"gender","tissue-classification":"tissue",
                                   "disease-classification":"disease"})
    
    # define a dataframe used as input for seaborn
    df_sns = pd.DataFrame({"cloneId":[], "isotype":[], "productivity":[], "region":[], "mutRate":[], "nAvaiReads":[]})

    # read the subfile and concatenate them, and then add the annotation info
    # e.g. gender, disease, age, tissue
    df = pd.read_csv(statfl, sep="\t")
    df_all = df[["cloneId","isotype","productivity","nAvaiReads","allMutRate","nProdReads","nUnprodReads"]]
    df_all["region"] = "allMutRate"
    df_all.rename(columns={"allMutRate":"mutRate"}, inplace=True)

    df_fr1 = df[["cloneId","isotype","productivity","nAvaiReads","fr1MutRate","nProdReads","nUnprodReads"]]
    df_fr1["region"] = "fr1MutRate"
    df_fr1.rename(columns={"fr1MutRate":"mutRate"}, inplace=True)

    df_cdr1 = df[["cloneId","isotype","productivity","nAvaiReads","cdr1MutRate","nProdReads","nUnprodReads"]]
    df_cdr1["region"] = "cdr1MutRate"
    df_cdr1.rename(columns={"cdr1MutRate":"mutRate"}, inplace=True)

    df_fr2 = df[["cloneId","isotype","productivity","nAvaiReads","fr2MutRate","nProdReads","nUnprodReads"]]
    df_fr2["region"] = "fr2MutRate"
    df_fr2.rename(columns={"fr2MutRate":"mutRate"}, inplace=True)

    df_cdr2 = df[["cloneId","isotype","productivity","nAvaiReads","cdr2MutRate","nProdReads","nUnprodReads"]]
    df_cdr2["region"] = "cdr2MutRate"
    df_cdr2.rename(columns={"cdr2MutRate":"mutRate"}, inplace=True)

    df_fr3 = df[["cloneId","isotype","productivity","nAvaiReads","fr3MutRate","nProdReads","nUnprodReads"]]
    df_fr3["region"] = "fr3MutRate"
    df_fr3.rename(columns={"fr3MutRate":"mutRate"}, inplace=True)

    df_temp = pd.concat([df_all, df_fr1, df_cdr1, df_fr2, df_cdr2, df_fr3])
    df_temp.index = [sample]*df_temp.shape[0]
    df_sns = pd.concat([df_sns, df_temp])
    
    # concatenate statistic dataframe with metadata for annotation
    df_sns = pd.concat([df_sns, df_meta], join="inner", axis=1)
    df_sns.to_csv("%s/%s_mut_freq_trans_anno.txt"%(outdir, sample), sep='\t')
    
if __name__ == "__main__":
    sample = sys.argv[1]
    statfl = sys.argv[2]
    outdir = sys.argv[3]
    main()
