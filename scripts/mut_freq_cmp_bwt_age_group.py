#!/usr/bin/env python
'''
    This script is for comparing regional mutation
    frequency between different age group.

    Usage: python mut_freq_cmp_bwt_age_group.py mut_comb_fl
'''

import sys, csv, re, os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
   
    # read the statistic file
    df = pd.read_csv(mut_comb_fl, sep="\t", index_col=0)
    
    # make general filter
    df_filter2 = df[df.tissue=="PBMC"]
    df_filter2 = df_filter2[df_filter2.productivity=="productive"]
    df_filter2 = df_filter2[df_filter2.disease=="Healthy"]
    df_filter2 = df_filter2[df_filter2.age!="Unknown"]
    df_filter2 = df_filter2[df_filter2.nProdReads>1]
    # delete variable to save disk space
    del df

    ######################### make different kinds of plots ##################################
    #  mutation rate between different different gender (separated by isotype)
    age_order = ["u 1", "1-14", "15-24", "25-44", "45-64", "65+", "Unknown"]
    
    for iso in ["IGHG", "IGHM", "IGHA", "IGHE", "IGHD"]:
        df_filter = df_filter2[df_filter2.isotype.str.contains(iso, na=False)]
        if df_filter.size != 0:
            df_filter.mutRate = df_filter.mutRate.astype(float)
            order = [x for x in age_order if x in list(df_filter.age.unique())]
            plt.figure()
            ax = sns.boxplot(x="region", y="mutRate",
                hue="age", palette=sns.color_palette("colorblind"),
                data=df_filter, hue_order=order)
            plt.ylim(0,0.7)
            plt.savefig("mut_freq_age_cmp_%s.png"%iso, dpi=300)
    
if __name__ == "__main__":
    mut_comb_fl = sys.argv[1]
    main()
