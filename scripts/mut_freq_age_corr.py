#!/usr/bin/env python
'''
    This script is for exploring the correlation between 
    age and mutation rate.
    
    Usage: python mut_freq_age_corr.py mut_comb_fl
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

def main():
    # read the statistic file
    df = pd.read_csv(mut_comb_fl, sep="\t", index_col=0)
    
    # filter to obtain qualified records
    df = df[df.productivity=="productive"]
    df = df[df.tissue=="PBMC"]
    df = df[df.disease=="Healthy"]
    df = df[df.isotype.notnull()]
    df = df[df.original_age!="Unknown"]
    df = df[df.nProdReads>1]
            
    # add new column for isotype for group
    df["iso"] = df.isotype.apply(lambda x:x[:4])
    
    # change dtype to circumvent error
    df.mutRate = df.mutRate.astype(float)
    df.original_age = df.original_age.astype("int")

    for iso, group in df.groupby(df.iso):
        # construct a new dataframe used as input of sns.lmplot
        mutRate = group.mutRate.groupby(group.index).describe()["mean"]
        group["index"] = group.index
        gender = group[["gender","original_age", "index"]].drop_duplicates()["gender"]
        age = group[["gender","original_age", "index"]].drop_duplicates()["original_age"]
        group = pd.concat([mutRate, gender, age], axis=1)
        group.rename(columns={"mean":"mutRate"}, inplace=True)
        
        # calculated regression model parameter
        for gender in list(group.gender.unique()):
            regr = linear_model.LinearRegression()
            regr.fit(group.original_age[group.gender==gender].values.reshape(-1,1), group.mutRate[group.gender==gender])
            [coef], intercept = regr.coef_, regr.intercept_
            r2 =regr.score(group.original_age[group.gender==gender].values.reshape(-1,1), group.mutRate[group.gender==gender])
            outfl.writerow([iso, gender, coef, intercept, r2])
        
        # make lmplot
        plt.figure() # open a new figure
        sns.lmplot(x="original_age", y="mutRate", hue="gender", data=group, 
                    height=3, aspect=1.33, fit_reg=True, scatter_kws={"s": 5})
        sns.despine(top=False, right=False)
        plt.savefig("mut_freq_age_regression_%s.png"%iso, dpi=300)
    
if __name__ == "__main__":
    mut_comb_fl = sys.argv[1]
    outfl = csv.writer(open("mut_freq_age_regression_model.txt", "w"), delimiter="\t")
    outfl.writerow(["isotype", "gender", "coef", "intercept", "R-squared"])
    main()
