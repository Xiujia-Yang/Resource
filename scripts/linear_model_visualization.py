#!/usr/bin/env python
'''
    This script is for visualizing the correlation between
    the number of public clones and the product of the
    total clones of the compared samples.

    python linear_model_visualization.py pairwise_rep_comp.tab 
'''

import sys, csv, re, os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
import seaborn as sns


def main(outfl):
    # import the statistic data
    df = pd.read_csv(pairwise_cmp_fl, sep="\t")
    # calculate the product of total clone number for two compared samples
    df["prod"] = df["nCloneA"]*df["nCloneB"]
    
    for threshold in [10000, 100000, 1000000, 2000000, 3000000, 4000000, 5000000]:
        # filter for healthy or diseased sample pairs
        df_filter = df[["prod", "nPubClone"]][(df.nCloneA>=threshold) & (df.nCloneB>=threshold)]
        # calculate regression equation
        regr = linear_model.LinearRegression()
        regr.fit(df_filter['prod'].values.reshape(-1,1), df_filter['nPubClone'])
        [coef], intercept = regr.coef_, regr.intercept_
        r2 =regr.score(df_filter['prod'].values.reshape(-1,1), df_filter['nPubClone'])
        outfl.writerow([threshold, coef, intercept, r2])
        # make scatter plot using seaborn
        plt.figure() # open a new figure
        ax = sns.lmplot(x="prod", y="nPubClone", data=df_filter, 
                    height=3, aspect=1.3, fit_reg=True, scatter_kws={"s": 5})
        sns.despine(top=False, right=False) # show the top and right spine
        plt.xlabel("Product of the number of total clones")
        plt.ylabel("Number of public clones")
        #plt.savefig("scatterPlot.size.gt."+str(threshold)+".png", dpi=300)
        plt.savefig("linear_model_scatter_plot_for_samples_more_than_"+str(threshold)+"_clones.png", dpi=300)

if __name__ == "__main__":
    pairwise_cmp_fl = sys.argv[1]
    fl = open("linear_model_params.txt", "w")
    outfl = csv.writer(fl, delimiter="\t")
    outfl.writerow(["CloneNumFilter", "Coef", "Intercept", "R2"])
    main(outfl)
    fl.close()
