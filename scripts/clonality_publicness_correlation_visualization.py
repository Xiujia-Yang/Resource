#!/usr/bin/env python

'''
    This script is for visualizing the relation between
    clonality and publicness.
    
    Usage: python clonality_publicness_correlation_visualization.py pub_record.tab metadata.tab
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
    # import the public clone records file 
    df = pd.read_csv(pub_record_fl, sep="\t",index_col=0)
    # calculated the percentage of each public clone
    df["pct"] = np.log10(df.cloneCount/df.nProdReads*100)
    # get clone tag
    df["cloneTag"] = df.allVHitsWithScore.str.split("*", expand=True)[0] + "-" + \
    df.allJHitsWithScore.str.split("*", expand=True)[0] + "-" + \
    df.targetSequences
    # Count cloneTag and assign the count to each cloneTag
    df_count = pd.DataFrame({"index":df.cloneTag.value_counts().index, "count":df.cloneTag.value_counts().tolist()})
    df_merge =  pd.merge(df, df_count, left_on=["cloneTag"], right_on=["index"])
    
    # Count cloneTag to obtain the number of donors and number of projects
    df_meta = pd.read_csv(metadata_fl, index_col=0)
    p = df_meta.BioProject.to_dict()
    d = df_meta.Donor.to_dict()
    df_merge["project"] = [p[x] for x in df_merge["sample"].tolist()]
    df_merge["donor"] = [d[x] for x in df_merge["sample"].tolist()]
    
    # Group by cloneTag to got the number of projects and donors share a certain cloneTag
    pn = {}; dn={}
    for [g, ps], [_, ds] in zip(df_merge.groupby("cloneTag")["project"], df_merge.groupby("cloneTag")["donor"]):
        pn.update({g:len(set(ps))})
        dn.update({g:len(set(ds))})
    df_merge["pn"] = [pn[x] for x in df_merge["cloneTag"].tolist()]
    df_merge["dn"] = [dn[x] for x in df_merge["cloneTag"].tolist()]
   
    # drop exception caused by mixcr's inconsistency between target sequences and cdr3 seqeucnes
    df_merge = df_merge[df_merge.pn!=1]
    df_merge = df_merge[df_merge.dn!=1]
 
    # obtain for each cloneTag the CDR3 length and insert it into the dataframe
    df_merge["cdr3len"] = df_merge.targetSequences.apply(lambda x:len(x))
    sdata = pd.DataFrame({"cdr3len":df_merge.cdr3len.groupby(df_merge["count"]).describe()["mean"].tolist(), "count":list(df_merge["count"].unique())})
    ddata = pd.DataFrame({"cdr3len":df_merge.cdr3len.groupby(df_merge["dn"]).describe()["mean"].tolist(), "count":list(df_merge["dn"].unique())})
    pdata = pd.DataFrame({"cdr3len":df_merge.cdr3len.groupby(df_merge["pn"]).describe()["mean"].tolist(), "count":list(df_merge["pn"].unique())})

    # visualization module
    plt.figure(figsize=(10,3))
    ax = sns.boxplot(x="dn", y="pct", data=df_merge)
    ax.set_xlabel("Donor sharing count")
    ax.set_ylabel("Normalized clone size (log10)")
    ax2 = ax.twinx()
    ax2.plot(ddata.sort_values("count").cdr3len.tolist())
    ax2.set_ylabel("Average CDR3 length (bp)")
    plt.savefig("Clonality_against_publicness_and_CDR3_length.png", dpi=300, order=df_merge["dn"].sort_values().unique(), bbox_inches="tight")
   
    
if __name__ == "__main__":
    pub_record_fl = sys.argv[1]
    metadata_fl = sys.argv[2]
    main()
