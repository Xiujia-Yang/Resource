from __future__ import division
import csv, sys, os, glob, re
import numpy as np
import subprocess
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse
from matplotlib import cm
from matplotlib import gridspec
from multiprocessing import Pool

def CountHeatmapProfileData(profile_inf, annoline_inf):
    mydf = pd.read_csv(profile_inf, sep="\t")
    myfamily=mydf["Family"]
    df = mydf.sort_values(by=["Clone_number"], ascending=False)
    #print (df.head(5))
    family=df.pop("Family")
    Germline_id=df.pop("Germline_id")
    Clone_num=df.pop("Clone_number")
    mycolor_bar = dict(zip(myfamily.unique(), ["#FF0000", "#0000FF", "#4CA64C", "#B452CD","#FF8C00","#1E90FF","#5A5A5A"]))
    #mycolor_bar = dict(zip(family.unique(), ["255, 0, 0","0, 0, 255","76, 166, 76","180, 82, 205","255, 140, 0","64, 224, 208","90, 90, 90"]))
    col_colors = family.map(mycolor_bar)
    profile_data=df
    clone_num=Clone_num.tolist()

    pos_avg_dict ={}
    for record in csv.reader(open(annoline_inf, 'r'), delimiter="\t"):
        if record[0] == "Family":
            continue
        else:
            pos_avg_dict[record[0]] = [float(i) for i in record[1:]]

    return profile_data, col_colors, clone_num, pos_avg_dict

def PlotComplexHeatmap(profile_inf, annoline_inf, isotype, outdir):
    profile_data, col_colors, clone_num, pos_avg_dict = CountHeatmapProfileData(profile_inf, annoline_inf)

    fig = plt.figure(figsize=(9, 9))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 8], width_ratios=[0.2, 8, 1])
    ax = plt.subplot(gs[1]) #top_annotation
    y_all=np.vstack([pos_avg_dict["All"]])
    x_all = np.arange(len(pos_avg_dict["All"]))
    ax.stackplot(x_all, y_all, colors="grey",alpha=0.7)
    ax.vlines([75,105,156,186],0,60,colors='k',linestyle="solid",linewidth=0.8)
    ax.set_xlim(-0.5,len(x_all)-0.5,1)
    ax.set_ylim(0,60)
    ax.set_yticks([30, 60])
    ax.set_xticks([])

    ax2=plt.subplot(gs[4])
    use_col = profile_data.columns[0:].to_list()
    profile_data[profile_data[use_col] == -100] = np.nan
    background = profile_data[profile_data[use_col] < 0].fillna(0)
    sns.heatmap(background[use_col], cmap="Blues", alpha=2, ax=ax2, linewidths=0,cbar=False)
    sns.heatmap(profile_data, cmap="Reds",alpha=0.6, ax=ax2, cbar=False,linewidths=0)
    ax2.set_yticks([])
    ax2.set_xticks([])

    ax3=plt.subplot(gs[5])
    y3=np.array([np.log10(i) for i in clone_num[::-1]])
    x3=np.arange(len(clone_num))
    ax3.barh(x3, y3, color="grey", align='center',height=1,linewidth=0)
    ax3.set_xlim(0,5)
    ax3.set_xticks([2, 4])
    ax3.set_ylim(-0.5, len(clone_num)-0.5,1)
    ax3.set_yticks([])

    ax4=plt.subplot(gs[3])
    x4=np.arange(len(clone_num))
    y4=[1]*len(clone_num)
    ax4.barh(x4, y4, color=col_colors[::-1], align='center',height =1, linewidth=0)
    ax4.set_ylim(-0.5,len(clone_num)-0.5,1)
    ax4.set_xlim(0,0.8)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.spines["top"].set_alpha(0)
    ax4.spines["bottom"].set_alpha(0)
    ax4.spines["left"].set_alpha(0)
    ax4.spines["right"].set_alpha(0)


    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig("%s/Fig.5-2.%s.Heatmap_profile.png"%(outdir, isotype),dpi=600)
    plt.close()

    return 0



def main():
    PlotComplexHeatmap(profile, anno, isotype, outdir)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Input the results about mutational profile")
    parser.add_argument('-p', '--profile', help="Input the file of mutational profile.")
    parser.add_argument('-a', '--anno', help="Input the file of annotation profile.")
    parser.add_argument('-type', '--isotype', help="Input the isotype of substitution.")
    parser.add_argument('-d', '--outdir', default="./", help='Output file directory. Default vaule is current directory.')
    args = parser.parse_args()
    profile = args.profile
    anno = args.anno
    outdir = args.outdir
    isotype = args.isotype
    main()
