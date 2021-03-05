from __future__ import division
import re,csv,glob,sys,os
import argparse
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import subprocess
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import seaborn as sns
import math
from matplotlib import gridspec

def get_motif_result(infile, isotype, outdir):
    motif_comp = []
    motif_frac = {}
    order_x = []
    germline_list = []
    for rec in csv.reader(open(infile, 'r'), delimiter="\t"):
        if rec[0] == "Motif" or rec[0] == "#Mut_type":
            continue
        else:
            order_x.append(rec[0])
            motif_comp.append([float(i) for i in rec[4:]])
            if float(rec[1]) == 0:
                motif_frac[rec[0]] = np.log10(float(rec[1]))
            else:
                motif_frac[rec[0]] = float(rec[1])
            germline_list.append(float(rec[3]))

    #order_x = ['WRCY', 'RGYW', 'WA', 'TW', 'SYC', 'GRS', 'A', 'C', 'G', 'T']
    myorder=['A','C','G','T']
    colors=[(0.00784313725490196, 0.6196078431372549, 0.45098039215686275), (0.00392156862745098, 0.45098039215686275, 0.6980392156862745), (0.8705882352941177, 0.5607843137254902, 0.0196078431372549), (0.8352941176470589, 0.3686274509803922, 0.0)]#['green','blue','orange','red']

    bar_dict = {}
    for i,j in enumerate(myorder):
        bar_dict[j]=[n[i] for n in motif_comp]

    line_data = []
    for i in order_x:
        line_data.append(motif_frac[i])

    fig=plt.figure(figsize=(6,5))
    gs=gridspec.GridSpec(2,1,height_ratios=[1,4])
    #fig, ax1 = plt.subplots(figsize=(6, 4))

    ax3 = plt.subplot(gs[0])
    sns.scatterplot(x=order_x, y=[np.log10(i) for i in germline_list], ax=ax3, color="black")
    ax3.plot(order_x,[np.log10(i) for i in germline_list],color="black",linestyle="solid",linewidth=1)
    ax3.set_xlim(-0.5,len(order_x)-0.5,1)
    my_x_ticks = np.arange(1, len(order_x), 1)
    plt.xticks(my_x_ticks)
    plt.xticks(order_x, my_x_ticks, rotation=90, ha="center")
    ax3.set_ylim([4,8])
    ax3.set_xticks([])

    ax1 = plt.subplot(gs[1])
    x = np.arange(len(order_x))
    bottom_arr=np.array([0]*len(order_x))
    for ind, ite in enumerate(myorder):
        if ind == 0:
            ax1.bar(x, bar_dict[ite],color=colors[ind],label=myorder[ind],linewidth=0)
        else:
            ax1.bar(x, bar_dict[ite],bottom=bottom_arr,color=colors[ind],label=myorder[ind],linewidth=0)
        bottom_arr = bar_dict[ite]+bottom_arr

    ax1.set_ylim([0,100])
    #print (line_data)
    ax2 = ax1.twinx()
    #sns.scatterplot(x, line_data, color="black", ax=ax2)
    ax4 = sns.scatterplot(x, line_data, color="black")
    ax2.plot(x,line_data,'black',linewidth=1)
    ax2.set_ylim([0,30])
    ax4.set_ylim([0,30])
    ax2.set_xlim(-0.5,len(order_x)-0.5,1)

    ax1.yaxis.set_ticks_position("right")
    ax2.yaxis.set_ticks_position("left")
    ax4.yaxis.set_ticks_position("left")

    my_x_ticks = np.arange(1, len(order_x), 1)
    plt.xticks(my_x_ticks)
    #plt.xticks(x, order_x, rotation=45, ha="center")
    plt.xticks(x, my_x_ticks, rotation=90, ha="center")
    ax1.spines['left'].set_linewidth(1.5)
    ax1.spines['bottom'].set_linewidth(1.5)
    ax1.spines['right'].set_linewidth(1.5)
    ax1.spines['top'].set_linewidth(1.5)
    plt.rcParams["xtick.major.size"] = 12
    plt.rcParams["ytick.major.size"] = 12
    ax1.set_xticklabels(order_x, rotation=90)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.legend(loc=[1.05, 0], frameon=False)
    plt.savefig("%s/Fig.6-1.%s.motif.2.png"%(outdir,isotype), bbox_inches='tight',dpi=600)
    plt.close()

    return 0



def main():
    get_motif_result(infile, isotype, outdir)



if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Input the information about motif")
    parser.add_argument('-i', '--input', help="Input the file of motif.")
    parser.add_argument('-type', '--isotype', help="Input the isotype of motif.")
    parser.add_argument('-d', '--outdir', default="./", help='Output file directory. Default vaule is current directory.')
    args = parser.parse_args()
    infile = args.input
    isotype = args.isotype
    outdir = args.outdir
    main()
