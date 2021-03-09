import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt
import csv, os, sys
from pandas import DataFrame
import warnings
from matplotlib import gridspec

def DDHeatmap(infile, Dlocal, output):
	DDtotaldf = pd.read_csv(infile, sep = "\t")
	Dlocaldf = pd.read_csv(Dlocal, header = None)
	Dlocaldf.columns = ["Dposition"]
	Dlist = Dlocaldf["Dposition"]
	DDtotaldf.dropna(how = 'any', inplace = True)
	DDtotaldf["D1_len"] = DDtotaldf['D1_seq'].apply(len)
	DDtotaldf["D2_len"] = DDtotaldf["D2_seq"].apply(len)
	DDtotaluse = DDtotaldf[(DDtotaldf["D1_len"] >= 12)&(DDtotaldf["D2_len"] >= 12)]
	DDfusion = pd.crosstab(DDtotaluse["D2_name"], DDtotaluse["D1_name"])
	emptydf = pd.DataFrame(np.zeros((len(Dlist), len(Dlist))), index = Dlist, columns = Dlist)
	DDfusion = DDfusion + emptydf
	DDfusion = DDfusion.loc[Dlist[::-1], Dlist]
	count1 = DDtotaluse.groupby(["D1_name"]).size().to_frame().reset_index()
	count1.columns = ["D1_name", "Size"]
	Dgap = pd.DataFrame({"D1_name" : Dlist, "Rank":np.arange(1,28)})
	DfillD1 = pd.merge(Dgap, count1, on = "D1_name", how = "left")
	count2 = DDtotaluse.groupby(["D2_name"]).size().to_frame().reset_index()
	count2.columns = ['D2_name', 'Size']
	D2gap = pd.DataFrame({"D2_name":Dlist, 'Rank':np.arange(1,28)[::-1]})
	DfillD2 = pd.merge(D2gap, count2, on = "D2_name", how = "left")
	plt.figure(figsize = (10, 10))
	gs = gridspec.GridSpec(2,2,height_ratios=[1,5], width_ratios=[5,1])
	ax = plt.subplot(gs[0])
	sns.barplot(x = "Rank", y = "Size", data = DfillD1, color = 'black')
	plt.xticks([])
	ax.set_xlim(-0.5,26.5)
	ax2 = plt.subplot(gs[2])
	sns.heatmap(DDfusion, cmap = 'Blues', vmin = 0, center  = 6000, vmax = 20000, cbar = False, ax = ax2, linewidths = 1, linecolor = 'white')
	plt.subplots_adjust(hspace = 0.01, wspace = 0.00)
	ax3 = plt.subplot(gs[3])
	sns.barplot(y = 'Rank', x = 'Size', data = DfillD2, color = 'black', orient = 'h')
	plt.yticks([])
	plt.ylabel('')
	plt.savefig(output)
	return DDtotaluse

def DDspan(df, output):
	DDspan = df[["D1_name", "D2_name"]]
	DDspan['D1_pos'] = DDspan['D1_name'].str.split("-",expand = True)[1].astype(np.int16)
	DDspan['D2_pos'] = DDspan['D2_name'].str.split('-', expand =  True)[1].astype(np.int16)
	DDspan['Rank'] = DDspan['D1_pos'] - DDspan["D2_pos"]
	Spanline = DDspan.groupby("Rank").size().to_frame().reset_index()
	Spanline.columns = ['Rank', "Size"]
	Spanline["Fre"] = Spanline["Size"]/Spanline["Size"].sum()
	plt.figure(figsize = (8,2))
	sns.lineplot(x = "Rank", y = "Fre", data = Spanline)
	large_span = Spanline[Spanline["Fre"] == np.max(Spanline["Fre"])]
	plt.scatter([-7.08], max(Spanline["Fre"]), c = 'red')
	plt.vlines(x = -7,ymin = -.01, ymax = 0.15, linestyle = 'dashed', alpha = 0.3, color = 'red')
	plt.ylim(-.01,0.18)
	plt.savefig(output)


def main(infile, Dlocal, output1, output2):
	DDusedf = DDHeatmap(infile, Dlocal, output1)
	#print(DDusedf)
	DDspan(DDusedf, output2)

if __name__ == "__main__":
	infile, Dlocal, output1, output2 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
	main(infile, Dlocal, output1, output2)
