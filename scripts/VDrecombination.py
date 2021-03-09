import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import multiprocessing as mp
import csv, os, sys


	
def draw_figure(df, gene1, gene2, output, mapping, coreV):
	df_remove1 = df[df["Number"] != 0]
	core_V = set([i[0] for i in csv.reader(open(coreV, 'r'), delimiter = "\t")])
	vlocal, vlist, vid = mapping[0], mapping[1], mapping[2]
	dlocal, dlist, did = mapping[3], mapping[4], mapping[5]
	jlocal, jlist, jid = mapping[6], mapping[7], mapping[8]
	core_dict = {}
	for V in vlocal:
		if V in core_V:
			core_dict.setdefault(vlocal[V], "Core")
	df_remove1.loc[:,"Vposition"] = list(map(lambda x:vlocal.get(x, np.nan), df_remove1["Vgene"]))
	df_remove1.loc[:,"Dposition"] = list(map(lambda x:dlocal.get(x, np.nan), df_remove1["Dgene"]))
	df_remove1.dropna(how = "any")
	df_remove1.loc[:,"Fre"] = df_remove1["Number"]/df_remove1["Number"].sum()
	df_remove1.loc[:,"rank"] = np.ceil(np.log10(df_remove1["Number"]))
	V_sum = df_remove1.groupby(["Vposition"])["Number"].apply(sum).reset_index()
	D_sum = df_remove1.groupby(["Dposition"])["Number"].apply(sum).reset_index()
	V_sum.loc[:,"Core"] = list(map(lambda x:core_dict.get(x, "Nocore"), V_sum["Vposition"]))
	D_sum_sort = D_sum.sort_values(by = "Dposition")["Number"].values
	V_gap = pd.DataFrame({"Vposition":np.arange(1,105)})
	V_usage_fill = pd.merge(V_sum, V_gap, how = "outer").sort_values(by = "Vposition")
	V_sum.loc[:,"Vlabel"] = len(vid) - V_sum["Vposition"] + 1
	V_sum_sort = V_usage_fill.sort_values(by = "Vposition")["Number"].fillna(0).values
	D_gap = pd.DataFrame({"Dposition":np.arange(0,27)})
	D_usage_fill = pd.merge(D_sum, D_gap, how = "outer").sort_values(by = "Dposition")
	plt.subplots(figsize = (6,18))
	gs = gridspec.GridSpec(3,2, height_ratios=[1,12,1],width_ratios=[4,1])
	ax = plt.subplot(gs[0])
	sns.scatterplot(x="Dposition", y = "Number", data = D_sum, ax = ax, color = "black")
	ax.plot(np.sort(D_sum["Dposition"].values), D_sum_sort, color = "black", linestyle = "dotted", linewidth = 1)
	plt.xticks(rotation = 90)
	ax.set_xticks([])
	ax.set_xlim(-0.5, 26.5)
	ax2 = plt.subplot(gs[2])
	sns.scatterplot(x = "Dposition", y = "Vposition", size = "Fre", hue = "rank",data = df_remove1, ax = ax2, palette=sns.light_palette("blue",n_colors=len(set(df_remove1["rank"]))))
	plt.legend(bbox_to_anchor = (2,1))
	ax2.set_ylim(0.5,104.5)
	ax2.set_xlim(-0.5,26.5)
	plt.yticks(np.arange(0,105), np.arange(0,105))
	plt.xticks(np.arange(0,27), np.arange(0,27))
	ax3 = plt.subplot(gs[3])
	sns.scatterplot(y = "Vposition", x = "Number", data = V_usage_fill, hue = "Core", palette = sns.color_palette(["#9b59b6","#2ecc71"]),ax=ax3)
	plt.plot(V_sum_sort, np.sort(V_usage_fill["Vposition"].fillna(0).values),color = "black", linestyle = "dotted", linewidth = 1)
	ax3.xaxis.set_ticks_position("top")
	plt.xticks([0,1000000], [" ", " "])
	plt.ylim(0,104.5)
	ax3.set_yticks([])
	plt.legend(bbox_to_anchor = (3, 1))
	plt.ylabel("")
	plt.subplots_adjust(hspace = 0, wspace = 0)
	plt.savefig(output, dpi = 300, bbox_inches = "tight")
	plt.close()


def get_location(Vconf, Dconf, Jconf):
	V_local, D_local, J_local = {}, {}, {}
	V_list, D_list, J_list = [], [], []
	V_id, D_id, J_id = [], [], []
	all_V_id = [i for i in csv.reader(open(Vconf, 'r'), delimiter = "\t")]
	for indV, V in enumerate(all_V_id):
		V_local.setdefault(V[0], len(all_V_id) - indV)
		V_list.append(V[0])
		V_id.append(len(all_V_id) - indV)
	for indD, D in enumerate(csv.reader(open(Dconf, "r"), delimiter = "\t")):
		D_local.setdefault(D[0], indD)
		D_list.append(D[0])
		D_id.append(indD)
	for indJ, J in enumerate(csv.reader(open(Jconf, 'r'), delimiter = "\t")):
		J_local.setdefault(J[0], indJ)
		J_list.append(J[0])
		J_id.append(indJ)
	return [V_local, V_list, V_id, D_local, D_list, D_id, J_local, J_list, J_id]

def main(File, gene1, gene2, Vpath, Dpath, Jpath, core_V, output):
	df = pd.read_csv(File, sep = "\t")
	pos_info = get_location(Vpath, Dpath, Jpath)
	draw_figure(df, gene1, gene2, output, pos_info, core_V)
	
			

if __name__ == "__main__":
	file,gene1, gene2 = sys.argv[1], sys.argv[2], sys.argv[3]
	V, D, J, coreV = sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]
	output = sys.argv[8]
	main(file, gene1, gene2, V, D, J, coreV, output)
