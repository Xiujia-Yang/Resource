import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import multiprocessing as mp
import csv, os, sys


def cal(arglist):
	path, gene1, gene2 = arglist[0], arglist[1], arglist[2]
	df_file = pd.read_csv(path, sep = "\t", usecols = ["produ", "Vgene", "Dgene", "Jgene"])
	srrid = path.split("/")[-1].split("_")[0]
	median = pd.crosstab([df_file["produ"], df_file[gene1]], df_file[gene2])
	profile = median.loc["productive"]
	#profinal = profile/profile.sum().sum()
	profinal = profile
	del profile
	produ = profinal.unstack().reset_index()
	produ.columns = [gene2, gene1, "Fre"]
	produ.loc[:,"sample"] = [srrid]*len(produ)
	if 'unproductive' in median.index:
		unprofile = median.loc["unproductive"]
		unprofinal = unprofile/unprofile.sum().sum()
		unprofinal = unprofile
		del unprofile
		unprodu = unprofinal.unstack().reset_index()
		unprodu.columns = [gene2, gene1, "Fre"]
		unprodu.loc[:,"sample"] = [srrid]*len(unprodu)
	else:
		unprodu = pd.DataFrame()
	return produ, unprodu

	
def draw_figure(df, gene1, gene2, output, mapping):
	df_remove1 = df[(df["Fre"] != 1)&(df["Fre"] != 0)]
	useful = df_remove1.groupby([gene2, gene1])["Fre"]
	vlocal, vlist, vid = mapping[0], mapping[1], mapping[2]
	dlocal, dlist, did = mapping[3], mapping[4], mapping[5]
	jlocal, jlist, jid = mapping[6], mapping[7], mapping[8]
	per50 = useful.apply(lambda x:np.percentile(x,50)).to_frame().reset_index()
	per50.loc[:,"position1"] = list(map(lambda x:jlocal.get(x, np.nan), per50[gene1]))
	per50.loc[:,"position2"] = list(map(lambda x:dlocal.get(x, np.nan), per50[gene2]))
	per50_final = per50.dropna(how = "any")
	df_use = df.dropna(how = "any")
	df_use.loc[:,"position1"] = list(map(lambda x:jlocal.get(x, np.nan), df_use[gene1]))
	J_usage = df_use.groupby(["sample", "position1"])["Fre"].apply(sum).reset_index()
	J_gap = pd.DataFrame({"position1":np.arange(6)})
	J_usage_fill = pd.merge(J_usage, J_gap, how = "outer").sort_values(by = "position1")
	J_usage_fill.loc[:,"Jlabel"] = 5 - J_usage_fill["position1"]
	plt.subplots(figsize = (6,18))
	gs = gridspec.GridSpec(3,2, height_ratios = [1,12,1], width_ratios = [4,1])
	ax = plt.subplot(gs[4])
	per50_final['fre'] = per50_final["Fre"]/per50_final["Fre"].sum()
	per50_final['rank'] = np.ceil(np.log10(per50_final['Fre']))
	per50_final.to_csv("Scatterplt.csv", sep = "\t")
	sns.scatterplot(x = "position2", y = "position1", size = "rank", palette = sns.light_palette('red', n_colors = len(set(per50_final['rank']))), data = per50_final, color = "red", ax = ax)
	ax.set_ylim(-0.5, 5.5)
	ax.set_yticks(np.arange(0,6))
	ax.set_xlim(-0.5, 26.5)
	plt.xticks(np.arange(0,27), np.arange(0,27))
	plt.legend(bbox_to_anchor = (2, 1))
	ax2 = plt.subplot(gs[5])
	final = J_usage_fill.groupby("Jlabel")["Fre"].apply(np.sum).to_frame().reset_index()
	final["fre"] = final["Fre"]/final["Fre"].sum()
	final['Size'] = np.repeat(3, final.shape[0])
	final.to_csv("Lineplot.csv", sep = "\t")
	ax2.plot(final['fre'], final['Jlabel'], linewidth = 1, linestyle = 'dotted', color = 'black')
	sns.scatterplot(y = 'Jlabel', x = 'fre', data = final, color = 'red', size = 'fre', size_norm = (0,1), ax = ax2)
	ax2.set_ylim(-0.5, 5.5)
	ax2.set_xlim(0, 1)
	plt.legend(bbox_to_anchor = (2,1))
	ax2.set_yticks([])
	plt.ylabel("")
	plt.subplots_adjust(hspace = 0, wspace = 0)
	plt.savefig(output, dpi = 300, bbox_inches = "tight")
	plt.close()

def get_location(Vconf, Dconf, Jconf):
	V_local, D_local, J_local = {}, {}, {}
	V_list, D_list, J_list = [], [], []
	V_id, D_id, J_id = [], [], []
	for indV, V in enumerate(csv.reader(open(Vconf, 'r'), delimiter = "\t")):
		V_local.setdefault(V[0], 104 - indV)
		V_list.append(V[0])
		V_id.append(104 - indV)
	for indD, D in enumerate(csv.reader(open(Dconf, "r"), delimiter = "\t")):
		D_local.setdefault(D[0], indD)
		D_list.append(D[0])
		D_id.append(indD)
	for indJ, J in enumerate(csv.reader(open(Jconf, 'r'), delimiter = "\t")):
		J_local.setdefault(J[0], indJ)
		J_list.append(J[0])
		J_id.append(indJ)
	return [V_local, V_list, V_id, D_local, D_list, D_id, J_local, J_list, J_id]

def main(confFile, gene1, gene2, Vpath, Dpath, Jpath, output, indir):
	samplelist = [i[0] for i in csv.reader(open(confFile, 'r'), delimiter = "\t")]
	pos_info = get_location(Vpath, Dpath, Jpath)
	pool = mp.Pool()
	srrlist = [("%s/%s_MIXCR.txt"%(indir, i), gene1, gene2) for i in samplelist]
	result = pool.map(cal, srrlist)
	unproduList = [i[1] for i in result]
	unpro_df = pd.concat(unproduList, axis = 0).sort_values(by = gene2)
	draw_figure(unpro_df, gene1, gene2, output, pos_info)

if __name__ == "__main__":
	file,gene1, gene2 = sys.argv[1], sys.argv[2], sys.argv[3]
	Vpath, Dpath, Jpath = sys.argv[4], sys.argv[5], sys.argv[6]
	output, indir = sys.argv[7], sys.argv[8]
	main(file, gene1, gene2, Vpath, Dpath, Jpath, output, indir)
