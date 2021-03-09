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


def main(infile, metafile, Clonefile, output):
	DD_new = pd.read_csv(infile, sep = "\t")
	meta = pd.read_csv(metafile)
	DD_new['D1_len'] = DD_new["D1_seq"].apply(len)
	DD_new['D2_len'] = DD_new['D2_seq'].apply(len)
	DDuse = DD_new[(DD_new['D1_len'] >= 12) & (DD_new["D2_len"] >= 12)]
	final = pd.merge(DDuse, meta, left_on = 'Sample', right_on = 'Run', how = 'left')
	Sizedf = DDuse.groupby(["Sample", 'Isotype']).size().to_frame().reset_index()
	Sizedf.columns = ['Run', 'variable', 'Size']
	Clonedf = pd.read_csv(Clonefile)
	Clonedf.rename(columns = {"Unnamed: 0":"Run"}, inplace = True)
	Clonedf = Clonedf.melt(id_vars = ['Run'])
	Cloneuse = Clonedf[Clonedf['value'] >= 5000]
	Cloneuse['Useid'] = Cloneuse['Run'] + '|' + Cloneuse['variable']
	result = pd.merge(Sizedf, Clonedf, on = ['Run', 'variable'], how = 'left')
	result['Fre'] = result['Size']/result['value']
	result = result[(result['value'] >= 5000)]
	figurefinal = pd.merge(result, meta[['Run', 'CellType']], on = 'Run')
	figurefinal['Legend'] = figurefinal['CellType'] + '|' + figurefinal['variable']
	figurefinal = figurefinal[(figurefinal['CellType'] != 'Unknown') & (figurefinal["CellType"] != 'Other')]
	uselist = ['Memory|IGHA',  'Memory|IGHG', 'Memory|IGHM', 'Naive|IGHD', "Naive|IGHM", 'Plasma|IGHA',
           'Plasma|IGHG','Plasma|IGHM', 'Unsorted|IGHA', 'Naive|IGHD', 
           'Unsorted|IGHG', 'Unsorted|IGHM']
	plt.figure(figsize = (5,2))
	sns.boxplot(x = 'Legend', y = 'Fre', data = figurefinal,showfliers = False, order = uselist, color = 'cornflowerblue', width = 0.8,
           boxprops = {"linewidth":1})
	plt.xticks(rotation = 90)
	plt.ylim(0, 0.0045)
	plt.yticks(np.arange(0,0.004, 0.001))
	plt.savefig(output, bbox_inches = 'tight')

if __name__ == "__main__":
	infile, metafile, Clonefile, output = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
	main(infile, metafile, Clonefile, output,)
