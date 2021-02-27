#!/usr/bin/env python

import csv, sys, os, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def main():
	df = pd.read_csv(infile, usecols=['Gene_order', 'Clone_fraction', 'Run'], sep='\t')
	df_group_gene = df.groupby(by=['Gene_order']).median()
	gene_orders, clone_fractions = df_group_gene.index.values, df_group_gene['Clone_fraction'].values
	slope_x, slope_y = [], []
	min_flag = False
	for index, order in enumerate(gene_orders):
		if index >= 5 and index < len(gene_orders) - 5:
			slope = (clone_fractions[index+5] - clone_fractions[index-5])/1100
			slope_x.append(index)
			slope_y.append(slope)
			if slope < 0.001 and min_flag==False:
				marked_x = index
				marked_y = slope
				min_flag = True
				out.writerow([index, slope])

	fig = plt.figure(figsize=(4, 3))
	ax1 = fig.add_subplot(111)
	sns.lineplot(x=df_group_gene.index.values, y=df_group_gene['Clone_fraction'], color='r', linewidth=1, ax=ax1)
	ax2 = plt.twinx()
	sns.lineplot(x=slope_x, y=slope_y, ax=ax2, linewidth=1)
	ax2.text(x=marked_x, y=marked_y, s='*')
	for ax in [ax1, ax2]:
		for line in ax.yaxis.get_ticklines():
			line.set_markeredgewidth(0.5)
	for ax in [ax1, ax2]:
		for line in ax.xaxis.get_ticklines():
			line.set_markeredgewidth(0.5)
	for ax in [ax1, ax2]:
		for axis in ['top', 'bottom', 'left', 'right']:
			ax.spines[axis].set_linewidth(0.5)
	plt.savefig('%s/%s_median_slope.png'%(outDir, os.path.basename(infile).split('.')[0]), dpi=800)


if __name__ == '__main__':
	infile = sys.argv[1]
	outDir = sys.argv[2]
	out = csv.writer(open('%s/%s_median_slope.txt'%(outDir, os.path.basename(infile).split('.')[0]), 'w'),delimiter='\t')
	main()

