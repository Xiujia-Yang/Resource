#!/usr/bin/env python

import csv, sys, os, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import mpl_toolkits.axisartist as AA
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec

def get_LL(infile_name):
	LL = []
	for rec in csv.reader(open(infile_name, 'r'),delimiter='\t'):
		LL.append(int(rec[-1]))
	return LL

def main():
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	spacing = 0.05

	rect_heatmap = [left, bottom, width, height]
	rect_left_line = [left + width + spacing, bottom, 0.08, height]
	rect_above_line = [left, bottom + height + spacing, 0.52, 0.05]

	# start with a rectangular Figure
	plt.figure(figsize=(4, 7))

	ax_heatmap = plt.axes(rect_heatmap)
	ax_left_line = plt.axes(rect_left_line)
	ax_above_line = plt.axes(rect_above_line)
	# heatmap
	matrix = pd.read_csv(matrixFile, sep='\t', index_col=0)
	all_genes = matrix.index.values
	sns.heatmap(matrix, cmap='Reds', ax=ax_heatmap, xticklabels=False)
	for line in ax_heatmap.yaxis.get_ticklines():
		line.set_markeredgewidth(0)
	ax_heatmap.set_yticks(np.arange(len(all_genes))+1)
	ax_heatmap.set_yticklabels(all_genes, fontsize=4)

	
	# right line plot
	sample_num_LL = get_LL(sampNumFile)
	ax_left_line.barh(np.arange(len(sample_num_LL)), sample_num_LL[::-1], color='black', edgecolor='black')
	ax_left_line.set_ylim(0,len(all_genes))
	ax_left_line.set_yticklabels([])
	ax_left_line.set_xticks(np.arange(0, 1805, 200))
	ax_left_line.set_xticklabels(np.arange(0,1805,200), fontsize=5, rotation=90)
	for tick in ax_left_line.get_xticklabels():
	    tick.set_rotation(90)

	# upper line plot
	gene_num_LL = get_LL(geneNumFile)
	ax_above_line.bar(np.arange(len(gene_num_LL)), gene_num_LL, color='black', edgecolor='black')
	ax_above_line.set_xticklabels([])
	ax_above_line.set_xlim(0,2152)

	plt.savefig('%s/%s_heatmap_bar.png'%(outDir, os.path.basename(matrixFile).split('.')[0]),dpi=800)



if __name__ == '__main__':
	matrixFile = sys.argv[1]
	sampNumFile = sys.argv[2]
	geneNumFile = sys.argv[3]
	outDir = sys.argv[4]
	main()

