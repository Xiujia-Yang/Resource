#!/usr/bin/env python

import csv, sys, os, glob
import numpy as np
import pandas as pd

def main():
	dfUsage = pd.read_csv(usageFile, sep='\t', index_col=0)
	# genes are sorted by the number of samples
	dfUsage['sampNum'] = dfUsage.count(axis=1)
	dfUsage = dfUsage.sort_values(['sampNum'], ascending=False)
	dfUsage[['sampNum']].to_csv('%s/%s_gene_sample_num.txt'%(outDir, geneType), sep='\t', header=False)
	dfUsage.drop(['sampNum'], axis=1, inplace=True)
	# samples are sorted by the number of reads
	dfUsage = dfUsage.transpose()
	dfUsage['Run'] = dfUsage.index.values
	dfUsage = dfUsage.merge(dfMeta, on=['Run'], how='inner')
	dfUsage = dfUsage.sort_values(['nProdReads'], ascending=True)
	dfUsage[['Run', 'nProdReads']].to_csv('%s/%s_gene_reads_num.txt'%(outDir, geneType), sep='\t', index=False, header=False)
	dfUsage.index = dfUsage['Run'].values
	dfUsage.drop(['Run', 'nProdReads'], axis=1, inplace=True)
	dfUsage['geneNum'] = dfUsage.count(axis=1)
	dfUsage[['geneNum']].to_csv('%s/%s_gene_num_per_sample.txt'%(outDir, geneType), sep='\t', header=False)
	dfUsage.drop(['geneNum'], axis=1, inplace=True)
	# calculate threshold which is equal to
	# the mean of the maximum for all samples
	thred = np.mean(dfUsage.max(axis=1))
	# output matrix for heatmap
	dfUsage = dfUsage.transpose()
	dfUsage.fillna(0, inplace=True)
	dfUsage[dfUsage > thred] = thred
	dfUsage.to_csv('%s/%s_gene_usage.txt'%(outDir, geneType), sep='\t')

if __name__ == '__main__':
	usageFile = sys.argv[1]
	geneType = sys.argv[2]
	outDir = sys.argv[3]
	main()
