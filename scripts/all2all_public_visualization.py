#!/usr/bin/env python

'''
    This script is for visualizing all-to-all public clone
    profile via seaborn heatmap.
    
    Usage: python all2all_public_visualization.py matrix.tab metadata.tab

'''

import pandas as pd
import seaborn as sns
import numpy as np
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
    # import original matrix
    df_matrix = pd.read_csv(matrix_fl, index_col=0, sep=",")

    # read the metadata append the nProdClones info to the last column of the matrix
    df_meta = pd.read_csv(meta_fl, usecols=["Run","BioProject","nProdClonesProVJ"], index_col=0, sep=",")
    order = df_meta.sort_values("nProdClonesProVJ", ascending=False).index.tolist()
    df_matrix = df_matrix[order].T[order].T
    
    # make logarithmic transformation for a better visualization
    df_matrix = np.log(df_matrix+1)

    # obtain color palette
    current_palette = ["black", "gray", "silver", "rosybrown", "firebrick", "red", \
                       "darksalmon", "sienna", "sandybrown", "bisque", "tan", "moccasin", "floralwhite", \
                       "gold", "darkkhaki", "lightgoldenrodyellow", "olivedrab", "chartreuse", "palegreen", \
                       "darkgreen", "seagreen", "mediumspringgreen", "lightseagreen", "paleturquoise", \
                       "darkcyan", "darkturquoise", "deepskyblue", "aliceblue", "slategray", "royalblue", \
                       "navy", "blue", "mediumpurple", "darkorchid", "plum", "m", "mediumvioletred", "palevioletred"]

    # assign a color for each project
    lut = dict(zip(df_meta["BioProject"].unique(), current_palette))
    # map colors for each sample
    row_colors = df_meta.sort_values("nProdClonesProVJ", ascending=False)["BioProject"].map(lut)
    col_colors = row_colors

    fig = sns.clustermap(df_matrix, row_cluster=False, col_cluster=False, cmap='Reds', \
                   row_colors=row_colors, col_colors=col_colors,\
                   yticklabels=False, xticklabels=False, figsize=(20,20))
    fig.savefig("all_to_all_comparison_heatmap.png", dpi=300)
    
if __name__ == "__main__":
    matrix_fl = sys.argv[1]
    meta_fl = sys.argv[2]
    main()
