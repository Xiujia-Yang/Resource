#!/zzh_gpfs/apps/python/bin/python
'''
    This script is for combining the statistics.
    
    Usage: python motif_mut_freq_merge_for_multiple_samples.py pathfile
'''
import pandas as pd
import glob, sys

def main():
    for i, path in enumerate(f):
        path = path.strip()
        df_temp = pd.read_csv(path, index_col=0, sep='\t')
        if i > 0:
            df = df + df_temp
        else:
            df = df_temp

    # reorder the row
    df = df.loc[['WRCY', 'RGYW', 'SYC', 'GRS', 'C', 'G', 'WA', 'TW', 'A', 'T'],]
    df.index.name = "Motif"
    df.rename(columns={"total":"Germline_num", "mut":"Mutated_num", "A":"A_mut", "C":"C_mut", "G":"G_mut", "T":"T_mut"}, inplace=True)
    df["Percentage"] = df.apply(lambda x:x.Mutated_num/x.Germline_num*100 if x.Germline_num != 0 else 0, axis=1)
    df_temp = df.apply(lambda x:x[["A_mut", "C_mut", "G_mut", "T_mut"]]/x[["A_mut", "C_mut", "G_mut", "T_mut"]].sum()*100, axis=1).fillna(0)
    df_temp.rename(columns={"A_mut":"A", "C_mut":"C", "G_mut":"G", "T_mut":"T"}, inplace=True)
    columns = ["Percentage","Mutated_num","Germline_num","A","C","G","T","A_mut","C_mut","G_mut","T_mut"]
    df = pd.concat([df, df_temp], axis=1)[columns].astype(float).round(2)
    df.to_csv("motif_mut_profile_merged.txt", sep="\t")

    
if __name__ == "__main__":
    path_file = sys.argv[1]
    f = open(path_file, "r")
    main()
    
