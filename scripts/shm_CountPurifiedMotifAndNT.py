from pyforest import *
import csv

def CountMotif(infiles, mutation_flag):
    i = 0
    for infile in infiles:
        i = i+1
        mydf = pd.read_csv(infile, sep="\t")
        if mutation_flag == "All":
            new_df = mydf
        else:
            new_df = mydf[mydf["Purified_flag"] == mutation_flag][mydf.columns]
        new_df.reset_index(drop=True)
        if i == 1:
            df_all = new_df
        else:
            df_all=pd.concat([df_all, new_df], ignore_index=True)

    total_df = df_all[[df_all.columns[0], df_all.columns[2]]+df_all.columns.tolist()[4:13]]
    #print (total_df)
    total_df.reset_index(drop=True)
    total_df.to_csv("Fig5.a-b.%s.motif.txt"%(mutation_flag), sep="\t", index=False)

    return total_df


def CountNTtransition(total_df, mutation_flag):
    mydf = total_df.groupby(["Ref"])
    mydict = dict(zip(total_df.columns.tolist()[2:], ['sum'] * (len(total_df.columns) - 2)))
    myref_df = mydf.agg(mydict)
    myref_df["Mutated_num"] = myref_df.apply(lambda x:x['S_A']+x['S_C']+x['S_G']+x['S_T']+x['R_A']+x['R_C']+x['R_G']+x['R_T'], axis=1)
    myref_df["Percentage"] = myref_df.apply(lambda x:x["Mutated_num"]*100/x["Cover_num"], axis=1)
    myref_df["Percent_A"] = myref_df.apply(lambda x:(x["S_A"]+x['R_A'])*100/x["Mutated_num"], axis=1)
    myref_df["Percent_C"] = myref_df.apply(lambda x:(x["S_C"]+x['R_C'])*100/x["Mutated_num"], axis=1)
    myref_df["Percent_G"] = myref_df.apply(lambda x:(x["S_G"]+x['R_G'])*100/x["Mutated_num"], axis=1)
    myref_df["Percent_T"] = myref_df.apply(lambda x:(x["S_T"]+x['R_T'])*100/x["Mutated_num"], axis=1)
    myref_df["A_mut"] = myref_df.apply(lambda x:(x["S_A"]+x['R_A']), axis=1)
    myref_df["C_mut"] = myref_df.apply(lambda x:(x["S_C"]+x['R_C']), axis=1)
    myref_df["G_mut"] = myref_df.apply(lambda x:(x["S_G"]+x['R_G']), axis=1)
    myref_df["T_mut"] = myref_df.apply(lambda x:(x["S_T"]+x['R_T']), axis=1)

    new_df=myref_df[["Percentage", "Mutated_num","Cover_num","Percent_A","Percent_C","Percent_G","Percent_T","A_mut","C_mut","G_mut","T_mut"]]
    new_df.columns=["Percentage","Mutated_num","Germline_num","A","C","G","T","A_mut","C_mut","G_mut","T_mut"]
    new_df.to_csv("Supp.%s.ACGT.txt"%(mutation_flag), sep="\t", index=True)

    nt_dict = new_df.to_dict()
    mydata = []
    for nt in ["A", "C", "G", "T"]:
        data = []
        data.append(nt)
        data.append(nt_dict["A_mut"][nt]*float(100)/nt_dict["Germline_num"][nt])
        data.append(nt_dict["C_mut"][nt]*float(100)/nt_dict["Germline_num"][nt])
        data.append(nt_dict["G_mut"][nt]*float(100)/nt_dict["Germline_num"][nt])
        data.append(nt_dict["T_mut"][nt]*float(100)/nt_dict["Germline_num"][nt])
        mydata.append(data)

    myfinal_df = pd.DataFrame(mydata, columns=["Motif"]+["A", "C", "G", "T"])
    myfinal_df.set_index("Motif", inplace=True)
    final_df = myfinal_df.T
    final_df.to_csv("Supp.%s.NT.transition.array.txt"%(mutation_flag), sep="\t")

    fig, ax1 = plt.subplots(figsize=(5, 5))
    sns.heatmap(final_df, cmap="Blues", cbar=False, vmin=0, vmax=8, robust=False)
    plt.savefig("Supp.%s.Heatmap_nt_transform.png" % (mutation_flag), dpi=600)
    plt.close()

    fig, ax1 = plt.subplots(figsize=(5, 5))
    sns.heatmap(final_df, cmap="Blues",vmin=0, vmax=8, robust=False)
    plt.savefig("Supp.%s.Heatmap_nt_transform.cbar.png" % (mutation_flag), dpi=600)
    plt.close()

    return 0


def CountMotifCalculation_motif(total_df, mutation_flag):
    mydf = total_df.groupby(["#Mut_type"])
    mydict = dict(zip(total_df.columns.tolist()[2:], ['sum'] * (len(total_df.columns) - 2)))
    myref_df = mydf.agg(mydict)
    myref_df["Mutated_num"] = myref_df.apply(lambda x:x['S_A']+x['S_C']+x['S_G']+x['S_T']+x['R_A']+x['R_C']+x['R_G']+x['R_T'], axis=1)
    myref_df["Percentage"] = myref_df.apply(lambda x:x["Mutated_num"]*100/x["Cover_num"], axis=1)
    myref_df["Percent_A"] = myref_df.apply(lambda x:(x["S_A"]+x['R_A'])*100/x["Mutated_num"], axis=1)
    myref_df["Percent_C"] = myref_df.apply(lambda x:(x["S_C"]+x['R_C'])*100/x["Mutated_num"], axis=1)
    myref_df["Percent_G"] = myref_df.apply(lambda x:(x["S_G"]+x['R_G'])*100/x["Mutated_num"], axis=1)
    myref_df["Percent_T"] = myref_df.apply(lambda x:(x["S_T"]+x['R_T'])*100/x["Mutated_num"], axis=1)
    myref_df["A_mut"] = myref_df.apply(lambda x:(x["S_A"]+x['R_A']), axis=1)
    myref_df["C_mut"] = myref_df.apply(lambda x:(x["S_C"]+x['R_C']), axis=1)
    myref_df["G_mut"] = myref_df.apply(lambda x:(x["S_G"]+x['R_G']), axis=1)
    myref_df["T_mut"] = myref_df.apply(lambda x:(x["S_T"]+x['R_T']), axis=1)

    new_df=myref_df[["Percentage", "Mutated_num","Cover_num","Percent_A","Percent_C","Percent_G","Percent_T","A_mut","C_mut","G_mut","T_mut"]]
    new_df.columns=["Percentage","Mutated_num","Germline_num","A","C","G","T","A_mut","C_mut","G_mut","T_mut"]
    motif_dict = new_df.to_dict()
    #print (motif_dict)
    motif_list = ['WRCY', 'RGYW', 'SYC', 'GRS', 'C', 'G', 'WA', 'TW', 'A', 'T']
    out = csv.writer(open("Fig6.a_b.%s.motif.txt"%(mutation_flag), 'w'), delimiter="\t")
    out.writerow(["#Mut_type","Percentage","Mutated_num","Germline_num","A","C","G","T","A_mut","C_mut","G_mut","T_mut"])
    for motif in motif_list:
        data = []
        data.append(motif)
        for k in motif_dict.keys():
            data.append(motif_dict[k].get(motif, 0))
        out.writerow(data)

    return 0

def main():
    infiles = glob.glob("z.mut_type.IGHV*.txt.flag")
    for mutation_flag in ["Purified_Syn", "Both", "Purified_Nonsyn", "All"]:
        total_df = CountMotif(infiles, mutation_flag)
        CountMotifCalculation_motif(total_df, mutation_flag)
        CountNTtransition(total_df, mutation_flag)



if __name__ == "__main__":
    main()
