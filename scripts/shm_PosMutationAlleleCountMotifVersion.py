import re,csv,glob,sys,os
from common_info import *
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter
from scipy.stats import poisson
import seaborn as sns

def get_MutType(ref_nt):
    ref_codon_list = [str(ref_nt)[i:i + 3] for i in range(0, len(str(ref_nt)), 3)]
    pos_dict = {}
    mut_dict = {}
    if len(ref_codon_list[-1]) != 3:
        ref_nt = ref_nt[:-3]
    else:
        ref_nt = ref_nt

    for i,j in enumerate(ref_nt):
        ref_codon = ref_codon_list[int(i/3)]
        mydict = {}
        mod_num = (i+1)%3
        for n in ['A', 'C', 'G', 'T']:
            if n == j:
                continue
            else:
                if mod_num == 0:
                    new_codon = ref_codon[:2]+n
                elif mod_num == 1:
                    new_codon = n+ref_codon[1:]
                else:
                    new_codon = ref_codon[0]+n+ref_codon[-1]

                #print ref_codon, new_codon
                ref_aa = dict_codon2aa[ref_codon][0]
                new_aa = dict_codon2aa[new_codon][0]

                if ref_aa == new_aa:
                    mut_feature = "Syn"
                else:
                    mut_feature = "Nonsyn"
                mydict[n] = mut_feature

        mut_dict[i] = mydict

        if j == "C":
            if i > 1:
                mystr_1 = ref_nt[i-2:i+2]
                mystr_2 = ref_nt[i-2:i+1]
                match_1 = re.search('[AT][AG]C[CT]', mystr_1)
                match_2 = re.search('[CG][CT]C', mystr_2)
                if match_1:
                    mut_type = "WRCY"
                elif match_2:
                    mut_type = "SYC"
                else:
                    mut_type = "C"
            else:
                mut_type = "C"

        elif j == "G":
            if i > 0:
                mystr_1 = ref_nt[i-1:i+3]
                mystr_2 = ref_nt[i:i+3]
                match_1 = re.search('[AG]G[CT][AT]', mystr_1)
                match_2 = re.search('G[AG][CG]', mystr_2)

                if match_1:
                    mut_type = "RGYW"
                elif match_2:
                    mut_type = "GRS"
                else:
                    mut_type = "G"

            else:
                mystr_2 = ref_nt[i:i + 3]
                match_2 = re.search('G[AG][CG]', mystr_2)

                if match_2:
                    mut_type = "GRS"
                else:
                    mut_type = "G"

        elif j == "A":
            if i > 0:
                mystr_1 = ref_nt[i-1:i+1]
                match_1 = re.search('[AT]A', mystr_1)

                if match_1:
                    mut_type = "WA"
                else:
                    mut_type = "A"
            else:
                mut_type = "A"

        else:
            mystr_1 = ref_nt[i:i+2]
            match_1 = re.search('T[AT]', mystr_1)

            if match_1:
                mut_type = "TW"
            else:
                mut_type = "T"

        pos_dict[i] = mut_type

    #mut_list = [i[1] for i in sorted(pos_dict.items(), key=lambda d:int(d[0]))]

    return pos_dict, mut_dict




def get_ref_dict(refile):
    ref_dict = {}
    ref_seq_dict = {}
    region_dict = {}
    for rec in SeqIO.parse(refile, 'fasta'):
        #print rec.seq
        if "N" not in rec.seq:
            rec_id = rec.description.split("\t")[0]
            rec_e = int(rec.description.split("|")[-1])
            rec_s = int(rec.description.split("\t")[-1].split("|")[0])
            ref_seq = str(rec.seq[rec_s:rec_e])
            pos_dict, mut_dict = get_MutType(ref_seq)
            ref_dict[rec_id] = [pos_dict, mut_dict]
            ref_seq_dict[rec_id] = ref_seq
            region_anno = rec.description.split("\t")[-1].split("|")
            region_list = ["FR1"]*(int(region_anno[1])-int(region_anno[0]))+["CDR1"]*(int(region_anno[2])-int(region_anno[1]))+ \
                ["FR2"]*(int(region_anno[3])-int(region_anno[2]))+["CDR2"]*(int(region_anno[4])-int(region_anno[3]))+ \
                ["FR3"]*(int(region_anno[5])-int(region_anno[4]))
            region_dict[rec_id] = region_list

    return ref_dict, ref_seq_dict, region_dict

def add_MutType_flag(infile, ref_dict, ref_seq_dict, region_dict):
    ref_id = ".".join(infile.split(".")[:2])
    pos_ref_dict, mut_ref_dict = ref_dict[ref_id]
    syn_mut, nonsyn_mut= {}, {}
    cover_dict, multiple_dict = {}, {}
    ref_seq = ref_seq_dict[ref_id]
    ref_len = len(ref_seq)
    #print len(ref_seq)
    #print len(mut_ref_dict)

    num = 0
    for rec in csv.reader(open(infile, 'r'), delimiter="\t"):
        if rec[0] == str(0):
            continue
        else:
            num += 1
            Vinfo = rec[3].split("|")[-2]
            if Vinfo:
                loc_list = re.findall(r"\d+\.?\d*", Vinfo)
                pos_list = [int(i) for i in loc_list if int(i) < ref_len]
                codon_pos = [int(j/3) for j in pos_list]
                codon_pos_counter = Counter(codon_pos)
                myseq = "".join(rec[6:13])

                for ind, ite in enumerate(pos_list):
                    #if codon_pos_counter[codon_pos[ind]] > 1:
                        #multiple_dict[ite] = multiple_dict.setdefault(ite, 0)+1
                    #else:
                        #print myseq[ite], ite
                    mut_type = mut_ref_dict[ite][myseq[ite]]
                    if myseq[ite] == "A":
                        myarr = np.array([1,0,0,0])
                    elif myseq[ite] == "C":
                        myarr = np.array([0,1,0,0])
                    elif myseq[ite] == "G":
                        myarr = np.array([0,0,1,0])
                    else:
                        myarr = np.array([0,0,0,1])

                    if mut_type == "Syn":
                        syn_mut[ite] = syn_mut.setdefault(ite, np.array([0,0,0,0]))+myarr
                    else:
                        nonsyn_mut[ite] = nonsyn_mut.setdefault(ite, np.array([0,0,0,0]))+myarr

    out = csv.writer(open("z.mut_type.%s"%(infile), 'w'), delimiter="\t")
    out.writerow(['#Mut_type', 'Region','Ref','Pos','Cover_num','R_A','R_C','R_G','R_T','S_A','S_C','S_G','S_T'])
    ref_region = region_dict[ref_id]
    for i in range(ref_len):
        data = []
        data.append(pos_ref_dict[i])
        data.append(ref_region[i])
        data.append(ref_seq[i])
        data.append(i)
        #data.append(num-multiple_dict.get(i,0))
        data.append(num)
        data = data+list(nonsyn_mut.get(i, [0,0,0,0]))
        data = data+list(syn_mut.get(i, [0,0,0,0]))
        out.writerow(data)

    #print num
    return 0


def main():
    flag_id = ".".join(infile.split(".")[:2])
    ref_dict, ref_seq_dict, region_dict = get_ref_dict(refile)
    add_MutType_flag(infile, ref_dict, ref_seq_dict, region_dict)



if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Extract the mutation information per allele from all reads")
    parser.add_argument('-i', '--input', help="Input allele alignments.txt.")
    parser.add_argument('-d', '--outdir', default=".", help='Output file directory. Default vaule is current directory.')
    parser.add_argument('-r', '--reference', help='Input the used germline sequence of V segment.')
    args = parser.parse_args()
    infile = args.input
    outdir = args.outdir
    refile = args.reference
    main()
