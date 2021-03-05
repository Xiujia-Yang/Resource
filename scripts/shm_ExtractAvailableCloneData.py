from __future__ import division
import csv, sys, os, glob, re
import numpy as np
import subprocess
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse

def group_clone(clone_file):
    mydict = {}
    isotype_dict = {}
    for rec in csv.reader(open(clone_file, 'r'), delimiter = "\t"):
        if rec[0].startswith("cloneId"):
            continue
        else:
            V_align = rec[5].split(",")[0].split("(")[0].strip().replace("*", ".").strip()
            isotype = rec[8].split(",")[0][:4]
            if isotype:
                #print isotype
                if "V" in V_align:
                    mydict[int(rec[0])] = V_align
                    isotype_dict.setdefault(isotype, []).append(int(rec[0]))

    return mydict, isotype_dict

def judge_productive_and_unproductive(rec):
    flag = "NA"
    if "," not in rec[3]:
        #print len(rec)
        if rec[6] and rec[7] and rec[8] and rec[9] and rec[10] and rec[11] and rec[12]:
            seq_aa = "".join(rec[13:20]).strip("_")
            Vinfo = rec[3].split("|")[-2]
            if Vinfo:
                if ("I" not in Vinfo) and ("D" not in Vinfo):
                    if ("*" not in seq_aa) and ("_" not in seq_aa):
                        flag = "productive"
                    else:
                        flag = "unproductive"
            else:
                if ("*" not in seq_aa) and ("_" not in seq_aa):
                    flag = "productive"
                else:
                    flag = "unproductive"

    return flag


def get_group_data(clone_file, align_file):
    read_abundance = {}
    productive_read_dict = {}
    clone_dict, isotype_dict = group_clone(clone_file)
    seq_set = set()
    clone_vr_dict = {}
    tmp_num = 0
    for rec in csv.reader(open(align_file, 'rU'), delimiter = "\t"):
        if rec[0].startswith("readId"):
            continue
        else:
            seq_nt = "".join(rec[6:13])
            flag = judge_productive_and_unproductive(rec)
            assign_v = rec[2].replace("*", ".").strip()
            clone_num = int(rec[1])
            if flag == "productive":
                if int(rec[1]) in clone_dict.keys() and assign_v == clone_dict[int(rec[1])]:
                    if seq_nt in seq_set and clone_num == clone_vr_dict.get(seq_nt,-2):
                        read_abundance[seq_nt] = read_abundance[seq_nt] + 1
                    elif seq_nt not in seq_set:
                        read_abundance[seq_nt] = 1
                        seq_set.add(seq_nt)
                        clone_vr_dict[seq_nt] = clone_num
                        tmp_num += 1
                        productive_read_dict.setdefault(int(rec[1]), []).append(rec)
                    else:
                        continue

    for k2,v2 in isotype_dict.items():
        #print k2,v2
        os.system("mkdir -p %s"%(k2))
        for n in v2:
            if productive_read_dict.get(n, False):
                v1 = productive_read_dict[n]
                mynum = 0
                mylist = []
                for x in v1:
                    vr = "".join(x[6:13])
                    mynum = read_abundance[vr] + mynum
                    new_record = [x[0]+"_"+str(read_abundance[vr])]+x[1:]
                    mylist.append(new_record)

                if mynum > 1:
                    mydf = pd.DataFrame(np.array(mylist))
                    mydf.to_csv("%s/clone_productive_%s.txt"%(k2,n), sep="\t", index=False)

    return 0


def main():
    work_dir = os.getcwd()
    os.chdir(work_dir)

    align_file = sys.argv[1]#"alignments.txt"
    clone_file = sys.argv[2]#"clones.txt"
    get_group_data(clone_file, align_file)


if __name__ == '__main__':
    main()
