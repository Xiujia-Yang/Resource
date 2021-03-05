import re,csv,glob,sys,os
from Bio import SeqIO
from common_info import *
import argparse

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

def unique_seq_remove_nomut_reads(reader_handle):
    uniq_dict = {}
    for i in reader_handle:
        myseq = "".join(i[6:13])
        if myseq in set(uniq_dict.keys()):
            continue
        else:
            uniq_dict[myseq] = i

    mylist = []
    for n in uniq_dict.values():
        mylist.append(n)

    return mylist

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
            ref_seq_dict[rec_id] = rec.seq
            region_anno = rec.description.split("\t")[-1].split("|")
            region_list = ["FR1"]*(int(region_anno[1])-int(region_anno[0]))+["CDR1"]*(int(region_anno[2])-int(region_anno[1]))+ \
                ["FR2"]*(int(region_anno[3])-int(region_anno[2]))+["CDR2"]*(int(region_anno[4])-int(region_anno[3]))+ \
                ["FR3"]*(int(region_anno[5])-int(region_anno[4]))
            region_dict[rec_id] = region_list

    return ref_dict, ref_seq_dict, region_dict


def add_MutType_flag(infile, ref_dict, ref_seq_dict, region_dict):
    ref_id = infile.split(".txt")[0]
    pos_ref_dict, mut_ref_dict = ref_dict[ref_id]
    reader_handle = csv.reader(open(infile, 'rU'), delimiter = "\t")
    pos_mut, mut_dict, synmut_dict = {}, {}, {}

    num = 0
    #record_list = unique_seq_remove_nomut_reads(reader_handle)
    #for rec in record_list:
    for rec in reader_handle:
        Vinfo = rec[3].split("|")[-2]
        num = num +1
        if Vinfo:
            rec_seq = "".join(rec[6:11])
            mut_list = re.findall(r"S[ATGC][0-9]{0,3}[ATGC]", Vinfo)
            loc_list = re.findall(r"\d+\.?\d*", Vinfo)

            for i in range(len(mut_list)):
                #print loc_list[i]
                #print len(rec_seq)
                if int(loc_list[i]) < len(rec_seq): #Just V segment
                    #print loc_list[i]
                    mut_char = mut_list[i][-1]
                    mut_type = mut_ref_dict[int(loc_list[i])][mut_char]
                    pos_mut.setdefault(int(loc_list[i]), []).append(mut_char)
                    mut_dict[int(loc_list[i])] = mut_dict.setdefault(int(loc_list[i]), 0)+1

                    if mut_type == "Syn":
                        synmut_dict[int(loc_list[i])] = synmut_dict.setdefault(int(loc_list[i]), 0)+1

    out = csv.writer(open("z."+infile.split(".txt")[0]+".pos.mut.txt", "w"), delimiter = "\t")
    out.writerow(["#%s"%(ref_id), "#Total number:%d"%(num)])
    out.writerow(['#Mut_type','region','ref','pos','A','T','C','G','Mut_number','Synonymous_number','Nonsynonymous_number','Mut_percent','Synonymous_percent','Nonsynonymous_percent'])

    #print pos_mut, mut_dict, synmut_dict
    ref_nt = ref_seq_dict[ref_id]
    region_list = region_dict[ref_id]
    for n in range(len(pos_ref_dict.keys())):
        data = []
        data.append(pos_ref_dict[n])
        data.append(region_list[n])
        data.append(ref_nt[n])
        data.append(n)
        if n in pos_mut.keys():
            data.append(pos_mut[n].count('A'))
            data.append(pos_mut[n].count('T'))
            data.append(pos_mut[n].count('C'))
            data.append(pos_mut[n].count('G'))
        else:
            data = data + [0]*4

        if n in mut_dict.keys():
            data.append(mut_dict[n])
            if n in synmut_dict.keys():
                data.append(synmut_dict[n])
                data.append(mut_dict[n] - synmut_dict[n])
            else:
                data.append(0)
                data.append(mut_dict[n])
        else:
            data.append(0)
            data.append(0)
            data.append(0)

        if n in mut_dict.keys():
            mut_percent = mut_dict[n]*float(100)/num
        else:
            mut_percent = 0
        data.append(mut_percent)

        if n in synmut_dict.keys():
            syn_percent = float(synmut_dict[n])*100/num
            data.append(syn_percent)
            data.append((mut_dict[n] - synmut_dict[n])*float(100)/num)
        else:
            data.append(0)
            if n in mut_dict.keys():
                data.append(mut_dict[n]*float(100)/num)
            else:
                data.append(0)

        #print data
        out.writerow(data)

    return 0


def main():
    flag_id = infile.split(".txt")[0]
    ref_dict, ref_seq_dict, region_dict = get_ref_dict(refile)
    if flag_id in ref_dict.keys():
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
