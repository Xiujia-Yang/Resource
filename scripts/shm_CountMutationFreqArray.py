import re,csv,glob,sys,os
from Bio import SeqIO
import argparse

def get_ref_region(reference, infiles):
    known_id = set()
    for infile in infiles:
        ref_id = infile.split("/")[-1].split(".pos")[0].strip("z.")
        known_id.add(ref_id)

    ref_dict = {}
    fr1, fr2, fr3, cdr1, cdr2 = [],[],[],[],[]
    for rec in SeqIO.parse(reference, 'fasta'):
        rec_id = rec.description.split("\t")[0]
        region_anno_list = rec.description.split("\t")[-1].split("|")
        region_anno = [int(i) for i in region_anno_list]

        if rec_id in known_id:
            fr1.append(int(region_anno[1])-int(region_anno[0]))
            cdr1.append(int(region_anno[2])-int(region_anno[1]))
            fr2.append(int(region_anno[3])-int(region_anno[2]))
            cdr2.append(int(region_anno[4])-int(region_anno[3]))
            fr3.append(int(region_anno[5])-int(region_anno[4]))
        ref_dict[rec_id] = region_anno

    region_max = [max(fr1), max(cdr1), max(fr2), max(cdr2), max(fr3)]
    #print region_max
    return ref_dict, region_max

def get_array_pos_mutation(infiles, reference):
    ref_dict, region_max = get_ref_region(reference, infiles)
    anno_region_list = ["FR1"] * (int(region_max[0])) + ["CDR1"] * (int(region_max[1])) + \
                  ["FR2"] * (int(region_max[2])) + ["CDR2"] * (int(region_max[3])) + \
                  ["FR3"] * (int(region_max[4]))

    out = csv.writer(open("Fig.5c_d.profile.mut.rate.txt", 'w'), delimiter="\t")
    out.writerow(["Germline_id","Family","Clone_number"]+anno_region_list)
    #out2=csv.writer(open("zz.profile/zz.clone.F.num.txt", 'w'), delimiter="\t")

    for infile in infiles:
        data = []
        ref_id = infile.split("/")[-1].split(".pos")[0].strip("z.")
        if ref_id not in ref_dict.keys():
            continue
        else:
            region_anno = ref_dict[ref_id]
            data.append(ref_id)
            data.append(ref_id.split("-")[0])
            region_list=[]
            for rec in csv.reader(open(infile, 'r'), delimiter="\t"):
                if rec[0].startswith("#IGHV"):
                    reads_num =  int(rec[1].split(":")[-1])
                    #out2.writerow([ref_id, reads_num])
                elif rec[0].startswith("#Mut_type"):
                    continue
                else:
                    region_list.append(int(rec[-6])*float(100)/reads_num)

        fr1_len = region_anno[1]-region_anno[0]
        cdr1_len = region_anno[2]-region_anno[1]
        fr2_len = region_anno[3]-region_anno[2]
        cdr2_len = region_anno[4]-region_anno[3]
        fr3_len = region_anno[5]-region_anno[4]

        data = data+[reads_num]+region_list[:region_anno[1]]+(region_max[0]-fr1_len)*[-100]+ \
            region_list[region_anno[1]:region_anno[2]]+(region_max[1]-cdr1_len)*[-100]+ \
            region_list[region_anno[2]:region_anno[3]]+(region_max[2]-fr2_len)*[-100]+ \
            region_list[region_anno[3]:region_anno[4]]+(region_max[3]-cdr2_len)*[-100]+ \
            region_list[region_anno[4]:region_anno[5]]+(region_max[4]-fr3_len)*[-100]

        out.writerow(data)

    return 0

def main():
    infiles = sorted(glob.glob("IGHV*.pos.mut.txt"))
    get_array_pos_mutation(infiles, refile)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Count the position mutation")
    parser.add_argument('-r', '--reference', help='Input the used germline sequence of V segment.')
    args = parser.parse_args()
    refile = args.reference
    main()
