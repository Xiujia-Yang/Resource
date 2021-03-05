from pyforest import *
from common_info import *

def CountPurifiedSynonymous(infile):
    mydf=pd.read_csv(infile, sep="\t")
    seqnt = mydf["Ref"].tolist()
    nt_list = ['A','C','G','T']
    flag_list = []

    mod_list = []
    for num in range(len(seqnt)):
        mod_list.append(num%3)
    mymod_list = mod_list[::-1]

    for i,j in enumerate(seqnt):
        #mylist = nt_list.remove(j)
        flag_set = set()
        for m in nt_list:
            if m == j:
                continue
            else:
                mod_num = mymod_list[i]
                if mod_num == 2:
                    new_codon = m+"".join(seqnt[i+1:i+3])
                    old_codon = "".join(seqnt[i:i+3])
                elif mod_num == 1:
                    if i-1 > -1:
                        new_codon = seqnt[i-1]+m+seqnt[i+1]
                        old_codon = "".join(seqnt[i-1:i+2])
                    else:
                        new_codon, old_codon = "NA", "NA"
                else:
                    if i-2 > -1:
                        new_codon = "".join(seqnt[i-2:i])+m
                        old_codon = "".join(seqnt[i-2:i+1])
                    else:
                        new_codon, old_codon = "NA", "NA"

                if new_codon != "NA":
                    #print (i,j)
                    #print (new_codon, old_codon)
                    if dict_codon2aa[new_codon] == dict_codon2aa[old_codon]:
                        flag_set.add("Syn")
                    else:
                        flag_set.add("Nonsyn")
                else:
                    flag_set.add("Valid")

        if "Valid" in flag_set:
            myflag = "Valid"
        else:
            if "Nonsyn" not in flag_set:
                myflag = "Purified_Syn"
            elif "Syn" not in flag_set:
                myflag = "Purified_Nonsyn"
            else:
                myflag = "Both"

        flag_list.append(myflag)
    mydf["Purified_flag"] = flag_list
    mydf.to_csv(infile+".flag", sep="\t", index=False)

    return 0


def main():
    infiles=glob.glob("z.mut_type.IG*.txt")
    for infile in infiles:
        CountPurifiedSynonymous(infile)





if __name__ == "__main__":
    main()