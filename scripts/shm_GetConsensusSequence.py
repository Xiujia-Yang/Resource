import re,csv,glob,sys,os
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from multiprocessing import Pool
from collections import Counter
import Levenshtein
import time

def get_consensus_seq(infile):
    mylist = []
    mydict = {}
    num = 0
    for rec in csv.reader(open(infile, 'rU'), delimiter = "\t"):
        num += 1
        if num == 1:
            continue
        else:
            myseq = "".join(rec[6:11])
            mylist.append(Seq(myseq))
            mydict.setdefault(myseq, []).append(rec)

    m = motifs.create(mylist)
    seq2 = str(m.consensus)

    myarr = []
    for i in mylist:
        seq1 = str(i)
        dis = Levenshtein.hamming(seq1, seq2)
        myarr.append(dis)

    pos = myarr.index(min(myarr))
    seq = str(mylist[pos])
    mydata = mydict[seq][0]

    return mydata


def main():
    start = time.clock()
    work_dir = os.getcwd()
    os.chdir(work_dir)
    sample_id = work_dir.split("/")[-2]
    flag_id = work_dir.split("/")[-1]

    infiles=glob.glob("clone_*.txt")
    out1 = csv.writer(open("%s-%s-consensus.txt"%(sample_id, flag_id), 'w'), delimiter = "\t")

    pool=Pool()
    for infile in infiles:
        res = pool.apply_async(get_consensus_seq, args=(infile,))
        out1.writerow(res.get())

    pool.close()
    pool.join()

    '''
    for infile in infiles:
        #print infile
        mydata = get_consensus_seq(infile)
        out1.writerow(mydata)
    '''

    end = time.clock()
    print (end - start)



if __name__ == "__main__":
    main()
