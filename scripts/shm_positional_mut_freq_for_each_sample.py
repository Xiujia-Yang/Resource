#!/usr/bin/env python

'''
    This script is for perform basic statistic analysis
    of somatic hypermutation for each sample. The invol
    -lved algorthim is based on weighted-matrix, which 
    considered all available reads within a clone.

    This is second version, in which only clone with 
    nProdReads great than 1 is considered.

    Usage: python positional_mut_freq_for_each_sample.py sample alignments.txt clones.txt isotype outdir
'''

import csv, sys, re, os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd

def mutation_tag_processing(allele, aligntagList, alleleDict, alleleInfoDf):
    length, seq = alleleInfoDf[alleleInfoDf["vAllele"]==allele].values.tolist()[0][6:]
    if allele in alleleDict.keys():
        alleleDict[allele]["nReads"] = map(lambda x:x+1, alleleDict[allele]["nReads"])
        readsNum = float(len(aligntagList))
        for aligntag in aligntagList:
            # extract mutation position and target nucleotide 
            # and store them into a tuple
            pos_tar_tup = [(int(x[2:-1]), x[-1]) for x in re.findall(r'S[ATCG]\d+[ATCG]', aligntag)]
            for pos, target in pos_tar_tup:
                ## deal with the situation in which pos exceeds the length
                ## for we consider only FR1-FR3
                if pos < length:
                    alleleDict[allele][target][pos] += 1.0/readsNum

    else:
        # initialize the allele information for later process
        alleleDict.update({allele: {"pos":range(length), "ref":[x for x in seq[:length]], \
                        "A":[0.0]*length, "C":[0.0]*length, "T":[0.0]*length, "G":[0.0]*length, \
                         "is_synA":["False"]*length, "is_synC":["False"]*length, \
                           "is_synT":["False"]*length, "is_synG":["False"]*length, \
                          "nReads":[0]*length}})
        offset = length%3
        for tar in ["A", "T", "C", "G"]:
            for pos in range(length):
                if pos < offset:
                    alleleDict[allele]["is_syn"+tar][pos] = "None"
                else:
                    if (pos-offset)%3 == 0:
                        codon = seq[pos:pos+3]
                        new_codon = tar + seq[pos+1:pos+3]
                    elif (pos-offset)%3 == 1:
                        codon = seq[pos-1:pos+2]
                        new_codon = seq[pos-1] + tar + seq[pos+1]
                    else:
                        codon = seq[pos-2:pos+1]
                        new_codon = seq[pos-2:pos] + tar
                    aa = Seq(codon, generic_dna).translate()
                    new_aa = Seq(new_codon, generic_dna).translate()
                    alleleDict[allele]["is_syn"+tar][pos] = "Yes" if aa == new_aa else "No"
        mutation_tag_processing(allele, aligntagList, alleleDict, alleleInfoDf)

def is_productive(s):
    if s[-1] == "_":
        s = s[:-1]
    if s.count('*') == 0 and s.count('_') == 0:
        return "Productive"
    else:
        return "Unproductive"
        
def main():
    # Create a directory for each sample to contain position-specific 
    # mutation information for each found allele if it does not exist before
    if not os.path.exists("%s/productive/%s-%s"%(outdir, sample, iso)):
        os.makedirs("%s/productive/%s-%s"%(outdir, sample, iso))
    if not os.path.exists("%s/unproductive/%s-%s"%(outdir, sample, iso)):
        os.makedirs("%s/unproductive/%s-%s"%(outdir, sample, iso))
    
    # define a dictionary to contain positional mutation information for
    # each allele. Key: allele. Value (dict): {pos:[], ref:[], A:[], C:[], \
    # T:[], G:[], is_synA:[], is_synC:[], is_synT:[], is_synG:[]}
    # Extracted aligntag were be processed and statistics will be stored
    # in this dict. When all reads were processed, we can tranform each
    # allele to a DataFrame object for easy output.
    alleleDictProd = {}
    alleleDictUnprod = {}
    # read the allele region delimitation information
    df_delim = pd.read_csv("data/v.allele.region.delimitation.txt", sep="\t")
    #
    df_pseudo = pd.read_csv("data/pseudo.allele.id.list.txt", sep="\t", header=None)
    p_list = df_pseudo[0].tolist()

    # read the clonefl and obtain vbesthit, cbesthit info
    df_clone = pd.read_csv(clonefl, sep="\t", usecols=["cloneId", "allVHitsWithScore", "allCHitsWithScore"])
    if df_clone.dropna(subset=["allCHitsWithScore"]).size == 0:
        os._exit(0)
    else:
        df_clone["bestVHit"] = df_clone.allVHitsWithScore.str.split("(", expand=True)[0]
        df_clone["bestCHit"] = df_clone.allCHitsWithScore.str.split("(", expand=True)[0]
    null_cloneId_list = df_clone.cloneId[df_clone.bestCHit.isnull()].tolist()
    # extract reads coming from a certain isotype
    df_clone2 = df_clone[df_clone["bestCHit"].str.contains(iso, na=False)]
    # remove pseudoallele containing clone
    df_clone2 = df_clone2[~df_clone2.bestVHit.isin(p_list)]
    # if there is no records remained, the script will exit
    if df_clone2.size == 0:
        os._exit(0)
        
        
    # read and filter alignment file for available reads for statistics
    df_align = pd.read_csv(alignfl, sep="\t") 
    # remove those clone with no isotype annotation
    df_align = df_align[~df_align.cloneId.isin(null_cloneId_list)]
    # extracted reads that can be merged
    df_align = df_align[~df_align.refPoints.str.contains(",", na=False)]
    # extracted reads that are with full-length variable region
    df_align = df_align.dropna(subset=["nSeqFR1", "nSeqCDR1", "nSeqFR2", "nSeqCDR2", "nSeqFR3", "nSeqCDR3", "nSeqFR4"])
    df_align = df_align.dropna(subset=["aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4"])
    # extracted reads with indel
    df_align = df_align[~df_align.bestVAlignment.str.contains('[DI]', regex=True)]
    # drop reads with different v and c allele annotation with its clone's
    s_id_allele = df_clone.bestVHit
    s_id_allele.index = df_clone.cloneId.astype(int).tolist()
    valleleDict = s_id_allele.to_dict()
    valleleDict.update({-1:""})
    vallele_list = [valleleDict[x] for x in df_align.cloneId.astype(int).tolist()]
    vallele_series = pd.Series(vallele_list)
    df_align.reset_index(inplace=True)
    df_align = df_align[df_align.bestVHit == vallele_series]  

    # assign the df_clone2 to df_clone
    # this step is for circumenting the keyerror raised in line 114
    df_clone = df_clone2
    
    # extracted reads that are without frame-shift and stop codon mutations (assumed to be productive)
    df_align["aaSeq"] = (df_align.aaSeqFR1 + df_align.aaSeqCDR1 + df_align.aaSeqFR2 + df_align.aaSeqCDR2 + df_align.aaSeqFR3 + df_align.aaSeqCDR3 + df_align.aaSeqFR4)
    df_align["status"] = df_align.aaSeq.apply(is_productive)
    if df_align.size == 0:
        os._exit(0)
    df_align = df_align[df_align.status=="Productive"]

    # remove reads with same variable regions but were assigned to different clones
    df_align["VR"] = df_align["nSeqFR1"]+df_align["nSeqCDR1"]+df_align["nSeqFR2"]+df_align["nSeqCDR2"]+df_align["nSeqFR3"]+df_align["nSeqCDR3"]+df_align["nSeqFR4"]
    vr_list = df_align["VR"].tolist()
    cloneId_list = df_align["cloneId"].tolist()
    vr_dict = {}
    sense_list = []
    for vr, cloneId in zip(vr_list, cloneId_list):
        if vr not in vr_dict:
            vr_dict.update({vr:cloneId})
            sense_list.append(True)
        else: 
            if vr_dict[vr] != cloneId:
                sense_list.append(False)
            else:
                sense_list.append(True)
    df_align = df_align[sense_list]

    # perform statistics analysis of the number of productive and
    # unproductive non-deduplicated reads for each clone
    d = {}
    for name, group in df_align.status.groupby(df_align.cloneId):
        p = group.tolist().count("Productive")
        up = group.tolist().count("Unproductive")
        d.update({str(name):[p, up]})

    # deduplicate reads
    df_align = df_align[~(df_align.nSeqFR1 + df_align.nSeqCDR1 + df_align.nSeqFR2 + df_align.nSeqCDR2 + \
    df_align.nSeqFR3 + df_align.nSeqCDR3 + df_align.nSeqFR4).duplicated()]
    
    # iterate over clone to extracted qualified reads for statistical analysis
    # for clone with remaining productive reads, the clone is regarded as prod
    # -uctive, and only productive reads were enrolled for analysis.
    # On the contary, if the remaining reads for 
    for cloneId in df_clone["cloneId"]:
        # obtain v allele hit
        vallele = df_clone[df_clone["cloneId"]==cloneId]["bestVHit"].tolist()[0]
        # obtain isotype
        isotype = df_clone[df_clone["cloneId"]==cloneId]["bestCHit"].tolist()[0]
        # obtain the number of productive and unproductive reads
        np, nu = (0,0) if str(cloneId) not in d.keys() else d[str(cloneId)]
        # extract reads that were assigned to a certain cloneId
        df_cloneId = df_align[df_align["cloneId"]==cloneId]
        # determine if the reads for a clone has been fully fitered out
        if df_cloneId.size != 0:
            # try to extract productive reads
            df_prod = df_cloneId[df_cloneId["status"]=="Productive"]
            # determine if exist productive reads for such clone
            if df_prod.shape[0] != 0 and np > 1:
                aligntagList = [x.split("|")[5] for x in df_prod.bestVAlignment.tolist()]
                mutation_tag_processing(vallele, aligntagList, alleleDictProd, df_delim)
            df_unprod = df_cloneId[df_cloneId["status"]=="Unproductive"]
            if df_unprod.shape[0] and nu > 1:
                aligntagList = [x.split("|")[5] for x in df_unprod.bestVAlignment.tolist()]
                mutation_tag_processing(vallele, aligntagList, alleleDictUnprod, df_delim)

    #### output module ####
    for allele in alleleDictProd.keys():
        allele_new_name = allele.replace('/', '.').replace('*', '.')
        df_allele = pd.DataFrame(alleleDictProd[allele])[["pos", "ref", "A", "T", "C", "G", "is_synA", "is_synT", "is_synC", "is_synG", "nReads"]]
        df_allele["A"] = df_allele["A"].apply(lambda x:"%.3f"%x)
        df_allele["T"] = df_allele["T"].apply(lambda x:"%.3f"%x)
        df_allele["C"] = df_allele["C"].apply(lambda x:"%.3f"%x)
        df_allele["G"] = df_allele["G"].apply(lambda x:"%.3f"%x)
        df_allele.to_csv("%s/productive/%s-%s/%s.pos.mut.type.stat"%(outdir,sample,iso,allele_new_name), sep='\t', index=0)
    for allele in alleleDictUnprod.keys():
        allele_new_name = allele.replace('/', '.').replace('*', '.')
        df_allele = pd.DataFrame(alleleDictUnprod[allele])[["pos", "ref", "A", "T", "C", "G", "is_synA", "is_synT", "is_synC", "is_synG", "nReads"]]
        df_allele["A"] = df_allele["A"].apply(lambda x:"%.3f"%x)
        df_allele["T"] = df_allele["T"].apply(lambda x:"%.3f"%x)
        df_allele["C"] = df_allele["C"].apply(lambda x:"%.3f"%x)
        df_allele["G"] = df_allele["G"].apply(lambda x:"%.3f"%x)
        df_allele.to_csv("%s/unproductive/%s-%s/%s.pos.mut.type.stat"%(outdir,sample,iso,allele_new_name), sep='\t', index=0)
                           
                           
if __name__ == "__main__":
    sample = sys.argv[1] #"Healthy_00001"
    alignfl =  sys.argv[2] #"alignments.txt"
    clonefl =  sys.argv[3] #"clones.txt"
    iso = sys.argv[4]
    outdir = sys.argv[5]
    main()
