#!/usr/bin/env python

'''
    This script is for perform basic statistic analysis
    of somatic hypermutation for each sample. The invol
    -lved algorthim is based on position weight matrix,
    which considered all available reads within a clone.
    
    Usage: python regional_mut_freq_for_each_sample.py sample alignment.txt clones.txt
'''

import csv, sys, re, os
import pandas as pd
import numpy as npy

def mutRateStat(vallele, alleleInfoDf, aligntagList):
    if aligntagList == []:
        return ["None", "None", "None", "None", "None", "None"]
    else:
        # obtain the number of valid reads
        nReads = len(aligntagList)
        # define a variable to record the number of point mutations
        # for all valid reads
        mutNumSum, fr1Sum, cdr1Sum, fr2Sum, cdr2Sum, fr3Sum = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        for aligntag in aligntagList:
            # obtain the positions of point mutation
            posList = [int(x[2:-1])+1 for x in re.findall(r'S[ATCG]\d+[ATCG]', aligntag)]
            # obtain mutation number for each valid reads
            mutNum = len(posList)
            mutNumSum += mutNum
            fr1, cdr1, fr2, cdr2, fr3 = regionMutCount(vallele, posList, alleleInfoDf)
            fr1Sum+=fr1; cdr1Sum+=cdr1; fr2Sum+=fr2; cdr2Sum+=cdr2; fr3Sum+=fr3
        # obtain the length
        _, fr1Begin, cdr1Begin, fr2Begin, cdr2Begin, fr3Begin, cdr3Begin, _ = alleleInfoDf[alleleInfoDf["vAllele"]==vallele].values.tolist()[0]
        fr1Len = cdr1Begin-fr1Begin
        cdr1Len = fr2Begin-cdr1Begin
        fr2Len = cdr2Begin-fr2Begin
        cdr2Len = fr3Begin-cdr2Begin
        fr3Len = cdr3Begin-fr3Begin 
        # calculated the average mutation rate across all region
        alleleLen = cdr3Begin
        mr = mutNumSum/(alleleLen*nReads)
        # obtain the mutation rate
        fr1mr = "NA" if fr1Len==0 else fr1Sum/(fr1Len*nReads)
        cdr1mr = cdr1Sum/(cdr1Len*nReads)
        fr2mr = fr2Sum/(fr2Len*nReads)
        cdr2mr = cdr2Sum/(cdr2Len*nReads)
        fr3mr = fr3Sum/(fr3Len*nReads)
        
        # mr means mutation rate
        return [mr, fr1mr, cdr1mr, fr2mr, cdr2mr, fr3mr]

def regionMutCount(vallele, posList, alleleInfoDf):
    '''
        This function give statistics of the number of 
        mutation for each functional region
    '''
    fr1, cdr1, fr2, cdr2, fr3 = [0, 0, 0, 0, 0]
    fr1Begin, cdr1Begin, fr2Begin, cdr2Begin, fr3Begin, _, seq = alleleInfoDf[alleleInfoDf["vAllele"]==vallele].values.tolist()[0][1:]
    for pos in posList:
        if pos < cdr1Begin:
            fr1 += 1
        elif pos < fr2Begin:
            cdr1 += 1
        elif pos < cdr2Begin:
            fr2 += 1
        elif pos < fr3Begin:
            cdr2 += 1
        else:
            fr3 += 1
    return [fr1, cdr1, fr2, cdr2, fr3]

def is_productive(s):
    if s[-1] == "_":
        s = s[:-1]
    if s.count('*') == 0 and s.count('_') == 0:
        return "Productive"
    else:
        return "Unproductive"

def make_sense(s):
    cloneId = s["cloneId"]
    VR = s["VR"]
    
        

def main():
    # read the allele region 
    df_delim = pd.read_csv("data/v.allele.region.delimitation.txt", sep="\t")

    # obtain pseudo allele list
    df_pseudo = pd.read_csv("data/pseudo.allele.id.list.txt", sep="\t", header=None)
    p_list = df_pseudo[0].tolist()
    
    # read the clonefl and obtain vbesthit, cbesthit info
    df_clone = pd.read_csv(clonefl, sep="\t", usecols=["cloneId", "cloneCount", "allVHitsWithScore", "allCHitsWithScore"])
    df_clone["bestVHit"] = df_clone.allVHitsWithScore.str.split("(", expand=True)[0]
    # deal with the exception that all clone were without C alignment info
    if df_clone[~df_clone.allCHitsWithScore.isnull()].shape[0] != 0:
        df_clone["bestCHit"] = df_clone.allCHitsWithScore.str.split("(", expand=True)[0]
    else:
        df_clone["bestCHit"] = npy.nan
    null_cloneId_list = df_clone.cloneId[df_clone.bestCHit.isnull()].tolist()
    s_id_allele = df_clone.bestVHit
    s_id_allele.index = df_clone.cloneId.astype(int).tolist()
    valleleDict = s_id_allele.to_dict()
    valleleDict.update({-1:""})
    df_clone = df_clone[~df_clone.bestVHit.isin(p_list)]
    
    # read and filter alignment file for available reads for statistics
    df_align = pd.read_csv(alignfl, sep="\t") 
    # remove those clone with no isotype annotation
    df_align = df_align[~df_align.cloneId.isin(null_cloneId_list)]
    # deal with the exception that some records is without refPoints info
    df_align = df_align[~df_align.refPoints.isnull()]
    # extracted reads that can be merged
    df_align = df_align[~df_align.refPoints.str.contains(",")]
    # extracted reads that are with full-length variable region
    df_align = df_align.dropna(subset=["nSeqFR1", "nSeqCDR1", "nSeqFR2", "nSeqCDR2", "nSeqFR3", "nSeqCDR3", "nSeqFR4"])
    # extracted reads without indel
    df_align = df_align[~df_align.bestVAlignment.str.contains('[DI]', regex=True)]
    # drop reads with different v and c allele annotation with its clone's
    vallele_list = [valleleDict[x] for x in df_align.cloneId.astype(int).tolist()]
    vallele_series = pd.Series(vallele_list)
    df_align.reset_index(inplace=True)
    df_align = df_align[df_align.bestVHit == vallele_series]  
    
    # extracted reads that are without frame-shift and stop codon mutations (assumed to be productive)
    df_align["aaSeq"] = (df_align.aaSeqFR1 + df_align.aaSeqCDR1 + df_align.aaSeqFR2 + df_align.aaSeqCDR2 + df_align.aaSeqFR3 + df_align.aaSeqCDR3 + df_align.aaSeqFR4)
    df_align["status"] = df_align.aaSeq.apply(is_productive)
    df_align = df_align[df_align.status=="Productive"]
    if df_align.size == 0:
        os._exit(0)

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
    d = {} # define a dictionary 
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
        # obtain clone count
        clonecount = df_clone[df_clone["cloneId"]==cloneId]["cloneCount"].tolist()[0]
        # extract reads that were assigned to a certain cloneId
        df_cloneId = df_align[df_align["cloneId"]==cloneId]
        # obtain the number of productive and unproductive reads
        np, nu = (0,0) if str(cloneId) not in d.keys() else d[str(cloneId)]
        # determine if the reads for a clone has been fully fitered out
        if df_cloneId.size == 0:
            prod = "unknown"
            nAvaiReads = 0
            outHandle.writerow([cloneId, clonecount, vallele, isotype, prod] + mutRateStat(vallele, df_delim, []) + [nAvaiReads, np, nu])
        else:
            df_prod = df_cloneId[df_cloneId.status=="Productive"]
            if df_prod.shape[0] != 0:
                prod = "productive"
                nAvaiReads = df_prod.shape[0]
                aligntagList = [x.split("|")[5] for x in df_prod.bestVAlignment.tolist()]
                outHandle.writerow([cloneId, clonecount, vallele, isotype, prod] + mutRateStat(vallele, df_delim, aligntagList) + [nAvaiReads, np, nu])
            df_unprod = df_cloneId[df_cloneId.status=="Unproductive"]
            if df_unprod.shape[0] != 0:
                prod = "unproductive"
                nAvaiReads = df_unprod.shape[0]
                aligntagList = [ x.split("|")[5] for x in df_unprod.bestVAlignment.tolist() ]
                outHandle.writerow([cloneId, clonecount, vallele, isotype, prod] + mutRateStat(vallele, df_delim, aligntagList) + [nAvaiReads, np, nu])

if __name__ == "__main__":
    sample = sys.argv[1]
    alignfl = sys.argv[2]
    clonefl = sys.argv[3]
    f = open(sample+"_mut_freq_stat.txt", "w")
    outHandle = csv.writer(f, delimiter="\t")
    outHandle.writerow(["cloneId", "cloneCount", "vAllele", "isotype", "productivity", "allMutRate", "fr1MutRate", \
                    "cdr1MutRate", "fr2MutRate", "cdr2MutRate", "fr3MutRate", "nAvaiReads", "nProdReads", "nUnprodReads"])
    main()
    f.close()
