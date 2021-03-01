#!/usr/bin/env python

'''
    This script is for calculating the frequency and
    mutation rate of different motifs, and the percent
    of target nucleotides.
    
    Usage: python motif_mut_freq_cal_for_single_sample.py sample alignments.txt clones.txt
'''

import csv, sys, re, os
import pandas as pd
import numpy as npy
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def extract_codon_and_variants(seq, cdr3Begin, locus):
    nt_list = [x for x in ["A", "C", "G", "T"] if x in list(set(["A", "C", "G", "T"])-set([seq[locus]]))]
    if ((locus+1)-cdr3Begin%3)%3 == 0:
        codon = seq[locus-2:locus+1]
        codon2, codon3, codon4 = [seq[locus-2]+seq[locus-1]+nt for nt in nt_list]
    elif ((locus+1)-cdr3Begin%3)%3 == 1:
        if locus+2 <= cdr3Begin:
            codon = seq[locus:locus+3]
            codon2, codon3, codon4 = [nt+seq[locus+1]+seq[locus+2] for nt in nt_list]
        else:
            return 0
    elif ((locus+1)-cdr3Begin%3)%3 == 2:
        if locus+1 <= cdr3Begin:
            codon = seq[locus-1:locus+2]
            codon2, codon3, codon4 = [seq[locus-1]+nt+seq[locus+1] for nt in nt_list]
        else:
            return 0

    aa = Seq(codon, generic_dna).translate()
    aa2 = Seq(codon2, generic_dna).translate()
    aa3 = Seq(codon3, generic_dna).translate()
    aa4 = Seq(codon4, generic_dna).translate()
    
    bar=""
    for AA in [aa2, aa3, aa4]:
        if AA != aa:
            bar=bar+"R"
        else:
            bar=bar+"S"
    return bar


def is_productive(s):
    if s[-1] == "_":
        s = s[:-1]
    if s.count('*') == 0 and s.count('_') == 0:
        return "Productive"
    else:
        return "Unproductive"


def motif_cal(vallele, aligntaglist): 
    if aligntaglist != []:
        # obtain the number of valid reads
        nReads = float(len(aligntaglist))
        '''
        # define a variable to record the number of point mutations
        # for all valid reads
        mutNumSum, fr1Sum, cdr1Sum, fr2Sum, cdr2Sum, fr3Sum = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        '''
        for aligntag in aligntaglist:
            # obtain the positions of point mutation
            posTarList = [[int(x[2:-1]), x[-1]] for x in re.findall(r'S[ATCG]\d+[ATCG]', aligntag)]
            
            
            posTarListOver = [x for x in posTarList if x[0] in wrcy[vallele]]
            mtf["WRCY"][0] += len(wrcy[vallele])/nReads
            mtf["WRCY"][1] += len(posTarListOver)/nReads
            mtf["WRCY"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["WRCY"][4] += len([x for x in posTarListOver if x[1] == "G"])/nReads
            mtf["WRCY"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in rgyw[vallele]]
            mtf["RGYW"][0] += len(rgyw[vallele])/nReads
            mtf["RGYW"][1] += len(posTarListOver)/nReads
            mtf["RGYW"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["RGYW"][3] += len([x for x in posTarListOver if x[1] == "C"])/nReads
            mtf["RGYW"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in wa[vallele]]
            mtf["WA"][0] += len(wa[vallele])/nReads
            mtf["WA"][1] += len(posTarListOver)/nReads
            mtf["WA"][3] += len([x for x in posTarListOver if x[1] == "C"])/nReads
            mtf["WA"][4] += len([x for x in posTarListOver if x[1] == "G"])/nReads
            mtf["WA"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in tw[vallele]]
            mtf["TW"][0] += len(tw[vallele])/nReads
            mtf["TW"][1] += len(posTarListOver)/nReads
            mtf["TW"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["TW"][4] += len([x for x in posTarListOver if x[1] == "G"])/nReads
            mtf["TW"][3] += len([x for x in posTarListOver if x[1] == "C"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in syc[vallele]]
            mtf["SYC"][0] += len(syc[vallele])/nReads
            mtf["SYC"][1] += len(posTarListOver)/nReads
            mtf["SYC"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["SYC"][4] += len([x for x in posTarListOver if x[1] == "G"])/nReads
            mtf["SYC"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in grs[vallele]]
            mtf["GRS"][0] += len(grs[vallele])/nReads
            mtf["GRS"][1] += len(posTarListOver)/nReads
            mtf["GRS"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["GRS"][3] += len([x for x in posTarListOver if x[1] == "C"])/nReads
            mtf["GRS"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in a[vallele]]
            mtf["A"][0] += len(a[vallele])/nReads
            mtf["A"][1] += len(posTarListOver)/nReads
            mtf["A"][3] += len([x for x in posTarListOver if x[1] == "C"])/nReads
            mtf["A"][4] += len([x for x in posTarListOver if x[1] == "G"])/nReads
            mtf["A"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in c[vallele]]
            mtf["C"][0] += len(c[vallele])/nReads
            mtf["C"][1] += len(posTarListOver)/nReads
            mtf["C"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["C"][4] += len([x for x in posTarListOver if x[1] == "G"])/nReads
            mtf["C"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in g[vallele]]
            mtf["G"][0] += len(g[vallele])/nReads
            mtf["G"][1] += len(posTarListOver)/nReads
            mtf["G"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["G"][3] += len([x for x in posTarListOver if x[1] == "C"])/nReads
            mtf["G"][5] += len([x for x in posTarListOver if x[1] == "T"])/nReads
            posTarListOver = [x for x in posTarList if x[0] in t[vallele]]
            mtf["T"][0] += len(t[vallele])/nReads
            mtf["T"][1] += len(posTarListOver)/nReads
            mtf["T"][2] += len([x for x in posTarListOver if x[1] == "A"])/nReads
            mtf["T"][4] += len([x for x in posTarListOver if x[1] == "G"])/nReads
            mtf["T"][3] += len([x for x in posTarListOver if x[1] == "C"])/nReads
            
    

def main():
    # obtain pseudo allele list
    df_pseudo = pd.read_csv("/zzh_gpfs02/yangxiujia/20190128-Resource-MiXCR/18-Somatic-hypermutation/06-positional-mutRate-v3/1.pseudo.allele.real.name.txt", sep="\t", header=None)
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
    #df_clone["bestCHit"] = df_clone["bestCHit"].astype(str).apply(lambda x:x if x.startswith("IGH") else "None")
    
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
    # extract reads coming from a certain isotype
    df_clone = df_clone[df_clone["bestCHit"].str.contains("IGHG", na=False)]
    # if there is no records remained, the script will exit
    if df_clone.size == 0:
        os._exit(0)
    
    # extracted reads that are without frame-shift and stop codon mutations (assumed to be productive)
    df_align["aaSeq"] = (df_align.aaSeqFR1 + df_align.aaSeqCDR1 + df_align.aaSeqFR2 + df_align.aaSeqCDR2 + df_align.aaSeqFR3 + df_align.aaSeqCDR3 + df_align.aaSeqFR4)
    df_align["status"] = df_align.aaSeq.apply(is_productive)
    df_align = df_align[df_align.status=="Productive"]
    # if there is no records remained, the script will exit
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
    for name, group in df_align["status"].groupby(df_align.cloneId):
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
    n = 0
    for cloneId in df_clone["cloneId"]:
        # obtain v allele hit
        vallele = df_clone[df_clone["cloneId"]==cloneId]["bestVHit"].tolist()[0]
        # obtain isotype
        isotype = df_clone[df_clone["cloneId"]==cloneId]["bestCHit"].tolist()[0]
        # extract reads that were assigned to a certain cloneId
        df_cloneId = df_align[df_align["cloneId"]==cloneId]
        # obtain the number of productive and unproductive reads
        np, nu = (0,0) if str(cloneId) not in d.keys() else d[str(cloneId)]
        # determine if the reads for a clone has been fully fitered out
        if df_cloneId.size != 0 and np > 1 :
            n += 1
            #print cloneId
            df_prod = df_cloneId[df_cloneId.status=="Productive"]
            if df_prod.shape[0] != 0:
                aligntagList = [x.split("|")[5] for x in df_prod.bestVAlignment.tolist()]
                motif_cal(vallele, aligntagList)
    # writer the result
    outHandle = csv.writer(open(outfl, "wb"), delimiter="\t")
    outHandle.writerow(["motif","total","mut","A","C","G","T"])
    motif_list = ["WRCY", "RGYW", "WA", "TW", "SYC", "GRS", "A", "C", "G", "T"]
    for motif in motif_list:
        outHandle.writerow([motif]+mtf[motif])
        
if __name__ == "__main__":
    sample = sys.argv[1]
    alignfl = sys.argv[2]
    clonefl = sys.argv[3]
    outfl = "%s_motif_stat.txt"%sample

    # define a dictionary to contain the statistic
    mtf = {"WRCY":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
           "RGYW":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "WA":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "TW":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "SYC":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "GRS":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "A":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "C":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "G":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "T":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}

    # define the pattern
    reDict = {
        "p_wrcy" : re.compile("(?=([AT][AG]C[CT]))"),
        "p_rgyw" : re.compile("(?=([AG]G[CT][AT]))"),
        "p_wa" : re.compile("(?=([AT]A))"),
        "p_tw" : re.compile("(?=(T[AT]))"),
        "p_syc" : re.compile("(?=([CG][CT]C))"),
        "p_grs" : re.compile("(?=(G[AG][CG]))"),
        "p_a" : re.compile("A"),
        "p_t" : re.compile("T"),
        "p_c" : re.compile("C"),
        "p_g" : re.compile("G")
    }

    # Import the reference sequence and region
    # delimitation info, and obtain the SSS type
    # motif/nucleotide locus
    locusDict = {}
    df = pd.read_csv("/zzh_gpfs02/yangxiujia/20190128-Resource-MiXCR/23-CR-Revision/02-SHM/reference_fasta_format/v.allele.region.delimitation.txt", sep="\t", index_col=0)
    for allele in df.index:
        CDR3Begin = df.loc[allele, "CDR3Begin"]
        seq = df.loc[allele, "Sequence"]
        for pat in reDict.keys():
            motif = pat[2:]
            if not locusDict.has_key(motif):
                locusDict.update({motif:{}})
            if motif in ["wrcy", "syc"]:
                offset = 2
            elif motif in ["rgyw", "wa"]:
                offset = 1
            else:
                offset = 0

            locusList = []
            for locus in [x+offset for x in [m.start() for m in reDict[pat].finditer(seq)] if (x+offset) < CDR3Begin]:
                typ = extract_codon_and_variants(seq, CDR3Begin, locus)
                if typ == "SSS":
                    locusList.append(locus)
                    #print locus
            locusDict[motif].update({allele:locusList})

    # Variable assignment for adaption to previous script
    wrcy = locusDict['wrcy']
    rgyw = locusDict['rgyw']
    syc = locusDict['syc']
    grs = locusDict['grs']
    wa = locusDict['wa']
    tw = locusDict['tw']
    a = locusDict['a']
    t = locusDict['t']
    c = locusDict['c']
    g = locusDict['g']

    # remove the overlapping locus in single nucleotide
    for allele in df.index:
        a[allele] = list(set(a[allele]) - set(wa[allele]))
        t[allele] = list(set(t[allele]) - set(tw[allele]))
        c[allele] = list(set(c[allele]) - set(wrcy[allele]) - set(syc[allele]))
        g[allele] = list(set(g[allele]) - set(rgyw[allele]) - set(grs[allele]))
        
    main()
