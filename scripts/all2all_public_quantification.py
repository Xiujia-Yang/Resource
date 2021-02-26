#!/usr/bin/env python

import csv
import pandas as pd

def cloneTagSampleDictAssign(cloneTagList, sampleId):
	for cloneTag in cloneTagList:
		try:
			cloneTagSampleDict[cloneTag].append(sampleId)
		except:
			cloneTagSampleDict.update({cloneTag:[sampleId]})

def main():
	for sampleId in df.index:
		print sampleId
		cloneflPath = df.loc[sampleId,"cloneflPath"]
		df_clone = pd.read_csv(cloneflPath, sep="\t", usecols=["targetSequences", \
		"allVHitsWithScore","allJHitsWithScore","aaSeqCDR3"])
		df_clone = df_clone[(~df_clone["aaSeqCDR3"].str.contains("\*")) & (~df_clone["aaSeqCDR3"].str.contains("_"))]
		cloneTagList = df_clone["allVHitsWithScore"].str.split(pat="*",expand=True)[0] + "_" +\
						df_clone["allJHitsWithScore"].str.split(pat="*",expand=True)[0] + "_" + \
						df_clone["targetSequences"]
		cloneTagList = cloneTagList[cloneTagList.notna()]
		cloneTagSampleDictAssign(cloneTagList, sampleId)

	# Initialize a dictionary object to contain the number of shared
	# clone between samples
	pubCloneDict = {}
	for sampleIdA in df.index:
		for sampleIdB in df.index:
			pubCloneDict.update({(sampleIdA,sampleIdB):0})
			
	# Iterate over cloneTagSampleDict to count the number of sample 
	for cloneTag, sampleIdList in cloneTagSampleDict.items():
		sampleIdSet = set(sampleIdList)
		sampleIdList = list(sampleIdSet)
		sampleIdListStr = ",".join(sampleIdList)
		pubStat1.writerow([cloneTag, sampleIdListStr, len(sampleIdList)])
		for sampleA in sampleIdList:
			for sampleB in sampleIdList:
				pubCloneDict[(sampleA,sampleB)]+=1
	
	pubStat2.writerow([""]+list(df.index))
	pubStat3.writerow([""]+list(df.index))
	for sampleA in df.index:
		sampleAProj = df.loc[sampleA,"project"]
		outLine2 = [sampleA] # original statistics
		outLine3 = [sampleA] # chimera-removed statistics
		for sampleB in df.index:
			sampleBProj = df.loc[sampleB,"project"]
			outLine2.append( 0 if sampleB == sampleA else pubCloneDict[(sampleA,sampleB)])
			outLine3.append( 0 if sampleBProj == sampleAProj else pubCloneDict[(sampleA,sampleB)] )
		pubStat2.writerow(outLine2)
		pubStat3.writerow(outLine3)

if __name__ == "__main__":
	cloneTagSampleDict = {}
	#df = pd.read_csv("cloneflpath", sep="\t", index_col=0)
	df = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
	pubStat2 = csv.writer(open("PublicClone_VJCDR3_Stat_Original.txt", "wb"), delimiter="\t")
	pubStat3 = csv.writer(open("PublicClone_VJCDR3_Stat_RmChimera.txt", "wb"), delimiter="\t")
	pubStat1 = csv.writer(open("Clone_VJCDR3_SharedSample_SampleNumber.txt", "wb"), delimiter="\t")
	main()
