# Resource

## Overview
This repository contains all scripts and preliminary statistics related to the key analyses in each section in Result of manuscript (MS) ***"Large-scale analysis of 2,152 Ig-seq datasets reveals key features of B cell biology and the antibody repertoire"***. These analyses can be classified into four categories, including gene usage, somatic recombination, somatic hypermutation and public clones. More details can be found in the following sections.

## Gene usage


## Somatic recombination


## Somatic hypermutation


## Public clone
### Preliminary statistics
* [PUB_all_to_all_repertoire_comparison_matrix.csv.gz](data/PUB_all_to_all_repertoire_comparison_matrix.csv.gz) 
	- The 2152*2152 matrix recording all-to-all repertoire comparisons (Fig. 5A in MS)
* [PUB_pairwise_repertoire_comparison.tab.gz](data/PUB_pairwise_repertoire_comparison.tab.gz) 
	- The table recording pairwise repertoire comparison together with the number of total clones for compared samples (Fig. 5B) 
* [PUB_public_clone_record_from_mixcr_output.tab.gz](data/PUB_public_clone_record_from_mixcr_output.tab.gz) 
	- The file containing all lines in the clones.txt file that records public clones (Fig 5C)
### Scripts
* All-to-all public clone quantification (Fig. 5A)
`python all2all_public_quantification.py clone_file_path.tab`
This python script implements the all-to-all public clone quantification among enrolled samples. It takes an annotation tabular file `clone_file_path.tab` as input. This file has three columns including sample id, project id and path to `clones.txt` that output by MiXCR and looks like below,
	```
	Sample	Project	Path
	DRR056252	PRJDB4353	DRR056252_MIXCR/clones.txt
	DRR056253	PRJDB4353	DRR056253_MIXCR/clones.txt
	DRR056254	PRJDB4353	DRR056254_MIXCR/clones.txt
	ERR1760498	PRJEB15295	ERR1760498_MIXCR/clones.txt
	ERR1812282	PRJEB18926	ERR1812282_MIXCR/clones.txt
	ERR1812283	PRJEB18926	ERR1812283_MIXCR/clones.txt
	```
* All-to-all public clone visualization (Fig. 5A)
* Linear model visualization (Fig. 5B)
* Quantification and visualization of the correlation between clonality and publicness of public clones (Fig. 5C)

## Dependencies
In-house scripts above were written in Python (v3.7) and MATLAB (v). For python, a series of modules are required, they include pandas, csv. For MATLAB, xxx is required. 