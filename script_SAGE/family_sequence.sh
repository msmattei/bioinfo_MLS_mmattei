#!/bin/bash

## Find nucleotide sequence from concatenated fasta file (all.headers_concat.ffn) for each member of all gene family
## Save each file in the family_sequence folder

cd /scratch/beegfs/monthly/mls_2016/group6/test/family

# First loop going through file containing the list of gene in a family
# second loop going through the gene in the file (gene family)

for file in *.txt
	do for gene in `cat "$file"`
		do grep -w -A1 "$gene" ../all.headers_concat.ffn
	 done > ../family_sequence/$file
 done
