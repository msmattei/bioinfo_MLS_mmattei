#!/bin/bash
##

## For each HGT candidate find all gene belonging to the same gene family 
## (using the mclOutput in the orthoMCL_output folder)

#change directory
cd /scratch/beegfs/monthly/mls_2016/group6/test

for gene in `cat hgt_candidate20170516_list.txt`; 
	do let "var++"; grep -w "$gene" ../../orthoMCL_output/mclOutput | tr '\t' '\n' > family/family$var.txt; 
done

