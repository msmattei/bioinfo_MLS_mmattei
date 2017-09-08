#!/bin/bash

## Remove gene and gene sequence of singular genes multiple times into the same gene family
## Output file = *.uni.txt used to create the *.phy files

cd /scratch/beegfs/monthly/mls_2016/group6/test/family_sequence
for file in *.txt
	do var=`echo $file |cut -f 1 -d '.'`
	 awk '/^>/{f=!d[$1];d[$1]=1}f' $file > $var.uni.txt
	
done

