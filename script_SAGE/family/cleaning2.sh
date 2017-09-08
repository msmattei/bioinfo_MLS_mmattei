#!/bin/bash

## remove all files containing less or equel then 3 genes:
cd /scratch/beegfs/monthly/mls_2016/group6/test/family_sequence

for file in *.txt
	do n=`grep ">" $file |wc -l`
	var=`echo $file |cut -f 1 -d '.'`
	if [ "$n" -le 3 ]
	then mv $var.phy short/
	mv $file short/
	fi 
done