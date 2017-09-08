#!/bin/bash

## Remove genes present twice and its correspondent sequence.

cd /scratch/beegfs/monthly/mls_2016/group6/test/family_sequence/temp
for file in *.txt
	do var=`echo $file |cut -f 1 -d '.'`
	 awk '/^>/{f=!d[$1];d[$1]=1}f' $file > $var.uni.txt
	
done

