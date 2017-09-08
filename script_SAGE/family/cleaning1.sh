#!/bin/bash

## Remove all identical gene family files (detected usig
cd /scratch/beegfs/monthly/mls_2016/group6/test/family_sequence
for file in `cat identical/identical.txt`
	do var=`echo $file |cut -f 1 -d '.'`
	mv $file identical/
	mv $var.phy identical/
done

