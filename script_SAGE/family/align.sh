#!/bin/bash

## Alignement of sequence for each gene family using MUSCLE
## Save file into as .phy (phyml compatible format)

cd /scratch/beegfs/monthly/mls_2016/group6/test/family_sequence
module add SequenceAnalysis/MultipleSequenceAlignment/muscle/3.8.31

for fastafile in *.uni.txt
	do var=`echo $fastafile |cut -f 1 -d '.'`
	 muscle -in $fastafile -phyiout $var.phy
done
