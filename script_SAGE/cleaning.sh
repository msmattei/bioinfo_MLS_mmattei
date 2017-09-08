#!/bin/bash

## find all identical gene family files and move them to the identical folder
## create a file identica.txt containing the list of all identical family file (used to remove them from the ../family_sequence folder)

cd /scratch/beegfs/monthly/mls_2016/group6/test/family

for file1 in *.txt; do 
	for file2 in *.txt; do
		if [ $file1 != $file2 ]
		then
			if cmp -s "$file1" "$file2"
			then
				mv $file2  identical/
				echo "$file2" >> identical/identical.txt
			fi
		fi
	done
done
