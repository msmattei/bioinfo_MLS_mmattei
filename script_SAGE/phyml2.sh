#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o ./phyml.out
#BSUB -e /phyml.err
#BSUB -M 10000000

module load Phylogeny/phyml/3.3.20170119
cd /scratch/beegfs/monthly/mls_2016/group6/test/family_sequence/temp
for i in 1
        do file=`sed -n ''$i'p' < best_method_result.txt | awk '{print $1}'`
        method=`sed -n ''$i'p' < best_method_result.txt | awk '{print $2}' |cut -f 1 -d '+'`
        met=`sed -n ''$i'p' < best_method_result.txt | awk '{print $2}'`
        if [ `echo $met | grep "+G" -c` -ne 0 ]
        then if [ `echo $met | grep "+I" -c` -ne 0 ]
                then phyml -i $file -m $method -b 100 -v e -a e
                else phyml -i $file -m $method -b 100 -a e
                fi
        else if [ `echo $met | grep "+I" -c` -ne 0 ]
                then phyml -i $file -m $method -b 100 -v e
                else phyml -i $file -m $method -b 100
                fi
        fi
done

