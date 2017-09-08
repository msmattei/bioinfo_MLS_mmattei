#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e error.txt
#BSUB -n 2

#for i in *.faa; do hmmbuild ${i%.faa}.hmm $i; hmmsearch --tblout ${i%.faa}.tblout -A ${i%.faa}.hmmaln --notextw ${i%.faa}.hmm ./all_firm5_genomes.faa; esl-reformat a2m ${i%.faa}.hmmaln > ${i%.faa}.a2m; done

for i in *.faa; do hmmsearch --tblout ${i%.faa}.tblout -A ${i%.faa}.hmmaln --notextw ${i%.faa}.hmm ./all_firm5_genomes.faa; esl-reformat a2m ${i%.faa}.hmmaln > ${i%.faa}.a2m; done

#for i in *hmm; do esl-reformat a2m ${i%.hmm}.hmmaln > ${i%.hmm}.a2m; done
