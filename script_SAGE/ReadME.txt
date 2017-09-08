# To linearize the sequences (multiline-FASTA to oneline-FASTA)
for i in *gbk.faa; do perl -pe '/^>/ ? print "\n" : chomp' $i > ${i%_headers.gbk.faa}_linear.gbk.faa; done

# To change the headers and get rid of whitespaces, otherwise they are cut at the genome name:
foreach i in *faa; do cat $i | tr " " "_" > ${i%..gbk.faa}_tr.gbk.faa; done

# merge them all together into a single file
cat *linear*tr.gbk.faa > one_file_to_rule_them_all_tr.gbk.faa

# Extract the sequences that match our hit list:
for gene in `cat hgt_candidate20170516_list.txt`; do grep "${gene}_" -A1 one_file_to_rule_them_all_tr.gbk.faa >> hgt_candidates_grep.faa ; done

# Run CD_HIT in vital-it, with different identity cutoffs (http://weizhongli-lab.org/lab-wiki/doku.php?id=cd-hit-user-guide)
# this will group together those genes that are from the same family... that is, if they have been "transferred" more than once.
# or those genes that might not be that similar between each other, but belong to a family that has a high transfer frequency!
module load UHTS/Analysis/cd-hit/4.6.1
cd-hit -i hgt_candidates_grep.faa -o clusters60.faa -d 80 -c 0.60 -T 10 -n 4 
cd-hit -i hgt_candidates_grep.faa -o clusters70.faa -d 80 -c 0.70 -T 10 -n 5
cd-hit -i hgt_candidates_grep.faa -o clusters80.faa -d 80 -c 0.80 -T 10 -n 5
cd-hit -i hgt_candidates_grep.faa -o clusters90.faa -d 80 -c 0.90 -T 10 -n 5
cd-hit -i hgt_candidates_grep.faa -o clusters95.faa -d 80 -c 0.95 -T 10 -n 5
cd-hit -i hgt_candidates_grep.faa -o clusters97.faa -d 80 -c 0.97 -T 10 -n 5
cd-hit -i hgt_candidates_grep.faa -o clusters99.faa -d 80 -c 0.99 -T 10 -n 5
# This creates two files per cutoff: one with 

# To count the number of sequences in each cluster:
awk '/^>Cluster/ {if (count) print count; print; count=0; next} {count++} END {print count}' clusters60.faa.clstr | tr '\r\n' ' ' | tr '>' '\n>'

#This perl oneliner is a quick way to explote which seqs are in a specific cluster (267 in this case)
#this PErL script excises all lines between ">Cluster X" and ">Cluster N+1",
#that is, all sequences that are members of that cluster X, from a *.clstr cd-hit output cluster file
#and since the headers are displayed, we can quickly visually explore their annotations.
perl -e 'while (<>){print if (/^>Cluster 267/../^>Cluster 268/);}' clusters60.faa.clstr

# Get all clusters at 60%ID, with format appropiate to split them corresponding to one-file per cluster:
cd-hit -i ../hgt_candidates_grep.faa -o ./clusters60.faa -d 0 -c 0.60 -T 10 -n 4 -g 1
make_multi_seq.pl ../hgt_candidates_grep.faa clusters60.faa.clstr multi-seq 1

# Align each cluster with MUSCLE
module load SequenceAnalysis/MultipleSequenceAlignment/muscle/3.8.31
for i in ./* ; do mv $i ${i}.fasta; done
for i in *.fasta; do muscle -in $i -out ${i%.fasta}.faa; done

# Create one HMMprofile per each cluster
module load SequenceAnalysis/HMM-Profile/hmmer/3.1b2
for i in *.faa; do hmmbuild ${i%.faa}.hmm $i; done

# Create one single database with all Firm5 genomes
cat ../../*.gbk.faa_tr.gbk.faa > ./all_firm5_genomes.faa

# Search every hmmProfile (HGT_putative_cluster) against all genomes, print an alignment of all seqs found, and convert it to FASTA-friendly...
for i in *hmm; do  hmmsearch --tblout ${i%.hmm}.tblout -A ${i%.hmm}.hmmaln --notextw $i ./all_firm5_genomes.faa; esl-reformat a2m ${i%.hmm}.hmmaln ${i%.hmm}.a2m; done

### We need to edit the header names YET AGAIN ###

# Use AnGST with all of them
/scratch/beegfs/monthly/mls_2016/group6/angst



module add Phylogeny/raxml/8.2.9
raxmlHPC -m GTRGAMMA -p 137 -s <seqfile> -n RAxML_phylo

angst: https://github.com/almlab/angst/blob/master/doc/manual.pdf
