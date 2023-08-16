#!/bin/sh -e

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# Document software versions used for publication
hostname
uname -a
fasda --version
pwd

hisat_dir=Results/13-fasda-abundance-hisat2
fasda_dir=Results/14-fasda-fc-hisat2

# Ex. sample16-cond1-rep1
for cond in $(seq 1 5); do
    printf "\n===\nNormalizing cond$cond...\n"
    ls $hisat_dir/sample*-cond$cond-rep*-abundance.tsv
    time fasda normalize \
	--output $fasda_dir/cond$cond-all-norm.tsv \
	$hisat_dir/sample*-cond$cond-rep*-abundance.tsv
done

for first in $(seq 1 4); do
    second_start=$(($first + 1))
    for second in $(seq $second_start 5); do
	printf "\n===\nComputing fold-changes cond$first vs cond$second...\n"
	time fasda fold-change \
	    --output $fasda_dir/cond$first-cond$second-FC.txt \
	    $fasda_dir/cond$first-all-norm.tsv \
	    $fasda_dir/cond$second-all-norm.tsv
    done
done
