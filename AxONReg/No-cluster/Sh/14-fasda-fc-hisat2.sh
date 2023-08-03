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

# Ex. sample16-time1-rep1
for time in $(seq 1 5); do
    printf "\n===\nNormalizing time$time...\n"
    ls $hisat_dir/sample*-time$time-rep*-abundance.tsv
    time fasda normalize \
	--output $fasda_dir/time$time-all-norm.tsv \
	$hisat_dir/sample*-time$time-rep*-abundance.tsv
done

for early in $(seq 1 4); do
    late_start=$(($early + 1))
    for late in $(seq $late_start 5); do
	printf "\n===\nComputing fold-changes time$early vs time$late...\n"
	time fasda fold-change \
	    --output $fasda_dir/time$early-time$late-FC.txt \
	    $fasda_dir/time1-all-norm.tsv \
	    $fasda_dir/time2-all-norm.tsv
    done
done
