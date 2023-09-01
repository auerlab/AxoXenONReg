#!/bin/sh -e

##########################################################################
#   Description:
#       Compute fold-changes and P-values for kallisto abundances.
##########################################################################

# Document software versions used for publication
hostname
uname -a
fasda --version
pwd

kallisto_dir=Results/09-kallisto-quant
fasda_dir=Results/10-fasda-kallisto

# Ex. sample16-cond1-rep1
for cond in $(seq 1 5); do
    printf "\n*** Normalizing cond$cond...\n\n"
    ls $kallisto_dir/sample*-cond$cond-rep*/abundance.tsv
    time fasda normalize \
	--output $fasda_dir/cond$cond-all-norm.tsv \
	$kallisto_dir/sample*-cond$cond-rep*/abundance.tsv
done

for first in $(seq 1 4); do
    second_start=$(($first + 1))
    for second in $(seq $second_start 5); do
	printf "\n*** Computing fold-changes cond$first vs cond$second...\n"
	time fasda fold-change \
	    --output $fasda_dir/cond$first-cond$second-FC.txt \
	    $fasda_dir/cond$first-all-norm.tsv \
	    $fasda_dir/cond$second-all-norm.tsv
    done
done
