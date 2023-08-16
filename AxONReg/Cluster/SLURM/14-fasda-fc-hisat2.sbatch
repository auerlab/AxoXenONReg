#!/bin/sh -e

#
#   Dependencies:
#       Requires kallisto abundances.  Run after *-kallisto-quant.sbatch.

#SBATCH --output=Logs/14-fasda-fc-hisat2/slurm-%A_%a.out
#SBATCH --error=Logs/14-fasda-fc-hisat2/slurm-%A_%a.err

# Document software versions used for publication
hostname
uname -a
fasda --version
pwd

kallisto_dir=Results/13-fasda-abundance-hisat2
fasda_dir=Results/14-fasda-fc-hisat2

# Ex. sample16-cond1-rep1
for cond in $(seq 1 5); do
    printf "\n===\nNormalizing cond$cond...\n"
    ls $kallisto_dir/sample*-cond$cond-rep*-abundance.tsv
    time fasda normalize \
	--output $fasda_dir/cond$cond-all-norm.tsv \
	$kallisto_dir/sample*-cond$cond-rep*-abundance.tsv
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
