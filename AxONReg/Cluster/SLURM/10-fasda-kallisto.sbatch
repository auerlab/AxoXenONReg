#!/bin/sh -e

#
#   Dependencies:
#       Requires kallisto abundances.  Run after *-kallisto-quant.sbatch.

#SBATCH --output=Logs/10-fasda-kallisto/slurm-%A_%a.out
#SBATCH --error=Logs/10-fasda-kallisto/slurm-%A_%a.err

# Document software versions used for publication
uname -a
fasda --version
pwd

kallisto_dir=Results/09-kallisto-quant
fasda_dir=Results/10-fasda-kallisto

# Ex. sample16-time1-rep1
for time in $(seq 1 5); do
    printf "\n===\nNormalizing time$time...\n"
    ls $kallisto_dir/sample*-time$time-rep*/abundance.tsv
    time fasda normalize \
	--output $fasda_dir/time$time-all-norm.tsv \
	$kallisto_dir/sample*-time$time-rep*/abundance.tsv
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
