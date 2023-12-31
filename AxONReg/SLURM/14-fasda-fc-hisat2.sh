#!/bin/sh -e

##########################################################################
#   Script description:
#       Run fasda on one hisat2 sample
##########################################################################

if which sbatch; then
    sbatch SLURM/14-fasda-fc-hisat2.sbatch
else
    # Debug
    # rm -f Results/14-fasda-fc-hisat2/*
    
    hw_threads=$(./get_hw_threads.sh)
    jobs=$(($hw_threads / 2))
    # Tried GNU parallel and ran into bugs.  Xargs just works.
    ls Results/04-kallisto-quant/*/abundance.tsv \
	| xargs -n 1 -P $jobs Xargs/14-fasda-fc-hisat2.sh
fi
