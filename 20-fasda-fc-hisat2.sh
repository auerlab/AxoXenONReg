#!/bin/sh -e

##########################################################################
#   Script description:
#
#   Dependencies:
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#
#   History:
#   Date        Name        Modification
#   2023-06-02  Jason Bacon Begin
##########################################################################

if which sbatch; then
    sbatch SLURM/20-fasda-fc-hisat2.sbatch
else
    # Debug
    # rm -f Results/20-fasda-fc-hisat2/*
    
    hw_threads=$(./get_hw_threads.sh)
    jobs=$(($hw_threads / 2))
    # Tried GNU parallel and ran into bugs.  Xargs just works.
    ls Results/04-kallisto-quant/*/abundance.tsv \
	| xargs -n 1 -P $jobs Xargs/20-fasda-fc-hisat2.sh
fi
