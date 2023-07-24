#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

if which sbatch; then
    sbatch SLURM/09-kallisto-quant.sbatch
else
    # Debug
    # rm -f Results/09-kallisto-quant/*
    
    # Consider both CPU cores and memory when selecting thread count
    # Use all cores for one job here to minimize contention and maximize
    # throughput
    threads_per_job=$(../Common/get-hw-threads.sh)
    hw_threads=$(../Common/get-hw-threads.sh)
    jobs=$(($hw_threads / $threads_per_job))
    # Tried GNU parallel and ran into bugs.  Xargs just works.
    ls Results/04-trim/*-R1.fastq.zst | \
	xargs -n 1 -P $jobs Xargs/09-kallisto-quant.sh $threads_per_job
fi
