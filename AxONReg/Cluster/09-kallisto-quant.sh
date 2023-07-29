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

    # Kallisto multithreading scales fairly well, so more jobs with
    # fewer cores each won't lead to much better performance
    # Using minimum number of jobs with at least 4 GB each to ensure
    # plenty of RAM per job
    hw_threads=$(../Common/get-hw-threads.sh)
    hw_mem=$(../Common/get-hw-mem.sh)
    
    # Find optimal # cores per job with at least 4 GB per job
    jobs=$hw_threads
    mem_per_job=$(($hw_mem / $jobs ))
    while [ $mem_per_job -lt 4000000000 ]; do
	jobs=$(( $jobs / 2 ))
	mem_per_job=$(($hw_mem / $jobs ))
    done
    threads_per_job=$(( $hw_threads / $jobs ))
    
    # Tried GNU parallel and ran into bugs.  Xargs just works.
    printf "Running $jobs jobs with $threads_per_job threads each...\n"
    ls Results/04-trim/*-R1.fastq.zst | \
	xargs -n 1 -P $jobs Xargs/09-kallisto-quant.sh $threads_per_job
fi
