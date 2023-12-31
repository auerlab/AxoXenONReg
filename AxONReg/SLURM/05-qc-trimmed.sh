#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on trimmed data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure and trimmed reads.
#       Run after *-organize.sh and 04-trim.sh.
##########################################################################

if which sbatch; then
    sbatch SLURM/05-qc-trimmed.sbatch
else
    hw_threads=$(../Common/get-hw-threads.sh)
    jobs=$(($hw_threads / 2))
    # Tried GNU parallel and ran into bugs.  Xargs just works.
    ls Results/04-trim/*.fastq.zst | \
	xargs -n 1 -P $jobs Xargs/05-qc-trimmed.sh
fi
