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
#   2023-06-02  Jason Bacon Begin
##########################################################################

if which sbatch; then
    sbatch SLURM/08-kallisto-index.sbatch
else
    # Debug
    # rm -f Results/08-kallisto-index/*
    
    hw_threads=$(./get_hw_threads.sh)
    jobs=$(($hw_threads / 2))
    # Tried GNU parallel and ran into bugs.  Xargs just works.
    ls Results/01-organize/Raw-renamed/*-R1.fastq.xz | \
	xargs -n 1 -P $jobs Xargs/08-kallisto-index.sh
fi
