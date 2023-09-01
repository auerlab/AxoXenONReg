#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
##########################################################################

if pwd | grep AxoXenOnReg/AxONReg/Cluster; then
    # hisat2 2.2.1
    mem=80g
elif pwd | grep AxoXenOnReg/XenONReg/Cluster; then
    # hisat2 2.2.1
    mem=8g
else
    cat << EOM

This script must be run from AxoXenOnReg/AxONReg/Cluster or from
AxoXenOnReg/XenONReg/Cluster.

EOM
    exit 1
fi

sbatch --mem=$mem SLURM/11-hisat2-index.sbatch
