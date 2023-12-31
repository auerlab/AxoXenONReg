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

if pwd | grep AxoXenOnReg/AxONReg/SLURM; then
    # hisat2 2.2.1
    mem=55g
elif pwd | grep AxoXenOnReg/XenONReg/SLURM; then
    # hisat2 2.2.1
    mem=6g
else
    cat << EOM

This script must be run from AxoXenOnReg/AxONReg/SLURM or from
AxoXenOnReg/XenONReg/SLURM.

EOM
    exit 1
fi

sbatch --mem=$mem SLURM/12-hisat2-align.sbatch
