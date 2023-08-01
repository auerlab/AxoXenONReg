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

if pwd | grep AxoXenOnReg/AxONReg/Cluster; then
    # STAR 2.7.10b
    mem=100g
elif pwd | grep AxoXenOnReg/XenONReg/Cluster; then
    # STAR 2.7.10b
    # mem=6g
else
    cat << EOM

This script must be run from AxoXenOnReg/AxONReg/Cluster or from
AxoXenOnReg/XenONReg/Cluster.

EOM
    exit 1
fi

sbatch --mem=$mem SLURM/11-star-index.sbatch
