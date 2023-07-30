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

usage()
{
    printf "Usage: $0 axo|xen\n"
    exit 1
}

if [ $# != 1 ]; then
    usage
fi
organism=$1

case $organism in
axo)
    mem=55g
    ;;

xen)
    # FIXME: Guess, verify
    mem=4g
    ;;

*)
    usage
    ;;

esac

sbatch --mem=$mem SLURM/12-hisat2-align.sbatch
