#!/bin/sh -e

##########################################################################
#   Script description:
#       Consolidate FastQC reports into one for easy viewing
#
#   Dependencies:
#       Requires raw FastQC results.  Run after *-qc-raw.sbatch.
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

##########################################################################
# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
# export LC_ALL=en_US.UTF8

LC_ALL=en_US.UTF-8
LANG=en_US.utf-8
export LC_ALL LANG

set -x
multiqc --outdir Results/03-multiqc-raw Results/02-qc-raw \
    2>&1 | tee Logs/03-multiqc-raw/multiqc.out
