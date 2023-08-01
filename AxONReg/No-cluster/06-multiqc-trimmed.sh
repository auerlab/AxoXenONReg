#!/bin/sh -e

##########################################################################
#   Script description:
#       Consolidate FastQC reports into one for easy viewing
#
#   Dependencies:
#       Requires trimmed FastQC results.  Run after *-trim.sh.
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

date=$(date +%Y-%m-%d-%H:%M)
set -x
multiqc --outdir Results/06-multiqc-trimmed Results/05-qc-trimmed \
    2>&1 | tee Logs/06-multiqc-trimmed/$date.out
