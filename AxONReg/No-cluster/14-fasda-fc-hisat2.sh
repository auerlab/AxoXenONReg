#!/bin/sh -e

##########################################################################
#   Script description:
#       Compute fold-changes and P-values for hisat2 data
#
#   Dependencies:
#       Requires hisat2 abundances.  Run after *-fasda-abundance-hisat2.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

date=$(date +%Y-%m-%d-%H:%M)
Sh/14-fasda-fc-hisat2.sh 2>&1 | tee Logs/14-fasda-fc-hisat2/$date.out
