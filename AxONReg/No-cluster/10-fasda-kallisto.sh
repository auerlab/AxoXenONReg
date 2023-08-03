#!/bin/sh -e

##########################################################################
#   Script description:
#       Compute fold-changes and P-values for kallisto abundances.
#
#   Dependencies:
#       Requires kallisto abundances.  Run after *-kallisto-align.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

date=$(date +%Y-%m-%d-%H:%M)
Sh/10-fasda-kallisto.sh 2>&1 | tee Logs/10-fasda-kallisto/$date.out
