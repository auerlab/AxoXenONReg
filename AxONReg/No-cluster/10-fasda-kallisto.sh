#!/bin/sh -e

##########################################################################
#   Script description:
#       Compute fold-changes and P-values for kallisto abundances.
#
#   Dependencies:
#       Requires kallisto abundances.  Run after *-kallisto-align.sh.
##########################################################################

# This script exists only to redirect the output of Sh/10-fasda-abundance.sh.
date=$(date +%Y-%m-%d-%H:%M)
Sh/10-fasda-kallisto.sh 2>&1 | tee Logs/10-fasda-kallisto/$date.out
