#!/bin/sh -e

##########################################################################
#   Script description:
#       Build kallisto transcriptome index.
#
#   Dependencies:
#       Requires trimmed reads and reference transcriptome.
#       Run after *-trim.sh and *-reference.sh.
##########################################################################

# This script exists only to redirect the output of Sh/08-kallisto-index.sh
date=$(date +%Y-%m-%d-%H:%M)
Sh/08-kallisto-index.sh 2>&1 | tee Logs/08-kallisto-index/$date.out
