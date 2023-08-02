#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#
#   Dependencies:
#       Requires trimmed reads and reference transcriptome.
#       Run after *-trim.sh and *-reference.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

date=$(date +%Y-%m-%d-%H:%M)
Sh/08-kallisto-index.sh 2>&1 | tee Logs/08-kallisto-index/$date.out
