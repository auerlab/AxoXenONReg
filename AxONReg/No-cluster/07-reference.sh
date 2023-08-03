#!/bin/sh -e

##########################################################################
#   Description:
#       Download/build genome and transcriptome references.
#       
#   Prerequisites:
#       None.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# This script exists only to redirect the output of Sh/07-reference.sh
date=$(date +%Y-%m-%d-%H:%M)
Sh/07-reference.sh 2>&1 | tee Logs/07-reference/$date.out
