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

Sh/07-reference.sh 2>&1 | tee Logs/07-reference/07-reference.out
