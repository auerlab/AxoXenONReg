#!/bin/sh -e

##########################################################################
#   Script description:
#
#   Dependencies:
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

date=$(date +%Y-%m-%d-%H:%M)
Sh/10-fasda-kallisto.sh 2>&1 | tee Logs/10-fasda-kallisto/$date.out
