#!/bin/sh -e

##########################################################################
#   Script description:
#       Generate hisat2 genome reference index.
#
#   Dependencies:
#       Requires genome reference.  Run after *-reference.sh.
#       For axolotl, requires nearly 80 GB RAM.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

Sh/11-hisat2-index.sh \
    2>&1 | tee Logs/11-hisat2-index/hisat2-index.out
