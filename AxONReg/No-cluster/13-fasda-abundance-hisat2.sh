#!/bin/sh -e

##########################################################################
#   Script description:
#       Compute abundances for each feature in hisat2 BAMs.
#
#   Dependencies:
#       Requires hisat2 BAMs.  Run after *-hisat2-align.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/12-hisat2-align/*.bam | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh \
    Sh/13-fasda-abundance-hisat2.sh
