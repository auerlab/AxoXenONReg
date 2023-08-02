#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# Debug
# rm -f Results/17-hisat2-align/*

hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))
# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/12-hisat2-align/*.bam | \
    xargs -n 1 -P $jobs Sh/13-fasda-abundance-hisat2.sh
