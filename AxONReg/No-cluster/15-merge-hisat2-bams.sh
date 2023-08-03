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

hw_threads=$(../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/01-organize/Raw-renamed/*-R1.fastq.xz | \
    xargs -n 1 -P $jobs Sh/15-merge-hisat2-bams.sh

