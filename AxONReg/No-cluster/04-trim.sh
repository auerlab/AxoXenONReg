#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters, poly-A tails, and low-quality bases from reads.
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/01-organize/Raw-renamed/*-R1.fastq.xz | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh Sh/04-trim.sh
