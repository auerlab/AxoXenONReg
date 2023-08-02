#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on trimmed data for comparison
#
#   Dependencies:
#       Requires directory structure and trimmed reads.
#       Run after *-organize.sh and 04-trim.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))
# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/04-trim/*.fastq.zst | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh Sh/05-qc-trimmed.sh
