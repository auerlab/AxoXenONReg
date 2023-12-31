#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on trimmed data for comparison
#
#   Dependencies:
#       Requires trimmed reads.  Run after 04-trim.sh.
##########################################################################

hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/04-trim/*.fastq.zst | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh Sh/05-qc-trimmed.sh
