#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on trimmed data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure and trimmed reads.
#       Run after *-organize.sh and 04-trim.sh.
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

hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))
# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/04-trim/*.fastq.zst | \
    xargs -n 1 -P $jobs Xargs/05-qc-trimmed.sh
