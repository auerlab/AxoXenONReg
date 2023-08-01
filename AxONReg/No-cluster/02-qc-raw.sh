#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure.  Run after 01-organize.sh.
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

# FastQC can utilize 2 cores, so divide total hardware cores by 2
# to determine the number of simultaneous jobs
hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/01-organize/Raw-renamed/*.fastq.xz | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh Sh/02-qc-raw.sh
