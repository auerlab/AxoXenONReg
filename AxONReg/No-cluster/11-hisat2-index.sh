#!/bin/sh -e

##########################################################################
#   Script description:
#       Generate hisat2 genome reference index.
#
#   Dependencies:
#       Requires genome reference.  Run after *-reference.sh.
#       For axolotl, requires nearly 80 GB RAM.
##########################################################################

# Check for sufficient RAM (80 GB)
mem=$(../../Common/get-hw-mem.sh)
mem_gb=$(( $mem / 1024 / 1024 / 1024 ))
if pwd | fgrep -q AxoXenOnReg/AxONReg && [ $mem_gb -lt 80 ]; then
    printf "You have only $mem_gb GiB RAM and you need 80 for axolotl.\n"
    exit 1
fi

# This script exists only to redirect the output of Sh/11-hisat2-index.sh.
date=$(date +%Y-%m-%d-%H:%M)
Sh/11-hisat2-index.sh 2>&1 | tee Logs/11-hisat2-index/$date.out
