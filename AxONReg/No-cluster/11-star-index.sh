#!/bin/sh -e

##########################################################################
#   Script description:
#       Generate star genome reference index.
#
#   Dependencies:
#       Requires genome reference.  Run after *-reference.sh.
#       For axolotl, requires more than 128 GB RAM (actual amount unknown).
#       For xenopus, requires nearly 60 GB RAM.
##########################################################################

# Check for sufficient RAM
if pwd | fgrep -q AxoXenOnReg/AxONReg; then
    # FIXME: Only a guess.  128 is not enough.
    mem_required=200
else
    mem_required=60
fi

mem=$(../../Common/get-hw-mem.sh)
mem_gb=$(( $mem / 1024 / 1024 / 1024 ))
if pwd | fgrep -q AxoXenOnReg/AxONReg && [ $mem_gb -lt $mem_required ]; then
    printf "You have only $mem_gb GiB RAM and you need $mem_required.\n"
    exit 1
fi

# This script exists only to redirect the output of Sh/11-star-index.sh.
date=$(date +%Y-%m-%d-%H:%M)
Sh/11-star-index.sh 2>&1 | tee Logs/11-star-index/$date.out
