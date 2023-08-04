#!/bin/sh -e

##########################################################################
#   Script description:
#       Run STAR alignment on all files
#
#   Dependencies:
#       Requires trimmed reads and STAR index.  Run after
#       *-trim.sh and *-star-index.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# Check for sufficient RAM
if pwd | fgrep -q AxoXenOnReg/AxONReg; then
    # FIXME: Only a guess.  Can't build index with 128 GB.
    mem_required=200
else
    # STAR 2.7.10b, measured 27 GB for first process.
    # --genomeLoad should greately reduce requirements for additional jobs.
    mem_required=30
fi

mem=$(../../Common/get-hw-mem.sh)
mem_gb=$(( $mem / 1024 / 1024 / 1024 ))
if pwd | fgrep -q AxoXenOnReg/AxONReg && [ $mem_gb -lt $mem_required ]; then
    printf "You have only $mem_gb GiB RAM and you need $mem_required.\n"
    exit 1
fi

jobs=2
ls Results/04-trim/*-R1.fastq.zst | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh Sh/12-star-align.sh

# Perhaps after configuring shared memory settings on nodes
# STAR --genomeLoad Remove
