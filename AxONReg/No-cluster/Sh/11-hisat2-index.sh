#!/bin/sh -e

##########################################################################
#   Description:
#       Build hisat2 index for reference genome.
#
#   Dependencies:
#       Requires a reference genome.  Run after *-reference.sbatch.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# FIXME: Check sufficient RAM (80 GB)

# If not running under SLURM, use all available cores
hw_threads=$(../../Common/get-hw-threads.sh)

# Document software versions used for publication
hostname
uname -a
hisat2 --version
samtools --version
pwd

# Run hisat2-build on a copy in 11-hisat2-index so it will put the .ht2
# files there
genome=genome-reference.fa
ln -f Results/07-reference/$genome Results/11-hisat2-index
printf "Using reference $genome...\n"

# FIXME: Do this in *-reference.sbatch
if [ ! -e $genome.8.ht2 ]; then
    printf "Building $genome.*.ht2...\n"
    # Unlike hisat2, hisat2-build does not support --threads
    # as a -p equivalent
    cd Results/11-hisat2-index
    set -x
    hisat2-build -p $hw_threads $genome $genome
fi
if [ ! -e $genome.fai ]; then
    printf "Building $genome.fai...\n"
    samtools faidx $genome
fi
ls -l
