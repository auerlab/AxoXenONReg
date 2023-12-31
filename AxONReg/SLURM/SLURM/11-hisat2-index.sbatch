#!/bin/sh -e

##########################################################################
#   Description:
#       Build hisat2 index for reference genome.
#
#   Dependencies:
#       Requires a reference genome.  Run after *-reference.sbatch.
##########################################################################

# Highest memory use observed was 78 G resident, 120 virtual for axolotl
#SBATCH --output=Logs/11-hisat2-index/slurm-%A.out
#SBATCH --error=Logs/11-hisat2-index/slurm-%A.err
#SBATCH --cpus-per-task=16

if [ -z "$SLURM_MEM_PER_NODE" ]; then
    printf "Usage: $0 --mem=X\n"
    printf "Use ../11-hisat2-index.sh rather than sbatch.\n"
    exit 1
fi

# If not running under SLURM, use all available cores
hw_threads=$(../../Common/get-hw-threads.sh)
: ${SLURM_CPUS_PER_TASK:=$hw_threads}

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
    hisat2-build -p $SLURM_CPUS_PER_TASK $genome $genome
fi
if [ ! -e $genome.fai ]; then
    printf "Building $genome.fai...\n"
    samtools faidx $genome
fi
ls -l
