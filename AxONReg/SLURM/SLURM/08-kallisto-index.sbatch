#!/bin/sh -e

##########################################################################
#   Description:
#       Build kallisto index for reference transcriptome.
#
#   Dependencies:
#       Requires reference transriptome.  Run after *-reference.sbatch.
##########################################################################

# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=4g
#SBATCH --cpus-per-task=16
#SBATCH --output=Logs/08-kallisto-index/slurm-%A.out
#SBATCH --error=Logs/08-kallisto-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
hostname
uname -a
kallisto version
samtools --version
pwd

ref_dir=Results/07-reference
transcriptome=transcriptome-reference.fa
index=${transcriptome%.fa}.index

printf "Using reference $transcriptome...\n"

# Needed for kallisto --genomebam
if [ ! -e $ref_dir/$transcriptome.fai ]; then
    printf "Building samtools index...\n"
    samtools faidx $ref_dir/$transcriptome
fi

printf "Building kallisto index...\n"
set -x
kallisto index --threads=$SLURM_CPUS_PER_TASK \
    --index=Results/08-kallisto-index/$index $ref_dir/$transcriptome
