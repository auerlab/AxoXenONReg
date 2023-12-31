#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
##########################################################################

# Set job array to number of samples.
# Each invocation of this script will run 2 fastqc processes, forward
# and reverse.
#SBATCH --array=1-15
#SBATCH --cpus-per-task=2
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=1g
#SBATCH --output=Logs/02-qc-raw/slurm-%A_%a.out
#SBATCH --error=Logs/02-qc-raw/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
hostname
uname -a
fastqc --version
pwd

sample=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
raw1=$(ls Results/01-organize/Raw-renamed/sample$sample-*R1*.fastq.xz)
raw2=$(ls Results/01-organize/Raw-renamed/sample$sample-*R2*.fastq.xz)

# Filename stems for fastqc output
stem_raw1=$(basename ${raw1%.fastq.xz})
stem_raw2=$(basename ${raw2%.fastq.xz})

# Background the first three to run 4 jobs in parallel
xzcat $raw1 | fastqc -o Results/02-qc-raw stdin:$stem_raw1 &
xzcat $raw2 | fastqc -o Results/02-qc-raw stdin:$stem_raw2

# Make sure backgrounded job completes before terminating script
wait
