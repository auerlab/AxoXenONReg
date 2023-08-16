#!/bin/sh -e

#############################################################################
#   Description:
#       Merge BAM files for the same cell line and all conditions
#
#   Dependencies:
#       Requires hisat2 alignments.  Run after *-hisat2-align.sbatch.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#############################################################################

# Set job array to number of samples / number of conditions.
# 3 time points, chondro and neuro for each
#SBATCH --array=000,026,050,055,120
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
# top showed a steady 25M virtual, 12 res with a spike to 42/30 at the
# end (samtools index?) on FreeBSD 13
#SBATCH --mem=50m
#SBATCH --output=Logs/15-merge-hisat2-bams/slurm-%A_%a.out
#SBATCH --error=Logs/15-merge-hisat2-bams/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=000}

# Document software versions used for publication
hostname
uname -a
samtools --version
pwd

input_dir=Results/12-hisat2-align
output_dir=Results/15-merge-hisat2-bams
cond=$(printf "%03d" $SLURM_ARRAY_TASK_ID)

input1=$(echo $input_dir/sample*-cond$cond-rep1.bam)
input2=$(echo $input_dir/sample*-cond$cond-rep2.bam)
input3=$(echo $input_dir/sample*-cond$cond-rep3.bam)

output=$output_dir/cond$cond.bam

set -x
samtools merge -o $output $input1 $input2 $input3
samtools index $output
