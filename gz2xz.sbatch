#!/bin/sh -e

#SBATCH --array=1-30
#SBATCH --output=Logs/00-gz2xz/slurm-%A_%a.out
#SBATCH --error=Logs/00-gz2xz/slurm-%A_%a.err

files=$(ls Raw/*/*.fq.gz)
my_file=$(echo $files | awk -v index=$SLURM_ARRAY_TASK_ID '{ print $index }')
set -x
gzcat $my_file | xz > ${my_file%.gz}.xz
